#!/usr/bin/env python3

# Written by Jasmin Frangenberg and released under the MIT license.
# See below for full license text.

from Bio import SeqIO
import pandas as pd
import argparse
import os
import re
from concurrent.futures import ProcessPoolExecutor
import warnings
from functions.filter import filter_bgc, cleanup_table, combine_tool, process_chunk, parallelization
from functions.deepbgc_workflow import deepbgc_workflow
from functions.gecco_workflow import gecco_workflow, getInterProID
from functions.antismash_workflow import antismash_workflow, prepare_multisample_input_antismash, parse_knownclusterblast, antismash_workflow

"""
===============================================================================
MIT License
===============================================================================

Copyright (c) 2023 Jasmin Frangenberg

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""


#########################################
# FUNCTION: Add sample metadata
#########################################
def sample_metadata_addition(merged_df, smetadata):
    """
    Adds the sample metadata only when turned on.
    Important to note: the first column must have the sample name and the should be tab seperated.
    """
    if smetadata == None :
        return merged_df
    else:
        # add the samples metadata 
        metadata_df = pd.read_csv(smetadata, sep='\t')
        metadata_df.rename(columns={metadata_df.columns[0]: 'sample_id'}, inplace=True)
        # merge it to the df using sample name 'name' as common
        df1 = merged_df.merge(metadata_df, on='sample_id', how='left')
        return df1



#########################################
# FUNCTION: Add contig metadata
#########################################
def contig_metadata_addition(merged_df, cmetadata):
    """
    Adds the contig metadata only when turned on.
    Important to note: the first column must have the sample name and the second name must have contig_id
                       the table should be tab seperated.
    """
    if cmetadata == None :
        return merged_df
    else:
        # add the samples metadata 
        merged_df = merged_df.drop("mmseqs_lineage_contig", axis=1)
        metadata_df = pd.read_csv(cmetadata, sep='\t')
        # rename first column : must contain the sample names
        metadata_df.rename(columns={metadata_df.columns[0]: 'sample_id'}, inplace=True)
        # rename second column _ must contain the contig names
        metadata_df.rename(columns={metadata_df.columns[1]: 'contig_id'}, inplace=True)
        # ensure that the dataframe is not empty or else KeyError arises
        if not merged_df.empty:
            # remove duplicate entries
            merged_df.drop_duplicates(inplace=True)
            # merge it to the df using sample name 'name' as common
            df2 = pd.merge(merged_df, metadata_df, on=['contig_id', 'sample_id'], how='left')
        else:
            df2 = merged_df
        return df2





########################
# MAIN
########################

def main():
    warnings.filterwarnings("ignore", category=FutureWarning, message="The behavior of DataFrame concatenation with empty or all-NA entries is deprecated.*")

    tool_version = "0.6.8"
    welcome = """\
                                                        
        ██████╗ ██████╗ ███╗   ███╗██████╗  ██████╗  ██████╗
        ██╔════╝██╔═══██╗████╗ ████║██╔══██╗██╔════╝ ██╔════╝
        ██║     ██║   ██║██╔████╔██║██████╔╝██║  ███╗██║     
        ██║     ██║   ██║██║╚██╔╝██║██╔══██╗██║   ██║██║     
        ╚██████╗╚██████╔╝██║ ╚═╝ ██║██████╔╝╚██████╔╝╚██████╗
         ╚═════╝ ╚═════╝ ╚═╝     ╚═╝╚═════╝  ╚═════╝  ╚═════╝
                                                        
                    ........................
                        * comBGC v.{version} *
                    ........................
        This tool aggregates the results of BGC prediction tools:
                    antiSMASH, deepBGC, and GECCO
             For detailed usage documentation please refer to 
                      https://nf-co.re/funcscan
        .........................................................""".format(
        version=tool_version
    )

    # Initialize parser
    parser = argparse.ArgumentParser(
        prog="comBGC",
        formatter_class=argparse.RawTextHelpFormatter,
        description=(welcome),
        add_help=True,
    )

    # Input options
    parser.add_argument(
        "-i",
        "--input",
        metavar="PATH(s)",
        dest="input",
        nargs="*",
        help="""path(s) to the required output file(s) of antiSMASH, DeepBGC and/or GECCO
    these can be:
    - antiSMASH: <sample name>.gbk and (optional) knownclusterblast/ directory
    - DeepBGC:   <sample name>.bgc.tsv
    - GECCO:     <sample name>.clusters.tsv
Note: Please provide files from a single sample only. If you would like to
summarize multiple samples, please see the --antismash_multiple_samples flag.""",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        metavar="PATH",
        dest="outdir",
        nargs="?",
        help="directory for comBGC output. Default: current directory",
        type=str,
        default=".",
    )
    parser.add_argument(
        "-a",
        "--antismash_multiple_samples",
        metavar="PATH",
        dest="antismash_multiple_samples",
        nargs="?",
        help="""directory of antiSMASH output. Should contain subfolders (one per
sample). Can only be used if --input is not specified.""",
        type=str,
    )
    parser.add_argument(
        "-vv", "--verbose", help="increase output verbosity", action="store_true"
    )
    parser.add_argument(
        "-v", "--version", help="show version number and exit", action="store_true"
    )
    parser.add_argument(
        "--cores", type=int, help="number of CPU cores to use"
    )
    parser.add_argument(
        "--min_length",
        type=int,
        metavar="LENGTH",
        dest="min_length",
        default=3000,
        help="Minimum length [bp] of BGC sequence to keep. Default: 3000"
    )
    parser.add_argument(
        "--contig_edge",
        type=int,
        metavar="BASES",
        dest="contig_edge",
        default=2,
        help="Exclude BGCs within X bases of the contig edge. Default: 2"
    )

    parser.add_argument(
        "--sample_metadata", 
        dest="samplemetadata", 
        metavar="[PATH]", 
        help="""Path to a TSV file containing sample metadata, e.g., 'path/to/sample_metadata.tsv'. 
The metadata table must have sample names in the first column.""",
        type=str, 
        default=None
    )

    parser.add_argument(
        "--contig_metadata", 
        dest="contigmetadata", 
        metavar="[PATH]",
        help="""Path to a TSV file containing contig metadata, e.g., 'path/to/contig_metadata.tsv'. 
The metadata table must have sample names in the first column and contig IDs in the second column.""",
        type=str, 
        default=None
    )

    # Get command line arguments
    args = parser.parse_args()

    # Assign input arguments to variables
    input = args.input
    dir_antismash = args.antismash_multiple_samples
    outdir = args.outdir
    verbose = args.verbose
    version = args.version
    cores = args.cores
    min_length = args.min_length
    contig_edge = args.contig_edge
    add_samplemetadata = args.samplemetadata
    add_contigmetadata = args.contigmetadata

    if version:
        exit("comBGC {version}".format(version=tool_version))
    if cores is None:
        parser.error("--cores is required. Please specify the number of CPU cores to use.")

    input_antismash = []
    input_deepbgc = []
    input_gecco = []

    # Assign input files to respective tools
    if input:
        for path in input:
            if path.endswith(".gbk") and not re.search(r"region\d\d\d\.gbk$", path): # Make sure to only fetch relevant GBK files, i.e. those containing all collective antiSMASH BGCs
                with open(path) as infile:
                    for line in infile:
                        if re.search("##GECCO-Data-START##", line):
                            input_gecco.append(path)
                            break
                        elif re.search("##antiSMASH-Data-START##", line):
                            input_antismash.append(path)
                            break
            elif path.endswith("bgc.tsv"):
                input_deepbgc = path
            elif path.endswith("clusters.tsv"):
                input_gecco.append(path)
            elif path.endswith("knownclusterblast/"):
                input_antismash.append(path)

    if input and dir_antismash:
        exit(
            "The flags --input and --antismash_multiple_samples are mutually exclusive.\nPlease use only one of them (or see --help for how to use)."
        )

    # Make sure that at least one input argument is given
    if not (input_antismash or input_gecco or input_deepbgc or dir_antismash):
        exit(
            "Please specify at least one input file (i.e. output from antismash, deepbgc, or gecco) or see --help"
        )

    if input_antismash:
        tools = {
            "antiSMASH": input_antismash,
            "deepBGC": input_deepbgc,
            "GECCO": input_gecco,
        }
    elif dir_antismash:
        tools = {"antiSMASH": dir_antismash}
    else:
        tools = {"deepBGC": input_deepbgc, "GECCO": input_gecco}

    tools_provided = {}

    for tool in tools.keys():
        if tools[tool]:
            tools_provided[tool] = tools[tool]

    if verbose:
        print(welcome)
        print("\nYou provided input for: " + ", ".join(tools_provided.keys()))

    # Aggregate BGC information into data frame
    summary_antismash = pd.DataFrame()
    summary_deepbgc = pd.DataFrame()
    summary_gecco = pd.DataFrame()

    for tool in tools_provided.keys():
        if tool == "antiSMASH":
            if dir_antismash:
                antismash_paths = prepare_multisample_input_antismash(dir_antismash)
                for input_antismash in antismash_paths:
                    summary_antismash_temp = antismash_workflow(input_antismash, verbose)
                    summary_antismash = pd.concat(
                        [summary_antismash, summary_antismash_temp]
                    )
            else:
                summary_antismash = antismash_workflow(input_antismash, verbose)
        elif tool == "deepBGC":
            summary_deepbgc = deepbgc_workflow(input_deepbgc, verbose)
        elif tool == "GECCO":
            summary_gecco = gecco_workflow(input_gecco, verbose)

    # Summarize and sort data frame
    summary_all = pd.concat([summary_antismash, summary_deepbgc, summary_gecco])
    summary_all.sort_values(
        by=["Sample_ID", "Contig_ID", "BGC_start", "BGC_length", "Prediction_tool"],
        axis=0,
        inplace=True,
    )

    # Rearrange and rename the columns in the summary df
    summary_all = summary_all.iloc[:, [0, 2, 1] + list(range(3, len(summary_all.columns)))]
    summary_all.rename(columns={'Sample_ID':'sample_id', 'Contig_ID':'contig_id', 'CDS_ID':'BGC_region_contig_ids'}, inplace=True)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # filter and standardize the summary table
    df = summary_all.copy()
    df.columns = df.columns.str.strip()

    df["BGC_length"] = pd.to_numeric(df["BGC_length"], errors="coerce", downcast="integer")
    df["BGC_start"] = pd.to_numeric(df["BGC_start"], errors="coerce", downcast="integer")
    df["BGC_end"] = pd.to_numeric(df["BGC_end"], errors="coerce", downcast="integer")
    df["BGC_probability"] = pd.to_numeric(df["BGC_probability"], errors="coerce", downcast="integer")

    df = filter_bgc(df, min_length, contig_edge)
    if df.empty:
        raise ValueError("No BGCs remain after filtering.")    
    df = cleanup_table(df)
    
    # dereplicated BGCs with the same identifier
    filtered_bgcs = parallelization(df, cores, verbose)
    filtered_bgcs["Product_class"] = filtered_bgcs["Product_class"].replace("", "Unknown")  # Replace empty strings
    filtered_bgcs["Product_class"] = filtered_bgcs["Product_class"].fillna("Unknown")
    filtered_bgcs = filtered_bgcs.drop("identifier", axis=1)
    
    # merge additional metadata if present
    sample_metadata_df = sample_metadata_addition(filtered_bgcs, add_samplemetadata)
    filtered_bgcs = sample_metadata_df
    contig_metadata_df = contig_metadata_addition(filtered_bgcs, add_contigmetadata)
    filtered_bgcs = contig_metadata_df


    # Write results to TSV
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    filtered_bgcs.to_csv(
        os.path.join(outdir, "combgc_summary.tsv"), sep="\t", index=False
    )
    print("Your BGC summary file is: " + os.path.join(outdir, "combgc_summary.tsv"))

if __name__ == "__main__":
    main()