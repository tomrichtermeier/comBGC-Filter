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

tool_version = "0.6.3"
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
    For detailed usage documentation please refer
    to https://nf-co.re/funcscan
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
    dest="[PATH]", 
    help="""Path to a TSV file containing sample metadata, e.g., 'path/to/sample_metadata.tsv'. 
The metadata table must have sample names in the first column.""",
    type=str, 
    default=None
)

parser.add_argument(
    "--contig_metadata", 
    dest="[PATH]", 
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

########################
# ANTISMASH FUNCTIONS
########################


def prepare_multisample_input_antismash(antismash_dir):
    """
    Prepare string of input paths of a given antiSMASH output folder (with sample subdirectories)
    """
    sample_paths = []
    for root, subdirs, files in os.walk(antismash_dir):
        antismash_file = "/".join([root, "index.html"])
        if os.path.exists(antismash_file):
            sample = root.split("/")[-1]
            gbk_path = "/".join([root, sample]) + ".gbk"
            kkb_path = "/".join([root, "knownclusterblast"])
            if os.path.exists(kkb_path):
                sample_paths.append([gbk_path, kkb_path])
            else:
                sample_paths.append([gbk_path])
    return sample_paths


def parse_knownclusterblast(kcb_file_path):
    """
    Extract MIBiG IDs from knownclusterblast TXT file.
    """

    with open(kcb_file_path) as kcb_file:
        hits = 0
        MIBiG_IDs = []

        for line in kcb_file:
            if line == "Significant hits: \n" and not hits:
                hits = 1  # Indicating that the following lines contain relevant information
            elif line == "\n" and hits:
                break
            elif line != "Significant hits: \n" and hits:
                MIBiG_ID = re.search(r"(BGC\d+)", line).group(1)
                MIBiG_IDs.append(MIBiG_ID)
    return MIBiG_IDs


def antismash_workflow(antismash_paths):
    """
    Create data frame with aggregated antiSMASH output:
    - Open summary GBK and grab relevant information.
    - Extract the knownclusterblast output from the antiSMASH folder (MIBiG annotations) if present.
    - Return data frame with aggregated info.
    """

    antismash_sum_cols = [
        "Sample_ID",
        "Prediction_tool",
        "Contig_ID",
        "Product_class",
        "BGC_probability",
        "BGC_complete",
        "BGC_start",
        "BGC_end",
        "BGC_length",
        "CDS_ID",
        "CDS_count",
        "PFAM_domains",
        "MIBiG_ID",
        "InterPro_ID",
    ]
    antismash_out = pd.DataFrame(columns=antismash_sum_cols)

    CDS_ID = []
    CDS_count = 0

    # Distinguish input files (i.e. GBK file and "knownclusterblast" folder)
    kcb_path = []
    for path in antismash_paths:
        if re.search("knownclusterblast", path):
            kcb_path = re.search(".*knownclusterblast.*", path).group()
        else:
            gbk_path = path

        kcb_files = []
        if kcb_path:
            kcb_files = [
                file
                for file in os.listdir(kcb_path)
                if file.startswith("c") and file.endswith(".txt")
            ]

        # Aggregate information
        Sample_ID = gbk_path.split("/")[-1].split(".gbk")[
            -2
        ]  # Assuming file name equals sample name
        if verbose:
            print("\nParsing antiSMASH file(s): " + Sample_ID + "\n... ", end="")

        with open(gbk_path) as gbk:
            for record in SeqIO.parse(
                gbk, "genbank"
            ):  # GBK records are contigs in this case
                # Initiate variables per contig
                cluster_num = 1
                antismash_out_line = {}
                Contig_ID = record.id
                Product_class = ""
                BGC_complete = ""
                BGC_start = ""
                BGC_end = ""
                BGC_length = ""
                PFAM_domains = []
                MIBiG_ID = "NA"

                for feature in record.features:
                    # Extract relevant infos from the first protocluster feature from the contig record
                    if feature.type == "protocluster":
                        if (
                            antismash_out_line
                        ):  # If there is more than 1 BGC per contig, reset the output line for new BGC. Assuming that BGCs do not overlap.
                            if not CDS_ID:
                                CDS_ID = ["NA"]
                            antismash_out_line = {  # Create dictionary of BGC info
                                "Sample_ID": Sample_ID,
                                "Prediction_tool": "antiSMASH",
                                "Contig_ID": Contig_ID,
                                "Product_class": ";".join(Product_class),
                                "BGC_probability": "NA",
                                "BGC_complete": BGC_complete,
                                "BGC_start": BGC_start,
                                "BGC_end": BGC_end,
                                "BGC_length": BGC_length,
                                "CDS_ID": ";".join(CDS_ID),
                                "CDS_count": CDS_count,
                                "PFAM_domains": ";".join(PFAM_domains),
                                "MIBiG_ID": MIBiG_ID,
                                "InterPro_ID": "NA",
                            }
                            antismash_out_line = pd.DataFrame([antismash_out_line])
                            antismash_out = pd.concat(
                                [antismash_out, antismash_out_line], ignore_index=True
                            )
                            antismash_out_line = {}

                            # Reset variables per BGC
                            CDS_ID = []
                            CDS_count = 0
                            PFAM_domains = []

                        # Extract all the BGC info
                        Product_class = feature.qualifiers["product"]
                        for i in range(len(Product_class)):
                            Product_class[i] = (
                                Product_class[i][0].upper() + Product_class[i][1:]
                            )  # Make first letters uppercase, e.g. lassopeptide -> Lassopeptide

                        if feature.qualifiers["contig_edge"] == ["True"]:
                            BGC_complete = "No"
                        elif feature.qualifiers["contig_edge"] == ["False"]:
                            BGC_complete = "Yes"

                        BGC_start = (
                            feature.location.start + 1
                        )  # +1 because zero-based start position
                        BGC_end = feature.location.end
                        BGC_length = feature.location.end - feature.location.start

                        # If there are knownclusterblast files for the BGC, get MIBiG IDs of their homologs
                        if kcb_files:
                            print(kcb_files)
                            kcb_file = "{}_c{}.txt".format(
                                record.id, str(cluster_num)
                            )  # Check if this filename is among the knownclusterblast files
                            if kcb_file in kcb_files:
                                MIBiG_IDs = ";".join(
                                    parse_knownclusterblast(
                                        os.path.join(kcb_path, kcb_file)
                                    )
                                )
                                if MIBiG_IDs != "":
                                    MIBiG_ID = MIBiG_IDs
                                cluster_num += 1

                    # Count functional CDSs (no pseudogenes) and get the PFAM annotation
                    elif (
                        feature.type == "CDS"
                        and "translation" in feature.qualifiers.keys()
                        and BGC_start != ""
                    ):  # Make sure not to count pseudogenes (which would have no "translation tag") and count no CDSs before first BGC
                        if (
                            feature.location.end <= BGC_end
                        ):  # Make sure CDS is within the current BGC region
                            if "locus_tag" in feature.qualifiers:
                                CDS_ID.append(feature.qualifiers["locus_tag"][0])
                            CDS_count += 1
                            if "sec_met_domain" in feature.qualifiers.keys():
                                for PFAM_domain in feature.qualifiers["sec_met_domain"]:
                                    PFAM_domain_name = re.search(
                                        r"(.+) \(E-value", PFAM_domain
                                    ).group(1)
                                    PFAM_domains.append(PFAM_domain_name)

                # Create dictionary of BGC info
                if not CDS_ID:
                    CDS_ID = ["NA"]
                antismash_out_line = {
                    "Sample_ID": Sample_ID,
                    "Prediction_tool": "antiSMASH",
                    "Contig_ID": Contig_ID,
                    "Product_class": ";".join(Product_class),
                    "BGC_probability": "NA",
                    "BGC_complete": BGC_complete,
                    "BGC_start": BGC_start,
                    "BGC_end": BGC_end,
                    "BGC_length": BGC_length,
                    "CDS_ID": ";".join(CDS_ID),
                    "CDS_count": CDS_count,
                    "PFAM_domains": ";".join(PFAM_domains),
                    "MIBiG_ID": MIBiG_ID,
                    "InterPro_ID": "NA",
                }

                if BGC_start != "":  # Only keep records with BGCs
                    antismash_out_line = pd.DataFrame([antismash_out_line])
                    antismash_out = pd.concat(
                        [antismash_out, antismash_out_line], ignore_index=True
                    )

                    # Reset variables per BGC
                    CDS_ID = []
                    CDS_count = 0
                    PFAM_domains = []

    if verbose:
        print("Done.")
    return antismash_out


########################
# DEEPBGC FUNCTIONS
########################


def deepbgc_workflow(deepbgc_path):
    """
    Create data frame with aggregated deepBGC output.
    """

    if verbose:
        print("\nParsing deepBGC file\n... ", end="")

    # Prepare input and output columns
    deepbgc_map_dict = {
        "sequence_id": "Contig_ID",
        "nucl_start": "BGC_start",
        "nucl_end": "BGC_end",
        "nucl_length": "BGC_length",
        "num_proteins": "CDS_count",
        "deepbgc_score": "BGC_probability",
        "product_class": "Product_class",
        "protein_ids": "CDS_ID",
        "pfam_ids": "PFAM_domains",
    }
    deepbgc_sum_cols = [
        "Sample_ID",
        "Prediction_tool",
        "Contig_ID",
        "Product_class",
        "BGC_probability",
        "BGC_complete",
        "BGC_start",
        "BGC_end",
        "BGC_length",
        "CDS_ID",
        "CDS_count",
        "PFAM_domains",
        "MIBiG_ID",
        "InterPro_ID",
    ]
    deepbgc_unused_cols = [
        "detector_version",
        "detector_label",
        "bgc_candidate_id",
        "num_domains",
        "num_bio_domains",
        "product_activity",
        "antibacterial",
        "cytotoxic",
        "inhibitor",
        "antifungal",
        "Alkaloid",
        "NRP",
        "Other",
        "Polyketide",
        "RiPP",
        "Saccharide",
        "Terpene",
        "bio_pfam_ids",
    ]

    # Grab deepBGC sample ID
    sample = os.path.basename(deepbgc_path).rsplit(".bgc", 1)[0]

    # Initiate dataframe
    deepbgc_out = pd.DataFrame(columns=deepbgc_sum_cols)

    # Add relevant deepBGC output columns per BGC
    deepbgc_df = (
        pd.read_csv(deepbgc_path, sep="\t")
        .drop(deepbgc_unused_cols, axis=1)
        .rename(columns=deepbgc_map_dict)
    )
    deepbgc_df["Sample_ID"] = sample
    deepbgc_df["Prediction_tool"] = "deepBGC"
    deepbgc_df["BGC_complete"] = "NA"
    deepbgc_df["MIBiG_ID"] = "NA"
    deepbgc_df["InterPro_ID"] = "NA"

    # Concatenate data frame to out w/o common index column (e.g. sample_id) due to duplicate row names
    if not deepbgc_df.empty:
      deepbgc_out = pd.concat([deepbgc_out, deepbgc_df], ignore_index=True, sort=False)

    # Return data frame with ordered columns
    deepbgc_out = deepbgc_out[deepbgc_sum_cols]
    if verbose:
        print("Done.")
    return deepbgc_out


########################
# GECCO FUNCTIONS
########################


def getInterProID(gbk_path):
    """
    Retrieve InterPro IDs from GECCO GBK file.
    """

    with open(gbk_path) as gbk:
        ip_ids = []
        id_pattern = r'InterPro\:(.*)"'

        for line in gbk:
            if line.find(r"InterPro:") != -1:
                new_id = re.search(id_pattern, line).group(1)
                ip_ids.append(new_id)
        ipids_str = ";".join(map(str, ip_ids))
    return ipids_str


def gecco_workflow(gecco_paths):
    """
    Create data frame with aggregated GECCO output.
    """

    if verbose:
        print("\nParsing GECCO files\n... ", end="")

    # GECCO output columns that can be mapped (comBGC:GECCO)
    map_dict = {
        "sequence_id": "Contig_ID",
        "bgc_id": "cluster_id",
        "type": "Product_class",
        "average_p": "BGC_probability",
        "start": "BGC_start",
        "end": "BGC_end",
        "domains": "PFAM_domains",
        "proteins": "CDS_ID",
    }
    summary_cols = [
        "Sample_ID",
        "Prediction_tool",
        "Contig_ID",
        "Product_class",
        "BGC_probability",
        "BGC_complete",
        "BGC_start",
        "BGC_end",
        "BGC_length",
        "CDS_ID",
        "CDS_count",
        "PFAM_domains",
        "MIBiG_ID",
        "InterPro_ID",
    ]
    unused_cols = [
        "max_p",
        "alkaloid_probability",
        "polyketide_probability",
        "ripp_probability",
        "saccharide_probability",
        "terpene_probability",
        "nrp_probability",
    ]

    tsv_path = ""
    gbk_paths = []

    for path in gecco_paths:
        if path.endswith(".tsv"):
            tsv_path = path
        else:
            gbk_paths.append(path)

    # Initiate dataframe
    gecco_out = pd.DataFrame(columns=summary_cols)

    # Add sample information
    sample = tsv_path.split("/")[-1].split(".")[0]
    gecco_df = (
        pd.read_csv(tsv_path, sep="\t")
        .drop(unused_cols, axis=1)
        .rename(columns=map_dict)
    )

    # Fill columns (1 row per BGC)
    gecco_df["Sample_ID"] = sample
    gecco_df["BGC_length"] = gecco_df["BGC_end"] - gecco_df["BGC_start"]
    gecco_df["CDS_count"] = [
        len(gecco_df["CDS_ID"].iloc[i].split(";")) for i in range(0, gecco_df.shape[0])
    ]  # Number of contigs in 'Annotation_ID'
    gecco_df["Prediction_tool"] = "GECCO"

    # Add column 'InterPro_ID'
    for gbk_path in gbk_paths:
        bgc_id = gbk_path.split("/")[-1][0:-4]
        gecco_df.loc[gecco_df["cluster_id"] == bgc_id, "InterPro_ID"] = getInterProID(
            gbk_path
        )

    # Add empty columns with no output from GECCO
    gecco_df["BGC_complete"] = "NA"
    gecco_df["MIBiG_ID"] = "NA"
    if not gecco_df.dropna(how="all").empty:
        gecco_out = pd.concat([gecco_out, gecco_df])

    # Fill all empty cells with NA
    for row in range(len(gecco_df["PFAM_domains"])):
        if gecco_out["PFAM_domains"].isnull().values[row]:
            gecco_out.loc[row, "PFAM_domains"] = "NA"

    # Return data frame with ordered columns
    gecco_out = gecco_out[summary_cols]

    if verbose:
        print("Done.")

    return gecco_out



#########################################
# FUNCTION: Filter bgcs
#########################################
def filter_bgc(table, min_length, contig_edge):
    """
    Filters a DataFrame of BGC (Biosynthetic Gene Cluster) data based on length, edge proximity, and tool-specific criteria.

    Args:
        table (pd.DataFrame): Input DataFrame containing BGC data with columns such as 'sample_id', 'BGC_length',
                            'Prediction_tool', 'BGC_probability', 'Product_class', 'contig_id', 'BGC_start', and 'BGC_end'.
        min_length (int): Minimum length of BGCs to retain.
        contig_edge (int): Threshold distance from contig edges to exclude BGCs.

    Returns:
        pd.DataFrame: Filtered DataFrame with:
            - Non-header rows (excludes rows where "sample_id" equals "sample_id").
            - BGCs meeting minimum length and probability criteria.
            - Excludes BGCs near contig edges and with "Product_class" as "Saccharide".
    """
    df = table.copy()
    # filter out headings, short bgcs by deepbgc, bgcs that could be cut off, and rows with column names
    df = df.loc[df["sample_id"] != "sample_id"]
    df = df.drop(df.loc[df["BGC_length"] < min_length].index)
    df = df.drop(
        df.loc[
            (df["Prediction_tool"] == "deepBGC") & (df["BGC_probability"] < 0.6)
        ].index
    )
    df = df.drop(df.loc[(df["Product_class"] == "Saccharide")].index)
    df = df.drop(
        df.loc[
            (
                df["contig_id"].str.extract(r"length_(\d+)")[0].astype(int)
                - df["BGC_end"]
            ).abs()
            <= contig_edge
        ].index
    )
    df = df.drop(df.loc[df["BGC_start"].isin(range(0, contig_edge + 1))].index)
    return df



#########################################
# FUNCTION: Clean/standardize table
#########################################
def cleanup_table(table):
    """
    Cleans up DataFrame by creating a unique identifier (sample_id + contig_id)for each row and standardizing
    the "Product_class" column.

    Args:
        table (pd.DataFrame): The input DataFrame containing columns "sample_id", "contig_id",
                              and "Product_class".

    Returns:
        pd.DataFrame: A cleaned and sorted DataFrame with updated "identifier" and "Product_class" columns.
    """

    df = table.copy()
    df["identifier"] = (
        df["sample_id"]
        .str.replace("-megahit", "", regex=False)
        .str.replace("-metaspades", "", regex=False)
        + "_"
        + df["contig_id"].str.replace(r"(NODE_\d+).*", r"\1", regex=True)
    )

    # standardize the product class
    df["Product_class"] = df["Product_class"].str.replace(
        "NRP;Polyketide", "NRP-Polyketide"
    )
    df["Product_class"] = df["Product_class"].str.replace("NRPS", "NRP")
    df["Product_class"] = df["Product_class"].str.strip()
    df["Product_class"] = df["Product_class"].str.replace("NRP-like", "NRP")
    df["Product_class"] = df["Product_class"].str.replace("NRPS-like", "NRP")
    df["Product_class"] = df["Product_class"].str.replace("NRP-Polyketide", "NRP")
    df["Product_class"] = df["Product_class"].str.replace("RiPP-like", "RiPP")
    df["Product_class"] = df["Product_class"].str.replace(
        "Polyketide-Terpene", "Polyketide"
    )
    df["Product_class"] = df["Product_class"].str.replace(
        r"^Lanthipeptide-class-\w+$", "Lanthipeptide", regex=True
    )
    df = df.sort_values(by="identifier")

    return df



#########################################
# FUNCTION: Combining prediction tools
#########################################
def combine_tool(table, bgc, representative):
    """
    Updates the given table DataFrame to mark the presence of a specific prediction tool.
    It searches for in the table for the line which equals to the representative's identifier and adds the tool
    for the bgcs that is given as an argument

    Args:
        table (pd.DataFrame): A DataFrame containing BGC data with an 'identifier' column.
        bgc (dict): A dictionary representing a BGC with a 'Prediction_tool' key.
        representative (dict): A dictionary representing the representative BGC with an 'identifier' key.

    Returns:
        The function modifies the 'table' DataFrame in place, setting the appropriate tool
        column ('deepBGC', 'GECCO', or 'antiSMASH') to "Yes" and updating the
        'Tool_representative'.
    """
    if bgc["Prediction_tool"] == "deepBGC":
        table.loc[table["identifier"] == representative["identifier"], "deepBGC"] = (
            "Yes"
        )
    elif bgc["Prediction_tool"] == "GECCO":
        table.loc[table["identifier"] == representative["identifier"], "GECCO"] = "Yes"
    elif bgc["Prediction_tool"] == "antiSMASH":
        table.loc[table["identifier"] == representative["identifier"], "antiSMASH"] = (
            "Yes"
        )
    table.loc[
        table["identifier"] == representative["identifier"], "Tool_representative"
    ] = representative["Prediction_tool"]



########################
# DEDUPLCATION FUNCTIONS
########################
def process_chunk(chunk):
    """
    Worker function for process_parallel.
    This function takes a chunk of bgcs (all on the same contig) and:
    1. If only one bgc is present, add it to the table.
    2. If multiple bgcs are present, find the longest one as representative and
       add it to the table. Then, combine all overlapping bgcs with the representative
       and add them to the table. If a bgc does not overlap, add it separately to the table.

    Args:
        chunk (pd.DataFrame): Part of the DataFrame that contains all bgcs on the same contig.

    Returns:
        pd.DataFrame: A cleaned and sorted DataFrame with updated "identifier" and "Product_class" columns.
    """
    filtered_bgcs_chunk = pd.DataFrame(
        columns=[
            "deepBGC",
            "GECCO",
            "antiSMASH",
            "merged",
            "Product_class",
            "mmseqs_lineage_contig",
        ]
    )  # create df for chunk

    identifier = chunk["identifier"].iloc[0]
    chunk = chunk.reset_index(drop=True)
    same_contig = {}
    for i in range(len(chunk)):
        # Assign a unique key for each entry
        bgc_identifier = f"{identifier}_{i}"
        same_contig[bgc_identifier] = chunk.iloc[i].to_dict()

    # make sure that id matches the key (_1)
    for key, value in same_contig.items():
        value["identifier"] = key

    # add individual bgc to table
    if len(same_contig) == 1:
        single_bgc = pd.DataFrame.from_dict(same_contig, orient="index")
        filtered_bgcs_chunk = pd.concat(
            [filtered_bgcs_chunk, single_bgc], ignore_index=True
        )
        combine_tool(filtered_bgcs_chunk, single_bgc.iloc[0], single_bgc.iloc[0])
        filtered_bgcs_chunk.loc[
            filtered_bgcs_chunk["identifier"] == single_bgc.iloc[0]["identifier"],
            "merged",
        ] = 1

    # combine overlapping bgcs on same contig
    elif len(same_contig) > 1:
        # find the longest bgc as representative
        representative = max(same_contig.values(), key=lambda x: x["BGC_length"])
        # create new table with all representatives and combined with other bgcs
        bgc_new = pd.DataFrame([representative])
        filtered_bgcs_chunk = pd.concat(
            [filtered_bgcs_chunk, bgc_new], ignore_index=True
        )
        duplicates = len(same_contig)
        product_classes = []
        outside_counter = 0

        # iterate over bgcs on the same contig
        for key, value in same_contig.items():
            if key != representative["identifier"]:  # dont compare to itself
                # bgc overlaps the representative
                if (
                    representative["BGC_start"] < value["BGC_end"]
                    and representative["BGC_end"] > value["BGC_start"]
                ):
                    combine_tool(filtered_bgcs_chunk, representative, representative)
                    combine_tool(filtered_bgcs_chunk, value, representative)
                    if (
                        pd.notna(representative["Product_class"])
                        and representative["Product_class"] not in product_classes
                    ):
                        product_classes.append(representative["Product_class"])
                    if (
                        pd.notna(value["Product_class"])
                        and value["Product_class"] not in product_classes
                    ):
                        product_classes.append(value["Product_class"])

                # bgc lays outside of representative
                else:
                    value = pd.DataFrame.from_dict(value, orient="index").T

                    value.loc[0, "contig_id"] = (
                        value.loc[0, "contig_id"] + f"_{chr(65 + outside_counter)}"
                    )
                    outside_counter += 1

                    filtered_bgcs_chunk = pd.concat(
                        [filtered_bgcs_chunk, value], ignore_index=True
                    )
                    duplicates -= 1
                    filtered_bgcs_chunk.loc[
                        filtered_bgcs_chunk["identifier"] == value.loc[0, "identifier"],
                        "merged",
                    ] = "1"
                    filtered_bgcs_chunk.loc[
                        filtered_bgcs_chunk["identifier"] == value.loc[0, "identifier"],
                        "Product_class",
                    ] = value.loc[0]["Product_class"]

                    combine_tool(filtered_bgcs_chunk, value.iloc[0], value.iloc[0])
                    combine_tool(filtered_bgcs_chunk, representative, representative)
        filtered_bgcs_chunk.loc[
            filtered_bgcs_chunk["identifier"] == representative["identifier"],
            "Product_class",
        ] = ", ".join(
            map(str, product_classes)
        )  # convert to str, so join works
        filtered_bgcs_chunk.loc[
            filtered_bgcs_chunk["identifier"] == representative["identifier"], "merged"
        ] = str(duplicates)
    return filtered_bgcs_chunk


def parallelization(table, cores):
    """
    Parallelizes the processing of a DataFrame by dividing it into chunks based on the same contig_id
    and processes each chunk in parallel using the specified number of CPU cores. The processed chunks
    are then combined, duplicates are removed, and the final DataFrame is reformatted.

    Args:
        table (pd.DataFrame): The input DataFrame containing data to be processed.
        cpu (int): The number of CPU cores to use for parallel processing.

    Returns:
        pd.DataFrame: The processed and reformatted DataFrame with duplicates removed and columns reordered.
    """
    df = table.copy()
    grouped_df = df.groupby("identifier")
    # Convert the grouped DataFrame to a list of groups
    chunks = [group for _, group in grouped_df]

    if verbose:
        print("\nDereplicating the files\n... ", end="")

    # Use ProcessPoolExecutor to parallelize the processing
    max_workers = cores  # number of cores
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = executor.map(
            process_chunk, chunks
        )  # each group (chunk) goes to the process_chunk function and gets saved in results

    # Combine results from all processes, drop duplicates from the overlapping chunks, and reformat the final table
    filtered_bgcs = pd.concat(results, ignore_index=True)
    filtered_bgcs = filtered_bgcs.drop(["Prediction_tool"], axis=1)
    filtered_bgcs = filtered_bgcs.drop_duplicates(subset=["identifier"], keep="last")
    filtered_bgcs["Product_class"] = filtered_bgcs["Product_class"].apply(
        lambda x: ", ".join(sorted(x.split(", "))) if isinstance(x, str) else x
    )
    filtered_bgcs = filtered_bgcs.reindex(
        [
            "identifier",
            "sample_id",
            "contig_id",
            "BGC_start",
            "BGC_end",
            "BGC_length",
            "deepBGC",
            "GECCO",
            "antiSMASH",
            "merged",
            "Product_class",
            "Tool_representative",
            "BGC_probability",
            "mmseqs_lineage_contig",
            "BGC_region_contig_ids",
            "MIBiG_ID",
            "InterPro_ID",
            "CDS_count",
            "BGC_complete",
            "PFAM_domains",
            "CDS_ID",
        ],
        axis=1,
    )
    if verbose:
        print("Done.")
    return filtered_bgcs



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

if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=FutureWarning, message="The behavior of DataFrame concatenation with empty or all-NA entries is deprecated.*")

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
                    summary_antismash_temp = antismash_workflow(input_antismash)
                    summary_antismash = pd.concat(
                        [summary_antismash, summary_antismash_temp]
                    )
            else:
                summary_antismash = antismash_workflow(input_antismash)
        elif tool == "deepBGC":
            summary_deepbgc = deepbgc_workflow(input_deepbgc)
        elif tool == "GECCO":
            summary_gecco = gecco_workflow(input_gecco)

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


    ## ADD THE FILTERING PART HERE
    df = summary_all.copy()

    df["BGC_length"] = pd.to_numeric(df["BGC_length"], errors="coerce", downcast="integer")
    df["BGC_start"] = pd.to_numeric(df["BGC_start"], errors="coerce", downcast="integer")
    df["BGC_end"] = pd.to_numeric(df["BGC_end"], errors="coerce", downcast="integer")
    df["BGC_probability"] = pd.to_numeric(df["BGC_probability"], errors="coerce", downcast="integer")

    df = filter_bgc(df, min_length, contig_edge)
    df = cleanup_table(df)
    filtered_bgcs = parallelization(df, cores)
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