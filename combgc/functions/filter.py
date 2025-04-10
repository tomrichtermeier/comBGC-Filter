import pandas as pd
import os
import re
from concurrent.futures import ProcessPoolExecutor



#########################################
# FUNCTION: Filter BGCs
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
            - Excludes header-like rows where ‘sample_id’ is literally the string ‘sample_id’.
            - BGCs meeting minimum length and probability criteria.
            - Excludes BGCs near contig edges and with "Product_class" as "Saccharide".
    """
    df = table.copy()
    df = df.loc[df["sample_id"] != "sample_id"] # remove header rows
    df = df[df["BGC_length"] >= min_length] # filter by length

    # Filter out BGCs by deepBGC with low probability (easy solution brought weird filter behavior)
    df_antiSMASH = df[df["Prediction_tool"] == "antiSMASH"]
    df_deepbgc = df[df["Prediction_tool"] == "deepBGC"]
    df_gecco = df[df["Prediction_tool"] == "GECCO"]
    # Apply probability filter to DeepBGC only
    df_deepbgc = df_deepbgc.drop(df_deepbgc.loc[df_deepbgc["BGC_probability"] < 0.6].index)
    # Recombine DeepBGC and GECCO before applying other filters
    df = pd.concat([df_deepbgc, df_gecco, df_antiSMASH], ignore_index=True)

    #Remove BGCs near contig edges
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

    df = df.drop(df.loc[(df["Product_class"] == "Saccharide")].index)
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


def parallelization(table, cores, verbose):
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
