from Bio import SeqIO
import pandas as pd
import re


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


def gecco_workflow(gecco_paths, verbose):
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