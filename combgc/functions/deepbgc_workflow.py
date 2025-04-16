from Bio import SeqIO
import pandas as pd
import os


########################
# DEEPBGC FUNCTIONS
########################


def deepbgc_workflow(deepbgc_path, verbose):
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
