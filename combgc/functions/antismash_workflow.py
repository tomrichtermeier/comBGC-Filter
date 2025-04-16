from Bio import SeqIO
import pandas as pd
import os
import re


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


def antismash_workflow(antismash_paths, verbose):
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
