# comBGC
<img src="combgc/images/com-bgc-logo.png" alt="Werner Siemens Foundation" width="250" />

comBGC is a tool that can be used to combine the results of the prediction tools deepBGC, GECCO, and antiSMASH for detecting bionsynthetic gene clusters (BGCs). It combines the results of all three tools and adds different filters, so only hight quality BGCs are included. Furthermore, a merging step is applied to the datasets. This removes the redundancy and combines the results of the prediction tool. This way no duplicate BGCs are included in the final table.
### Installation
---
To use comBGC, you first need to create an environment using **conda**:
    
    conda create -n combgc

Activate the environment

    conda activate combgc

and install the required packages using **pip**

    pip install -r requirements.txt

### Using comBGC
---
To run comBGC and parse the results of multiple BGC prediction tools into a single table you need to specify diffent parameters which can be displayed with using `-h` :

    python3 combgc/comBGC.py --help

Using --helpt gives you an overview of all parameters that can be modified, of which `-i`, `-o`, and `--cores` are required.

    arguments:
        -h, --help            show this help message and exit
        -i [PATH(s) [PATH(s) ...]], --input [PATH(s) [PATH(s) ...]]
                                path(s) to the required output file(s) of antiSMASH, DeepBGC and/or GECCO
                                these can be:
                                - antiSMASH: <sample name>.gbk and (optional) knownclusterblast/ directory
                                - DeepBGC:   <sample name>.bgc.tsv
                                - GECCO:     <sample name>.clusters.tsv
                                Note: Please provide files from a single sample only. If you would like to
                                summarize multiple samples, please see the --antismash_multiple_samples flag.
        -o [PATH], --outdir [PATH]
                                directory for comBGC output. Default: current directory
        -a [PATH], --antismash_multiple_samples [PATH]
                                directory of antiSMASH output. Should contain subfolders (one per
                                sample). Can only be used if --input is not specified.
        -vv, --verbose        increase output verbosity
        -v, --version         show version number and exit
        --cores CORES         number of CPU cores to use
        --min_length LENGTH   Minimum length [bp] of BGC sequence to keep. Default: 3000
        --contig_edge BASES   Exclude BGCs within X bases of the contig edge. Default: 2

**Example Usage**

    python3 combgc/comBGC.py \
    -i tests/antismash_example/ECO010-megahit/ECO010-megahit.gbk \
    tests/deepbgc_example/ECO010-megahit/ECO010-megahit.bgc.tsv \
    tests/gecco_example/ECO010-megahit/ECO010-megahit.clusters.tsv \
    -o ./output/ --cores 10


### Authors
---
This tool was initiated by Jasmin Frangenberg (@jasmezz) and expanded by Tom Richtermeier (@tomrichtermeier).

### Funding
---
This project was funded by Werner Siemens Foundation grant ‘Palaeobiotechnology’

<img src="combgc/images/wss.svg" alt="Werner Siemens Foundation" width="250" />

