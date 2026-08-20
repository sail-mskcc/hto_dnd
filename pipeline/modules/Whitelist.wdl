version 1.0

import "Translate.wdl" as Translate
import "Exceptions.wdl" as Exceptions

task Prepare {

    ### --- DOCS/INPUTS --- ###
    meta {
        description: "Prepare cell barcode whitelist by unzipping, removing '-1' suffix and selecting first column"
        author: "Tobias Krause"
        email: "krauset@mskcc.org"
    }

    parameter_meta {
        path_whitelist: { help: "Input file should have one barcode per line, optionally gzipped. For example, use `filtered_feature_bc_matrix/barcodes.tsv.gz` from 10x Genomics Cell Ranger output." }
    }

    input {
        File path_whitelist
        Map[String, String] docker_versions
    }

    ### --- COMMAND --- ###
    command <<<
        set -euo pipefail

        mkdir -p results

        # unzip if necessary and select first column
        if [[ "~{path_whitelist}" == *.gz ]]; then
            zcat ~{path_whitelist} | cut -f1 > results/whitelist.txt
        else
            cat ~{path_whitelist} | cut -f1 > results/whitelist.txt
        fi

        # write a boolean if '-1' is in barcodes (first line check)
        if head -n 1 results/whitelist.txt | grep -q -- '-1$'; then
            echo true > results/check_suffix.txt
        else
            echo false > results/check_suffix.txt
        fi

        # remove '-1' suffix if present
        sed 's/-1$//' results/whitelist.txt > results/whitelist_nosuffix.txt

        # count number of barcodes
        echo $(wc -l < results/whitelist.txt) > results/number_barcodes.txt
    >>>

    output {
        File whitelist = "results/whitelist.txt"
        File whitelist_nosuffix = "results/whitelist_nosuffix.txt"
        Int number_barcodes = read_int("results/number_barcodes.txt")
        Boolean has_1_suffix = read_boolean("results/check_suffix.txt")
    }

    runtime {
        docker: "docker.io/sailmskcc/hto-dnd-utils:" + docker_versions["hto-dnd-utils"]
    }
}
