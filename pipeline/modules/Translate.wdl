version 1.0

task TranslateBarcodes {
    ### --- DOCS/INPUTS --- ###
    meta {
        description: "Translate cell barcode whitelist from TotalSeq B to 10x format"
    }

    parameter_meta {
        chemistry: { help: "Chemistry used for gene expression. Supported: '10x V3.1', '10x V4'" }
        whitelist: { help: "File containing cell barcode whitelist to be translated. Ideally output from Whitelist.Prep task." }
    }

    input {
        String chemistry
        File whitelist
        Map[String, String] docker_versions
    }

    ### --- COMMAND --- ###
    command <<<
        set -euo pipefail

        mkdir -p results

        hto-dnd-utils translate-whitelist \
            --input ~{whitelist} \
            --chemistry "~{chemistry}" \
            --output results/whitelist.txt
    >>>

    output {
        File barcodes_translated = "results/whitelist.txt"
    }

    runtime {
        docker: "docker.io/sailmskcc/hto-dnd-utils:" + docker_versions["hto-dnd-utils"]
    } 
}
