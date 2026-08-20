version 1.0

task HtoDnd {
    ### --- DOCS/INPUTS --- ###
    meta {
        description: "Demultiplex HTO/ADT data using HTO-DND: https://github.com/sail-mskcc/hto_dnd/"
    }

    parameter_meta {
        adata_raw: { help: "h5ad-File containing raw ADT/HTO count matrix." }
        adata_filtered: { help: "h5ad-File containing filtered ADT/HTO count matrix." }
        adata_gex_raw: { help: "h5ad-File or h5-File containing raw gene expression count matrix." }
    }

    input {
        File adata_raw
        File adata_filtered
        File adata_gex_raw
        Map[String, String] docker_versions
    }

    ### --- COMMAND --- ###
    command <<<
        set -euo pipefail

        mkdir -p results

        hto demultiplex \
            --adata-hto ~{adata_filtered} \
            --adata-hto-raw ~{adata_raw} \
            --adata-gex ~{adata_gex_raw} \
            --adata-out results/adata_demux.h5ad \
    >>>

    output {
        File adata_demux = "results/adata_demux.h5ad"
    }

    runtime {
        docker: "docker.io/sailmskcc/hto-dnd:" + docker_versions["hto-dnd"]
    }

}