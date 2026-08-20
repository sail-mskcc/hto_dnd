version 1.0

task BasicQC {
    ### --- DOCS/INPUTS --- ###
    meta {
        description: "Basic QC report for HTO/ADT data using HTO-Basic-QC"
        author: "Tobias Krause"
        email: "krauset@mskcc.org"
    }

    parameter_meta {
        sample_name: { help: "Sample name" }
        adata_demux: { help: "h5ad-File containing demultiplexed ADT/HTO count matrix." }
        cell_stats: { help: "File containing cell features obtained from aleving-fry (featureDump.txt)." }
        template_html: { help: "HTML template file for the report." }
    }

    input {
        String sample_name
        File adata_demux
        File cell_stats
        File alignment_info
        String template_html
        Map[String, String] docker_versions
    }

    ### --- COMMAND --- ###
    command <<<
        set -euo pipefail

        mkdir -p results

        python /opt/render/cli.py render ~{template_html} \
                --sample-name ~{sample_name} \
                --path-h5ad ~{adata_demux} \
                --path-stats ~{cell_stats} \
                --path-alignment-info ~{alignment_info} \
                --path-output results/hashtag_report.html
    >>>

    output {
        File html_report = "results/hashtag_report.html"
    }

    runtime {
        docker: "docker.io/sailmskcc/hto-basic-qc:" + docker_versions["hto-basic-qc"]
    }
}