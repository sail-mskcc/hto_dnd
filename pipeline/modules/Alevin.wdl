version 1.0

# Quantify ADT/HTO data using an Alevin pipeline with the following steps:
# 
# 1. Prepare taglist
# 2. Build index
# 3. Quantify
# 4. Correct cellbarcodes using provided whitelists
# 5. Collapse cellbarcodes
# 6. Create cell x gene matrix
# 7. Create AnnData object
# 8. Filter AnnData object
# 
# The steps are roughly following this tutorial: https://combine-lab.github.io/alevin-fry-tutorials/2021/af-feature-bc/

task Quantify {
    ### --- DOCS/INPUTS --- ###
    meta {
        description: "Quantify ADT/HTO data using Salmon Alevin."
        author: "Tobias Krause"
        email: "krauset@mskcc.org"
    }

    parameter_meta {
        path_r1_files: { help: "Array of paths to R1 FASTQ files." }
        path_r2_files: { help: "Array of paths to R2 FASTQ files." }
        path_taglist: { help: "Path to taglist file defining ADT/HTO tags. Must be tab-seperated with <tag_name> \\t <tag_sequence>." }
        path_whitelist_raw: { help: "Path to raw cell barcode whitelist file." }
    }

    input {
        # inputs
        Array[File] path_r1_files
        Array[File] path_r2_files
        File path_taglist
        File path_whitelist_raw
        String read_geometry
        String bc_geometry
        String umi_geometry
        String orientation
        String library_type
        Map[String, String] docker_versions
    }

    ### --- COMMAND --- ###
    # command
    command <<<
        set -euo pipefail

        # log the inputs:
        echo ">>INPUT: Read Geometry: ~{read_geometry}"
        echo ">>INPUT: BC Geometry: ~{bc_geometry}"
        echo ">>INPUT: UMI Geometry: ~{umi_geometry}"
        echo ">>INPUT: Orientation: ~{orientation}"
        echo ">>INPUT: Library Type: ~{library_type}"

        # setup paths
        mkdir -p results
        cd results
        
        # prepare t2g map: (<tag_name> \t <tag_name>)
        awk '{print $1"\t"$1;}' ~{path_taglist} > tag2gene.txt

        # build index
        salmon index \
            --transcripts ~{path_taglist} \
            --index index \
            --features \
            -k7

        # quantify
        # --libType: For example "ISR" - inward, stranded, reverse. more info: https://salmon.readthedocs.io/en/latest/library_type.html
        # --sketch: only align, no downstream quantification
        salmon alevin \
            --index index \
            --output quant_salmon \
            --libType ~{library_type} \
            --read-geometry ~{read_geometry} \
            --bc-geometry ~{bc_geometry} \
            --umi-geometry ~{umi_geometry} \
            -1 ~{sep=" " path_r1_files} \
            -2 ~{sep=" " path_r2_files} \
            --sketch

        # call cells
        # --expected-ori fw: forward orientation
        alevin-fry generate-permit-list \
            --expected-ori ~{orientation} \
            --input quant_salmon \
            --output-dir permit_list \
            --valid-bc ~{path_whitelist_raw}

        # collate
        alevin-fry collate \
            --rad-dir quant_salmon \
            --input-dir permit_list
        
        # quantify
        alevin-fry quant \
            --tg-map tag2gene.txt \
            --input-dir permit_list \
            --output-dir cxg \
            --resolution cr-like \
            --use-mtx

        # rename feature dump (for nice name only)
        cp cxg/featureDump.txt cxg/cell_stats.txt

        # alignment logs
        mkdir -p logs
        cp cxg/quant.json logs/
        cp permit_list/generate_permit_list.json logs/
        cp permit_list/collate.json logs/
        cp quant_salmon/cmd_info.json logs/
        cp index/info.json logs/
        cp index/versionInfo.json logs/

    >>>

    output {
        File mtx = "results/cxg/alevin/quants_mat.mtx"
        File barcodes = "results/cxg/alevin/quants_mat_rows.txt"
        File features = "results/cxg/alevin/quants_mat_cols.txt"
        File cell_stats = "results/cxg/cell_stats.txt"
        File alignment_info = "results/quant_salmon/aux_info/meta_info.json"
        Array[File] logs = glob("results/logs/*")
    }

    runtime {
        docker: "docker.io/sailmskcc/hto-dnd-utils:" + docker_versions["hto-dnd-utils"]
    }
}



task ToAdata {
    ### --- DOCS/INPUTS --- ###
    meta {
        description: "Convert cell x gene matrix in MTX format to AnnData object and filter using provided whitelist."
    }

    parameter_meta {
        mtx: { help: "Path to cell x gene matrix in MTX format." }
        barcodes: { help: "Path to barcodes file corresponding to the MTX matrix." }
        features: { help: "Path to features file corresponding to the MTX matrix." }
        path_whitelist_filtered: { help: "Path to filtered cell barcode whitelist file." }
        has_1_suffix: { help: "Boolean indicating if barcodes have '-1' suffix." }
    }

    input {
        File mtx
        File barcodes
        File features
        File path_whitelist_filtered
        Boolean has_1_suffix
        Map[String, String] docker_versions
    }

    ### --- COMMAND --- ###
    command <<<
        set -euo pipefail

        mkdir -p results

        # add '-1' suffix back if necessary
        if [ "~{has_1_suffix}" = true ]; then
            sed -i 's/$/-1/' ~{barcodes}
        fi

        # to adata
        hto-dnd-utils alevin-to-anndata \
            --path-mtx ~{mtx} \
            --path-barcodes ~{barcodes} \
            --path-features ~{features} \
            --path-output results/adata_raw.h5ad

        # filter adata
        hto-dnd-utils subset-anndata \
            --path-adata results/adata_raw.h5ad \
            --path-whitelist ~{path_whitelist_filtered} \
            --path-output results/adata_filtered.h5ad
    >>>

    output {
        File adata_raw = "results/adata_raw.h5ad"
        File adata_filtered = "results/adata_filtered.h5ad"
    }

    runtime {
        docker: "docker.io/sailmskcc/hto-dnd-utils:" + docker_versions["hto-dnd-utils"]
    }
}


