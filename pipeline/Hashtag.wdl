version 1.0

import "modules/Whitelist.wdl" as Whitelist
import "modules/Alevin.wdl" as Alevin
import "modules/Translate.wdl" as Translate
import "modules/Demultiplex.wdl" as Demultiplex
import "modules/BasicQC.wdl" as BasicQC
import "modules/Exceptions.wdl" as Exceptions
import "modules/Geometries.wdl" as Geometries

workflow Hashtag {

    input {

        String sample_name
        String chemistry
        File path_taglist
        File path_whitelist_raw
        File path_whitelist_filtered
        File path_gex_raw
        Array[File] path_r1_files
        Array[File] path_r2_files

        Map[String, String] docker_versions = {
            "hto-dnd-utils": "2.0.0",
            "hto-basic-qc": "1.1.0",
            "hto-dnd": "1.2.0"
        }
    }

    # ASSERTIONS
    # whitelist has no header
    # r1 and r2 files must have matching order
    # taglist has correct format
    # maybe even that tags can be found in the first x bases of read2?
    
    # GET GEOMETRY
    call Geometries.GetGeometry as Geometry { input: chemistry = chemistry }

    ### --- PREPARE WHITELISTS --- ###
    call Whitelist.Prepare as WhitelistRaw {
        input:
            path_whitelist = path_whitelist_raw,
            docker_versions = docker_versions,
    }

    call Whitelist.Prepare as WhitelistFiltered {
        input:
            path_whitelist = path_whitelist_filtered,
            docker_versions = docker_versions,
    }

    if (Geometry.barcode_translation) {
        call Translate.TranslateBarcodes as WhitelistTranslated {
            input:
                whitelist = WhitelistRaw.whitelist,
                chemistry = Geometry.barcode_translation_list,
                docker_versions = docker_versions,
        }
    }
    File barcodes_translated = select_first([WhitelistTranslated.barcodes_translated, WhitelistRaw.whitelist_nosuffix])

    ### --- ALIGN & QUANTIFY --- ###
    call Alevin.Quantify {
        input:
            path_r1_files = path_r1_files,
            path_r2_files = path_r2_files,
            path_taglist = path_taglist,
            path_whitelist_raw = barcodes_translated,
            read_geometry = Geometry.read_geometry,
            bc_geometry = Geometry.bc_geometry,
            umi_geometry = Geometry.umi_geometry,
            orientation = Geometry.orientation,
            library_type = Geometry.library_type,
            docker_versions = docker_versions,
    }

    ### --- TRANSLATE AND CREATE H5AD --- ###
    if (Geometry.barcode_translation) {
        call Translate.TranslateBarcodes as BarcodesTranslated {
            input:
                whitelist = Quantify.barcodes,
                chemistry = Geometry.barcode_translation_list,
                docker_versions = docker_versions,
        }
    }
    File barcodes = select_first([BarcodesTranslated.barcodes_translated, Quantify.barcodes])

    call Alevin.ToAdata {
        input:
            mtx = Quantify.mtx,
            barcodes = barcodes,
            features = Quantify.features,
            path_whitelist_filtered = WhitelistFiltered.whitelist,  # should include "-1" if it's also in whitelist_raw
            has_1_suffix = WhitelistRaw.has_1_suffix,  # adds "-1" back if necessary
            docker_versions = docker_versions,
    }

    ### --- HTO-DND --- ###
    call Demultiplex.HtoDnd {
        input:
            adata_raw = ToAdata.adata_raw,
            adata_filtered = ToAdata.adata_filtered,
            adata_gex_raw = path_gex_raw,
            docker_versions = docker_versions,
    }


    ### --- BASIC QC --- ###
    call BasicQC.BasicQC {
        input:
            sample_name = sample_name,
            adata_demux = HtoDnd.adata_demux,
            cell_stats = Quantify.cell_stats,
            alignment_info = Quantify.alignment_info,
            template_html = "hashtag_report.html",
            docker_versions = docker_versions,
    }

    # NOTE: See here for saturation calculation: https://github.com/COMBINE-lab/simpleaf/issues/157


    output {
        File html_report = BasicQC.html_report
        File cell_stats = Quantify.cell_stats
        File alignment_info = Quantify.alignment_info

        File whitelist_gex_raw = path_whitelist_raw

        File adata_hto_demux = HtoDnd.adata_demux
        File adata_hto_raw = ToAdata.adata_raw
    }
}
