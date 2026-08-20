version 1.0

workflow GetGeometry {
    input {
        String chemistry
    }

    Map[String, Map[String, String]] geometry_dict = {
        "Biolegend TotalSeq A V3": {
            "read_geometry": "2[1-16]",
            "bc_geometry": "1[1-16]",
            "umi_geometry": "1[17-28]",
            "orientation": "fw",
            "barcode_translation": "false",
            "barcode_translation_list": "",
            "library_type": "ISR"
        },
        "Biolegend TotalSeq B V2": {
            "read_geometry": "2[10-25]",
            "bc_geometry": "1[1-16]",
            "umi_geometry": "1[17-28]",
            "orientation": "fw",
            "barcode_translation": "true",
            "barcode_translation_list": "10x v3.1",
            "library_type": "ISR"
        },
        "Biolegend TotalSeq B V3": {
            "read_geometry": "2[10-25]",
            "bc_geometry": "1[1-16]",
            "umi_geometry": "1[17-28]",
            "orientation": "fw",
            "barcode_translation": "true",
            "barcode_translation_list": "10x v4",
            "library_type": "ISR"
        },
        "Biolegend TotalSeq C V3": {
            "read_geometry": "2[10-25]",
            "bc_geometry": "1[1-16]",
            "umi_geometry": "1[17-28]",
            "orientation": "fw",
            "barcode_translation": "false",
            "barcode_translation_list": "",
            "library_type": "ISR"
        }
    }

    output {
        String read_geometry = geometry_dict[chemistry]["read_geometry"]
        String bc_geometry = geometry_dict[chemistry]["bc_geometry"]
        String umi_geometry = geometry_dict[chemistry]["umi_geometry"]
        String orientation = geometry_dict[chemistry]["orientation"]
        String library_type = geometry_dict[chemistry]["library_type"]
        Boolean barcode_translation = (geometry_dict[chemistry]["barcode_translation"] == "true")
        String barcode_translation_list = geometry_dict[chemistry]["barcode_translation_list"]
    }
}
