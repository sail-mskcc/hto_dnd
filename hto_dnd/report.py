def create_report():
    # Save outputs (try catch to return the adata object even if saving fails)
    try:
        # create paths (first, so that we can save the paths in the metadata)
        paths = {}
        path_viz = os.path.join(os.getcwd(), "dsb_viz.png")
        if path_adata_out is not None:
            os.makedirs(os.path.dirname(path_adata_out), exist_ok=True)
            path_viz = os.path.join(
                os.path.dirname(path_adata_out),
                os.path.basename(path_adata_out).split(".")[0] + "_dsb_viz.png",
            )
            paths["adata_denoised"] = path_adata_out
        if create_viz:
            os.makedirs(os.path.dirname(path_viz), exist_ok=True)
            paths["viz"] = path_viz
        adata = add_meta(adata, step="paths", **paths)

        # save adata
        if path_adata_out is not None:
            adata.write_h5ad(path_adata_out)
            logger.info(f"AnnData object saved to '{path_adata_out}'")

        # save visualization
        if create_viz:
            create_visualization(adata, path_viz)
            logger.info(f"Visualization plot saved to '{path_viz}'")

    except Exception as e:
        logger.error(f"Failed to save outputs: '{e}'")
