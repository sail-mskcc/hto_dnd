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



def numpy_to_python(obj):
    """Convert NumPy types to native Python types.

    Args:
        obj (object): The object to convert.

    Returns:
        object: The converted object.
    """
    if isinstance(obj, np.generic):
        return obj.item()
    elif isinstance(obj, dict):
        return {k: numpy_to_python(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [numpy_to_python(i) for i in obj]
    else:
        return obj


def write_stats(result_df, metrics, output_file="stats.yml"):
    """Write statistics and metrics to a YAML file.

    Args:
        result_df (pd.DataFrame): The result DataFrame containing the hashID and Doublet_Info columns.
        metrics (dict): A dictionary containing the metrics for each HTO.
        output_file (str): The output file path. Default is 'stats.yml'.

    Returns:
        None
    """
    stats = result_df.groupby(by="hashID").size().to_dict()
    stats["Total"] = len(result_df)

    # Convert NumPy values to native Python types
    metrics = numpy_to_python(metrics)

    output_dict = {"stats": stats, "metrics": metrics}

    # Write stats and metrics to the YAML file
    with open(output_file, "wt") as fout:
        yaml.dump(output_dict, fout, sort_keys=False, default_flow_style=False)
