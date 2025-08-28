"""Run external tool: GMM-Demux.

Source: https://github.com/CHPGenetics/GMM-Demux
"""

import os
import subprocess
import tempfile

import numpy as np
import pandas as pd

from hto._exceptions import UserInputError
from hto._logging import get_logger


def gmm_demux(
    df: pd.DataFrame,
    hash_id: str = "hash_id",
    to_counts: bool = None,
    params: str = "",
    verbose: int = 1,
    **kwargs,
) -> pd.DataFrame:
    """Run GMM-Demux on filtered HTO. GMM-Demux uses raw counts. However, the input to this function is normalised, log-transformed and denoised.

    Use exponentiation to get an approximation of the cleaned counts.

    Args:
        df (pd.DataFrame): Dataframe with HTO counts
        hash_id (str): Hash id for output
        to_counts (bool): Whether to convert logged values to counts. If none, try to infer from the data.
        params (str): Parameters for GMM demux
        verbose (int): Verbosity level
        **kwargs: Additional arguments

    """
    # assertions
    prohibited_params = ["-f ", "--full ", "-o ", "--output ", "--csv ", "-c "]
    for param in prohibited_params:
        if param in params:
            raise UserInputError(
                f"Prohibited parameter '{param}' found in params. Please remove it."
            )

    # set tmp path
    path_tmp_dir = os.environ.get(
        "DND_TMPDIR", os.path.join(os.path.expanduser("~"), ".dnd_tmp")
    )
    os.makedirs(path_tmp_dir, exist_ok=True)
    path_tmp = tempfile.TemporaryDirectory(dir=path_tmp_dir, delete=False)
    path_tmp_str = path_tmp.name
    logger = get_logger("gmm-demux", level=verbose)

    # set to_counts
    # - if integers, then False
    # - if rounded floats equal to integers, then False
    # - if floats, then True
    # - else, raise error
    if to_counts is None:
        if all([np.issubdtype(t, np.integer) for t in df.dtypes]):
            to_counts = False
        if np.all(df.iloc[:1000].round(0) == df.iloc[:1000]):
            to_counts = False
        elif all([np.issubdtype(t, np.floating) for t in df.dtypes]):
            to_counts = True
        else:
            raise UserInputError(
                f"Cannot infer 'to_counts' from data types: {df.dtypes}. Please set 'to_counts' explicitly."
            )

    # transform to counts
    if to_counts:
        df = np.exp(df) - 1

    # write dataframe
    path_df = os.path.join(path_tmp_str, "hto.csv")
    df.to_csv(path_df, index=True, header=True)

    # get hto_array
    hto_array = kwargs.get("hto_array", ",".join(df.columns))

    # run gmm
    cmd = f"cd {path_tmp_str} && GMM-demux {path_df} {hto_array} --output {path_tmp_str} --full FULL {params} --csv"
    try:
        _ = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        logger.info("GMM demux completed successfully.")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"GMM demux failed with return code {e.returncode}.\n"
            f"Command: '{cmd}'\n"
            f"STDOUT: {e.stdout}\n"
            f"STDERR: {e.stderr}"
        )

    # read assignment
    # - mark doublets with : instead of -
    # - note that celltypes may include -, so first replace those with <dash>, then replace - with : and finally replace <dash> with -
    df_config = pd.read_csv(
        os.path.join(path_tmp_str, "FULL/GMM_full.config"),
        index_col=None,
        sep=", ",
        header=None,
        names=["Cluster_id", hash_id],
        engine="python",
    )
    for c in df.columns:
        df_config[hash_id] = df_config[hash_id].str.replace(c, c.replace("-", "<dash>"))
    df_config[hash_id] = df_config[hash_id].str.replace("-", ":")
    df_config[hash_id] = df_config[hash_id].str.replace("<dash>", "-")

    # join assignment
    df_full = pd.read_csv(os.path.join(path_tmp_str, "FULL/GMM_full.csv"), index_col=0)
    df_out = df_full.join(
        df_config, on="Cluster_id", how="left", lsuffix="_full", rsuffix="_config"
    )

    # rename index to "barcode"
    df_out.index.name = "barcode"

    # clean up
    path_tmp.cleanup()

    return df_out[[hash_id]]
