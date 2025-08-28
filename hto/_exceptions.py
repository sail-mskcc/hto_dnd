import numpy as np


class UserInputError(Exception):
    pass


class AnnDataFormatError(Exception):
    def __init__(self, key, *args, **kwargs):
        message = self.get_message(key, *args, **kwargs)
        super().__init__(message)

    def get_message(self, key, *args, **kwargs):
        _func = {
            "adata_raw_missing_cells": self._adata_raw_missing_cells,
            "adata_no_overlapping_names": self._adata_filtered_too_few_cells,
            "adata_not_float": self._adata_not_float,
        }.get(key, self._default)
        return _func(*args, **kwargs)

    def _default(self, *args, **kwargs):
        args = ", ".join(args)
        kwargs = ", ".join([f"{k}={v}" for k, v in kwargs.items()])
        return f"Unknown error occurred. Args: {args}, kwargs: {kwargs}"

    def _adata_raw_missing_cells(self, empty_barcodes):
        return (
            f"Too few empty droplets detected: '{len(empty_barcodes)}'. "
            f"Make sure that 'adata_raw' contains cells that are not present in 'adata_filtered'."
        )

    def _adata_filtered_too_few_cells(self, n_cells):
        return f"Too few cells in 'adata_filtered': '{n_cells}'. "

    def _adata_not_float(self, x):
        msg = f"Input matrix must contain float values, got {x.dtype} ({', '.join(np.unique(x).astype(str)[:10])})."
        return msg
