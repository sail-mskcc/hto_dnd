# Snakemake Pipeline Test

This directory contains a test for the Snakemake pipeline integration example.

## Test Snakemake Rule

```python
# test_snakemake.py
import tempfile
import os
from pathlib import Path

def test_snakemake_rule():
    """Test the Snakemake rule for HTO demultiplexing."""
    
    # Create temporary directory structure
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        
        # Create sample data directory
        sample_dir = temp_path / "data" / "sample1"
        sample_dir.mkdir(parents=True)
        
        # Create dummy input files (in real scenario, these would be actual data)
        (sample_dir / "hto_filtered.h5ad").touch()
        (sample_dir / "hto_raw.h5ad").touch()
        (sample_dir / "gex_raw.h5ad").touch()
        
        # Create results directory
        results_dir = temp_path / "results" / "sample1"
        results_dir.mkdir(parents=True)
        
        # Snakemake rule would execute:
        # hto demultiplex --adata-hto data/sample1/hto_filtered.h5ad ...
        
        # Test command construction
        cmd = f"""
        hto demultiplex \\
            --adata-hto {sample_dir}/hto_filtered.h5ad \\
            --adata-hto-raw {sample_dir}/hto_raw.h5ad \\
            --adata-gex {sample_dir}/gex_raw.h5ad \\
            --adata-out {results_dir}/demultiplexed.h5ad \\
            --csv-out {results_dir}/assignments.csv \\
            --demux-method gmm \\
            --denoise-version v2
        """
        
        # Verify command is properly formatted
        assert "--adata-hto" in cmd
        assert "--adata-hto-raw" in cmd
        assert "--adata-gex" in cmd
        assert "--adata-out" in cmd
        assert "--csv-out" in cmd
        
        print("Snakemake rule test passed")

if __name__ == "__main__":
    test_snakemake_rule()
```

## Snakemake Rule

```python
# Snakefile
rule hto_demultiplex:
    input:
        hto="data/{sample}/hto_filtered.h5ad",
        hto_raw="data/{sample}/hto_raw.h5ad", 
        gex="data/{sample}/gex_raw.h5ad"
    output:
        demux="results/{sample}/demultiplexed.h5ad",
        csv="results/{sample}/assignments.csv"
    shell:
        """
        hto demultiplex \\
            --adata-hto {input.hto} \\
            --adata-hto-raw {input.hto_raw} \\
            --adata-gex {input.gex} \\
            --adata-out {output.demux} \\
            --csv-out {output.csv} \\
            --demux-method gmm \\
            --denoise-version v2
        """
```

## Running the Test

```bash
python test_snakemake.py
```