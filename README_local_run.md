# Local Run Guide

From the `Spatial_Transcriptome_Studio` directory:

```bash
pip install -r requirements.txt
python app_cli.py --input C:\path\to\sample.h5ad --output results\sample1 --analysis frame26 --gene-mode fixed
```

Run multiple samples from a manifest:

```bash
python app_cli.py --manifest examples\sample_manifest.csv --output results\batch --analysis frame26 --gene-mode fixed
```

Use a custom config:

```bash
python app_cli.py --input C:\path\to\sample.h5ad --output results\sample1_istz --analysis istz --gene-mode custom --config examples\example_config.yaml
```

Distance is currently spot-based. Physical distance calibration fields are present in the config for future implementation.
