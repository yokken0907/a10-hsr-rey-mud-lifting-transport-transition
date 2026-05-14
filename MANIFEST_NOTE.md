# Manifest policy

This v0.1.1-public-gate archive uses a manifest-excluding-self policy.

`FILE_MANIFEST.csv` and `FILE_MANIFEST.json` intentionally exclude the manifest files themselves. This avoids self-referential hash instability and prevents false mismatch reports after manifest regeneration.

The selected source-material directory was normalized to an ASCII-safe path:

`source_materials/selected_from_uploaded_hsr_project/hsr_simulation_rareearth_selected/`
