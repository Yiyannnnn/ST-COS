# External-data layout

Place downloaded reference data under this directory without changing the
subdirectory names below. These files are intentionally not committed.

```text
data_external/
  DLPFC/
    layers/
      151507_layers.txt
      ...
      151676_layers.txt
    151507/
      151507_filtered_feature_bc_matrix.h5
      spatial/
    ...
  MERFISH/
    count/
      X0_04_count.csv
      X0_09_count.csv
      X0_14_count.csv
      X0_19_count.csv
      X0_24_count.csv
    mapping/
      X0_04_mapping.csv
      X0_09_mapping.csv
      X0_14_mapping.csv
      X0_19_mapping.csv
      X0_24_mapping.csv
  HER2/
  STEREOSEQ/
```

The DLPFC analyses use sections 151507-151510, 151669-151672, and
151673-151676. The MERFISH confirmation uses fields X0_04, X0_09, X0_14,
X0_19, and X0_24. See the cited data sources in the manuscript for accession
and preprocessing details.
