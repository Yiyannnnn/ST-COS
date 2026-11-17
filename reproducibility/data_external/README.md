# External data

Public reference data are not committed. Download them from the sources listed
in `../EXTERNAL_INPUTS.txt`, convert each dataset to the shared reference-based
format (`counts.csv` and `coords.csv`; see `../README.md`), and place the
prepared files here:

```text
data_external/prepared/
  DLPFC_151507/counts.csv, coords.csv
  ...                                    # 12 DLPFC sections
  MERFISH_X0_04/counts.csv, coords.csv
  ...                                    # fields X0_04, X0_09, X0_14, X0_19, X0_24
  MOB_StereoSeq/counts.csv, coords.csv
```

scCube additionally needs its pretrained model checkpoint and the ordered
training-gene list for each dataset; these model artifacts are not
redistributed.
