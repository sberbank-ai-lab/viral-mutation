# Data

## Reference sequence

This file contains Spike's reference sequence:

```
cov2_spike_wt.fasta
```

---

## Mutation escapes

The following files contain stats on mutation escapes for different types of antibodies:

```
AZ_cocktail_raw_data.txt
Ellebedy_invivo_raw_data.txt
MAP_paper_antibodies_raw_data.txt
REGN_and_LY-CoV016_raw_data.txt
human_sera_raw_data.txt
```

Important columns:

- `site`: mutation position (1-based, corresponds to the reference sequence from above);
- `wildtype`: reference amino acid;
- `mutation`: mutation amino acid;
- `mut_escape`: mutation escape value;
- `condition`: antibodies type.
