# Data

## Reference sequence

This file contains Spike's reference sequence:

```
cov2_spike_wt.fasta
```

---

## Mutation escapes

### Raw files

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


### Aggregated table

This file contains mutation escape values for all available types of antibodies:

```
aggregated_mut_escapes.csv
```

Columns:

- `site`, `wildtype`, `mutation`: same as above;
- All other column names correspond to the available antibodies types. So each row contains mutation escape values for all available antibodies.

  **NOTE:** If there is no data for a particular antibodies type, then `-9` is placed!
