# NMR novelty check

Query: `c1ccccc1`  
Generated: 2026-11-12T00:05:55+00:00 (nmr-novelty-check 0.1.0)

## Verdict

**reported**, confidence 0.95

Experimental NMR data for this exact structure are present in the indexed public sources.

- 13 experimental record(s) with an identical InChIKey (base 0.70)
- 11 record(s) fit the predicted spectrum (+0.10)
- 3 independent documents (+0.10)
- found in more than one database (+0.05)

## Target

| field | value |
| --- | --- |
| canonical SMILES | `c1ccccc1` |
| InChIKey | UHOVQNZJYSORNB-UHFFFAOYSA-N |
| formula | C6H6 |
| atom numbering | `[cH:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1` |

Records found:

| nucleus | exact, experimental | distinct DOIs | exact, calculated | stereo/isotope variants |
| --- | --- | --- | --- | --- |
| 13C | 4 | 1 | 0 | 0 |
| 1H | 9 | 2 | 0 | 0 |

### Evidence, 13C

| evidence | match | source | reference | solvent | agreement | MAE / ppm | peaks / ppm |
| --- | --- | --- | --- | --- | --- | --- | --- |
| experimental | exact | nmrshiftdb2 | nmrshiftdb2 molecule 7901 | CDCl3 | 1.00 | 0.00 | 128.50 |
| experimental | exact | nmrshiftdb2 | nmrshiftdb2 molecule 7901 | - | 1.00 | 0.00 | 128.50 |
| experimental | exact | nmrshiftdb2 | nmrshiftdb2 molecule 7901 | - | 0.97 | 0.09 | 128.59 |
| experimental | exact | nmrexp | [10.1021/acscatal.1c03375](https://doi.org/10.1021/acscatal.1c03375) | - | 0.96 | 0.12 | 128.62 |

### Evidence, 1H

| evidence | match | source | reference | solvent | agreement | MAE / ppm | peaks / ppm |
| --- | --- | --- | --- | --- | --- | --- | --- |
| experimental | exact | nmrshiftdb2 | nmrshiftdb2 molecule 7901 | CDCl3 | 1.00 | 0.00 | 7.36 (6H) |
| experimental | exact | nmrshiftdb2 | nmrshiftdb2 molecule 7901 | CD3COCD3 | 1.00 | 0.00 | 7.36 (6H) |
| experimental | exact | nmrshiftdb2 | nmrshiftdb2 molecule 7901 | DMSO-d6 | 0.97 | 0.01 | 7.37 (6H) |
| experimental | exact | nmrshiftdb2 | nmrshiftdb2 molecule 7901 | CD3CN | 0.97 | 0.01 | 7.37 (6H) |
| experimental | exact | nmrshiftdb2 | nmrshiftdb2 molecule 7901 | CD3OD | 0.91 | 0.03 | 7.33 (6H) |
| experimental | exact | nmrexp | [10.1021/acscatal.1c03375](https://doi.org/10.1021/acscatal.1c03375) | - | 0.77 | 0.08 | 7.28 (s, 6H) |
| experimental | exact | nmrshiftdb2 | nmrshiftdb2 molecule 7901 | CCl4 | 0.65 | 0.13 | 7.23 (6H) |
| experimental | exact | nmrshiftdb2 | nmrshiftdb2 molecule 7901 | C6D6 | 0.50 | 0.21 | 7.15 (6H) |
| experimental | exact | nmrexp | [10.1021/ja202900h](https://doi.org/10.1021/ja202900h) | CDCl3 | 0.00 | 6.39 | 0.82 (t, 3H), 0.97 (d, 3H) |

### Assignment, 13C

| atoms | observed | predicted | delta | nuclei | env. radius | training hits |
| --- | --- | --- | --- | --- | --- | --- |
| C1,2,3,4,5,6 | 128.50 | 128.50 | 0.00 | 6 | 4 | 18 |

### Assignment, 1H

| atoms | observed | predicted | delta | nuclei | env. radius | training hits |
| --- | --- | --- | --- | --- | --- | --- |
| C1,2,3,4,5,6 | 7.36 | 7.36 | 0.00 | 6 | 4 | 44 |

The predictor was trained on data that include this compound; observed and predicted values are not independent.

## Data sources

| source | file | spectra indexed | licence | md5 | citation |
| --- | --- | --- | --- | --- | --- |
| nmrshiftdb2 | nmrshiftdb2withsignals.sd | 76440 | nmrshiftdb2 Data License (open content) | d8fc88de76f1 | Kuhn S., Schloerer N.E., Magn. Reson. Chem. 2015, 53, 582 |
| nmrexp | NMRexp_10to24_1_1004.parquet | 3351159 | CC BY 4.0 | 2f9ed8bc5333 | Wang J.-J. et al., Sci. Data 2025; Zenodo 10.5281/zenodo.17296666 |

Index built 2026-11-11T23:52:38+00:00, 1547488 compounds.
Predictor hold-out error (10% of nmrshiftdb2 compounds): 13C: MAE 2.701 ppm over 43611 signals; 1H: MAE 0.28 ppm over 12560 signals.

## Degradation log

| code | subject | detail |
| --- | --- | --- |
| source_not_indexed | `nmrbank` | NMRBank.csv was not available when the index was built |
| prediction_not_independent | `c1ccccc1` | this compound is part of the predictor's training data, so agreement with nmrshiftdb2 records is partly circular |

## Limitations

- Coverage is limited to the indexed sources. NMRexp is drawn from supporting information of 2010-2024 papers; compounds published elsewhere, in patents, or in subscription databases (SciFinder, Reaxys) are not seen.
- NMRexp was extracted automatically from PDFs. Structures can be paired with the wrong paragraph and peak lists can contain recognition errors.
- Structures are compared by InChIKey. Salts, solvates, tautomers and different charge states are treated as different compounds.
- Predicted shifts come from an environment look-up trained on nmrshiftdb2. Solvent, concentration, temperature and stereochemistry are ignored, and exchangeable protons are unreliable.
- Peak assignment is by nearest predicted shift. It is a proposal, not a 2D-NMR-grade assignment.
- No user spectrum is involved: the tool checks whether data have been published, it does not verify a sample.
