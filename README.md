# nmr-novelty-check

Has the NMR spectrum of this compound already been published?

`nmrcheck` takes a structure (or a reaction) as SMILES, looks it up in open
collections of experimental NMR data, and writes a report that says whether
the compound has been reported, how sure the tool is, and what the answer
rests on. You do not need any spectra or datasets of your own: on first use
the tool downloads the public data, builds a local index and then answers
from that index.

```
$ nmrcheck check c1ccccc1
reported (confidence 0.95)
Experimental NMR data for this exact structure are present in the indexed public sources.
json: reports/c1ccccc1.json
md: reports/c1ccccc1.md
```

The full report for this run is in [examples/benzene.md](examples/benzene.md).

## Automatic download verified

`nmrcheck check <SMILES>` on a fresh clone with no local data will download NMRexp from Zenodo and nmrshiftdb2 from SourceForge, verify the NMRexp checksum, build the local index, and run the check. This was verified end-to-end on 2026-11-12: benzene (`c1ccccc1`) returned `reported` with confidence 0.95 after an automatic download and index build.

If the download step fails on your network (firewall, proxy, or a changed upstream URL), fall back to manual download and `nmrcheck build --from-dir "<path>" --offline`.

## Install

Python 3.10 or newer.

```
git clone https://github.com/yingjieshanjrh-coder/nmr-novelty-check
cd nmr-novelty-check
pip install -e .
```

## Use

### One compound

```
nmrcheck check "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
```

The first call downloads about 820 MB of source data into `./data` and builds
the index. With the complete NMRexp set that took 12 minutes on a two-core
machine and produced a 1.1 GB SQLite file. Later calls answer in under a
second and need no network.

For a quick trial, `nmrcheck build --limit 200000` indexes nmrshiftdb2 plus
only the first 200 000 NMRexp records (under two minutes). The report then
carries a `source_truncated` entry.

### A reaction

Give a reaction SMILES (`reactants>agents>products`), or name the other
species with flags. The first product is the target. Every other species is
looked up as well, and the report lists where their signals fall on top of
the target's and which of their signals stay clear of it.

```
nmrcheck check "CC(=O)Cl.OCc1ccccc1>CCN(CC)CC>CC(=O)OCc1ccccc1" \
    -i "CC(=O)O" -i "CC(=O)OC(C)=O" -s CDCl3
```

`-r` adds a reactant, `-g` a reagent, `-i` a possible impurity, `-s` the NMR
solvent (its residual signals are checked against the target). Output of this
command: [examples/benzyl_acetate_reaction.md](examples/benzyl_acetate_reaction.md).

### Web form

```
nmrcheck serve
```

then open http://127.0.0.1:8000. The same report is available as JSON at
`/api/check?target=...`.

### Data you downloaded yourself

If the machine has no network access, download the two files listed under
[Data sources](#data-sources) by hand, put them in one folder and run

```
nmrcheck build --from-dir /path/to/folder --offline
```

`--data-dir` (or the `NMRCHECK_DATA` environment variable) moves the data
folder somewhere else.

## What the report contains

Each run writes a JSON file and a Markdown rendering of it.

| section | content |
| --- | --- |
| verdict | status, confidence in 0..1, and the reasons that built the confidence |
| target | canonical SMILES, InChIKey, formula, atom numbering |
| evidence | every matching record: experimental or calculated, exact or variant match, source, DOI or database id, solvent, peak list, agreement with the predicted spectrum |
| assignment | observed peak, atoms, predicted shift, difference, and how much training data stood behind the prediction |
| other species | the same look-up for reactants, reagents and impurities |
| overlap | signals of other species and of the solvent that coincide with the target's |
| data sources | file, checksum, licence, citation and record count per source; predictor hold-out error |
| degradation log | every point where the run had to settle for less (see below) |
| limitations | what the result cannot tell you |

### Verdicts

| status | meaning |
| --- | --- |
| `reported` | experimental data for the identical InChIKey exist and at least one peak list fits the predicted spectrum |
| `reported_unverified` | records for the identical InChIKey exist but none fits the prediction; the pairing of structure and spectrum in the source may be wrong |
| `reported_as_variant` | only another stereoisomer or isotopologue (same first InChIKey block) has data |
| `listed_without_experiment` | the structure is stored only with calculated spectra |
| `not_found` | nothing in the index. This does not mean the compound is new |

### Confidence

The number is a fixed, readable rule rather than a fitted probability:

- identical InChIKey with experimental data: 0.70
- at least one of those peak lists agrees with the prediction (score >= 0.5): +0.10
- two or more independent documents: +0.10
- present in more than one database: +0.05
- no peak list agrees with the prediction: -0.25
- variant match only: 0.50, +0.10 if the peak lists agree
- not found: 0.50 x (core sources loaded / core sources expected). It never
  goes higher, because absence from two databases says little about the rest
  of the literature.

The agreement score between a predicted and an observed peak list is
`exp(-MAE / s) x (share of observed peaks that were matched)`, with
s = 3 ppm for 13C and 0.3 ppm for 1H.

### Degradation log

| code | raised when |
| --- | --- |
| `source_unavailable` | a download failed, or offline mode was on and the file was absent |
| `source_not_indexed` | a known source was missing when the index was built |
| `source_truncated` | the index was built with `--limit` |
| `source_unreadable` | a source file could not be parsed |
| `checksum_mismatch` | a downloaded file differs from the published checksum |
| `structures_skipped` | SMILES in a source that RDKit rejected |
| `predictor_untrained` | no assigned spectra were available to train the predictor |
| `prediction_fallback` | an atom had no trained environment at all; a fixed atom-type value was used |
| `prediction_low_radius` | signals predicted from an environment narrower than two bonds |
| `prediction_not_independent` | the compound itself is in the predictor's training data |
| `assignment_predicted_only` | no experimental record to assign, predicted shifts listed instead |
| `records_truncated` | more records than the report lists |
| `invalid_structure` | a SMILES could not be parsed |
| `multiple_products` | a reaction had several products; the first became the target |
| `solvent_unknown` | no residual-signal entry for the given solvent |

## How it works

1. **Fetch.** Raw files are downloaded with resume support and, where the
   publisher gives one, verified against the MD5 checksum.
2. **Index.** Every record is reduced to InChIKey, canonical SMILES, formula,
   nucleus, solvent, peak list and a pointer to its origin, and stored in
   SQLite. Structure normalisation runs in parallel worker processes.
3. **Predict.** The assigned spectra in nmrshiftdb2 train a look-up predictor:
   each atom is described by its circular environments of radius 0 to 4 bonds
   (Morgan identifiers), and the median shift of the widest environment seen
   in training is the prediction. This is the HOSE-code idea with RDKit
   environments. Ten percent of the compounds are held out during the build
   to measure the error, which is stored in the index and printed in every
   report. On the data described below it was 2.70 ppm (13C, 43 611 signals)
   and 0.28 ppm (1H, 12 560 signals) overall, and 1.10 ppm / 0.14 ppm where a
   radius-4 environment was available.
4. **Compare.** Predicted and observed peak lists are paired with the
   Hungarian algorithm; predicted signals left over share the nearest
   observed peak, since accidental overlap is normal in real spectra.
5. **Decide.** The verdict and confidence follow the rules above.

The comparison in step 4 matters because NMRexp was extracted from PDFs
automatically. In the benzene example one NMRexp record pairs benzene with a
triplet at 0.82 ppm and a doublet at 0.97 ppm; that record gets an agreement
score of 0.00 and does not count as support.

## Validation

| Compound | In training set? | Confidence | Verdict | What it demonstrates |
|---|---|---|---|---|
| benzene (`c1ccccc1`) | Yes | 0.95 | `reported` | In-distribution match; report flags `prediction_not_independent`, so not evidence of generalization |
| ibuprofen (`CC(C)Cc1ccc(cc1)C(C)C(=O)O`) | Yes | 0.95 | `reported` | Same caveat; also exposes 0.83 ppm error on the exchangeable COOH proton |
| `CC1=C2C(=O)N(c3ccc(Cl)cc3)C(=S)N2CC(O)O1` | No | 0.45 | `reported_unverified` | Out-of-distribution; predicted spectrum does not fit the single experimental record; tool downgrades instead of overclaiming |

> `reported_unverified` means a record with an identical InChIKey exists, but none of its peak lists fits the predicted spectrum. This should **not** be read as "the structure is wrong". It means either the literature extraction paired the structure with the wrong peak list, or the predictor is unreliable for this skeleton. The JSON report's `degradations` field (shown as "Degradation log" in the Markdown report) currently only flags predictor unreliability via `prediction_low_radius`; it cannot yet distinguish extraction errors, so human review is required.

## Data sources

| source | file | size | licence | fetched automatically |
| --- | --- | --- | --- | --- |
| [nmrshiftdb2](https://nmrshiftdb.nmr.uni-koeln.de) | [`nmrshiftdb2withsignals.sd`](https://sourceforge.net/projects/nmrshiftdb2/files/data/) | 159 MB | nmrshiftdb2 Data License (open content) | yes |
| [NMRexp](https://zenodo.org/records/17296666) | `NMRexp_10to24_1_1004.parquet` | 661 MB | CC BY 4.0 | yes |
| [NMRBank](https://doi.org/10.1039/D4SC08802F) | `NMRBank.csv` | - | CC BY-NC 3.0 | no |

Indexed from these files: 76 440 spectra from nmrshiftdb2 and 3 351 159 from
NMRexp, covering 1 547 488 distinct compounds.

NMRBank is optional. It has no stable direct download link and its licence is
non-commercial, so the tool does not fetch it. If you place a CSV at
`data/nmrbank/NMRBank.csv` it is read by a generic reader that looks for a
SMILES column and a 13C shift column; this reader has only been tested on
synthetic input, not on the published NMRBank file. When the file is absent,
every report says so in the degradation log.

No source data are redistributed in this repository.

Please cite the data when you use results:

- Kuhn S., Schlörer N. E. *Magn. Reson. Chem.* 2015, 53, 582 (nmrshiftdb2)
- Wang J.-J. et al. NMRexp: A database of 3.3 million experimental NMR
  spectra. *Sci. Data* 2025; Zenodo 10.5281/zenodo.17296666
- Wang Q. et al. *Chem. Sci.* 2025, 16, 11548 (NMRBank)

## Limitations

- `not_found` means "not in these two databases". NMRexp covers supporting
  information of papers from 2010 to 2024; older literature, patents and
  commercial databases such as SciFinder or Reaxys are outside its reach. Use
  the result as a first screen, not as a substitute for a literature search.
- Matching is by InChIKey. Salts, solvates, tautomers and other charge states
  count as different compounds.
- The predictor ignores solvent, concentration, temperature and
  stereochemistry, and is poor for exchangeable protons and for environments
  that nmrshiftdb2 does not cover. Check the "env. radius" and "training hits"
  columns before trusting a predicted value.
- For compounds that are themselves in nmrshiftdb2 the prediction is largely
  a recall of the stored assignment. The report flags this.
- Assignments are proposals based on the nearest predicted shift.
- Only 13C and 1H are predicted and compared. Records for 19F, 31P, 11B and
  29Si are indexed but not shown in the report yet.
- Input is SMILES only; compound names are not resolved.
- The tool never sees a spectrum of your sample. It answers "has this been
  published", not "is my product what I think it is".

## Related work

Tools that start from a measured spectrum and search for the structure:
[Sherlock](https://doi.org/10.3390/molecules28031448),
[MixONat](https://acs.figshare.com/collections/MixONat_a_Software_for_the_Dereplication_of_Mixtures_Based_on_sup_13_sup_C_NMR_Spectroscopy/5026139),
[VersaDB](https://doi.org/10.1002/cmtd.202300020), lucy-ng, and the search
form of nmrshiftdb2 itself. This project goes the other way: it starts from
the structure you intend to make and asks what is already on record.

## Development

```
pip install -e ".[dev]"
pytest
```

The tests build a small index from synthetic records written in the same
formats as the real sources, and exercise the download code against a local
HTTP server, so they run without network access.

## Licence

MIT, see [LICENSE](LICENSE). The licence covers the code only; the data
sources keep their own licences.
