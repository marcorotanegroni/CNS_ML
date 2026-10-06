# Binary concordance: agreement composition and threshold sensitivity

The analysis addresses threshold sensitivity and the composition of the original
patient-state agreement. It uses the same 5,881 matched TCGA patients and 58
signatures under four rules. No matched patients were excluded for missingness.

## Definitions

For each signature pair, n11 counts active-active patients, n00 inactive-inactive,
and n10 and n01 discordant classifications. Let A = n11 + n00 be the number of
concordant patients and N the total number of matched patients.

- Original patient-state Jaccard: A / (2N - A).
- Active-active share of observed agreements: n11 / A.
- Inactive-inactive share of observed agreements: n00 / A.
- Active-state Jaccard: n11 / (n11 + n10 + n01).
- Expected active-active patients under independence: N * p1 * p2, where p1 and
  p2 are the active prevalences; reported as observed/expected n11 and as
  Cohen's kappa, which uses the same marginal prevalences.

The two shares sum to one when A > 0. Their denominator is the concordant
patients for that pair, not all 5,881 patients. If A = 0, both shares are recorded
as NA and patient-state Jaccard is zero. Shares must be shown alongside the
original Jaccard and counts: they describe what agreements consist of, not how
frequently patients agree. Shared inactivity is an observed classification,
not automatically an artefact or evidence of biological equivalence.

## Thresholds

The original rule uses each signature's median in its full compendium
(Drews: 6,335; Steele: 9,699; Tao: 10,370), before restriction to matched patients.
If the median is zero, activity is x > 0; otherwise it is x >= median. The zero
exception prevents zero exposure from being classified as active. Inactivity
means falling below this operational rule and need not mean zero exposure or
absence of copy-number alterations in the tumor.

Three alternatives use the 25th, 50th, and 75th percentiles of strictly positive
exposures in the same full-compendium reference populations. These perturb
sparse signatures whose unconditional quantiles can remain zero. Quantiles use
R type 7 and include positive threshold ties; zero exposures stay inactive.
Thresholds and achieved active prevalences are reported in thresholds.csv.

## Results

Original zero/median rule; percentages below refer only to concordant patients:

| Pair | Patient-state J | Concordant patients A | Active-active n11 | Inactive-inactive n00 | Active share | Inactive share |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| CN1 / Sig1 | 0.765 | 5,098 | 2,014 | 3,084 | 39.5% | 60.5% |
| CN1 / Sig2 | 0.716 | 4,908 | 1,855 | 3,053 | 37.8% | 62.2% |
| CX1 / CN1 | 0.576 | 4,299 | 1,898 | 2,401 | 44.1% | 55.9% |
| CN1 / CN2 | 0.150 | 1,531 | 487 | 1,044 | 31.8% | 68.2% |

For CN4/CN10, all 5,829 observed agreements are inactive-inactive: the inactive
share is 100% and patient-state Jaccard is 0.982. In contrast, CN1/Sig1 and
CN1/Sig2 contain both coactivation and shared inactivity. The 68.2% inactive
share for CN1/CN2 accompanies low overall agreement; the share alone must not
be interpreted as high concordance.

Active-state Jaccard (original rule): CN1/Sig1 0.720, CN1/Sig2 0.656,
CX1/CN1 0.545, CN1/CN2 0.101, CN4/CN10 0. The 76 pairs with patient-state
Jaccard above 0.90 all involve sparse Steele signatures, have more than 99% of
their agreements inactive-inactive and an active-state Jaccard below 0.09.

Threshold sensitivity, summarized on the 48 cross-compendium pairs among the 12
selected signatures (`principal_pair_ranks.csv`). Rank among the 48 pairs
(1 = most concordant), patient-state / active-state Jaccard:

| Pair | Original | Pos. Q25 | Pos. Q50 | Pos. Q75 |
| --- | ---: | ---: | ---: | ---: |
| CN1 / Sig1 | 1 / 1 | 1 / 1 | 2 / 10 | 2 / 11 |
| CN1 / Sig2 | 2 / 2 | 2 / 2 | 1 / 6 | 1 / 10 |
| CX1 / CN1 | 6 / 5 | 19 / 14 | 15 / 13 | 17 / 7 |

Median active prevalence of the 12 selected signatures: 49%, 50%, 27% and 11%
under the four rules. At the stringent rules the patient-state Jaccard is
inflated by shared inactivity (median inactive share 0.88 and 0.99).

The pair table also contains expected active-active counts under independence
(`expected_n11`, `observed_expected_n11`) and Cohen's kappa for reference; they
are not used in the manuscript or the response.

## Reproduction and outputs

Run from the repository root:

```bash
Rscript code/revision/01_binary_concordance.R
```

The final section of code/01_signature_exploration.Rmd runs the same analysis.
Use --no-plots to generate tables only.

- state_agreement_components.png: patient-state Jaccard, inactive-inactive share of agreements and active-state Jaccard under the original rule, for all 58 and the 12 selected signatures.
- threshold_sensitivity_principal_pairs.png: patient-state and active-state Jaccard of the 48 cross-compendium selected pairs under the four rules, with the pairs discussed in the manuscript and their ranks.
- principal_pair_ranks.csv: rank of each pair discussed in the manuscript among the 48 selected cross-compendium pairs, per rule, for patient-state and active-state Jaccard.
- pair_metrics.csv: raw counts, number of concordant patients, original Jaccard, agreement shares and marginal active prevalences for all 6,612 pair/rule combinations.
- principal_pairs.csv: the four manuscript comparisons and the sparse example CN4/CN10 under each rule, with all counts and metrics.
- thresholds.csv: exposure thresholds, reference populations, ties and achieved prevalences.
- provenance.txt and session_info.txt: input checksum, calculation conventions and software versions.

All figures use a fixed blue scale from 0 to 1 and fixed signature ordering.
The color mapping distinguishes low values rather than using the white plateau
below 0.25 in the original Figure 4. The original figure calculation is retained.

Validation reproduces the original binary matrix exactly and independently
checks the counts, patient-state Jaccard and agreement shares, including the
case with no concordant patients. Final manuscript integration and supplementary
figure/table numbering remain editorial steps.
