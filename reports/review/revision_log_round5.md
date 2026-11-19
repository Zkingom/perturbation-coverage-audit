# Revision log, round 5 (2026-10-28): convergence to manuscript v1.0

Scope frozen: no new analysis, no new claim, no change to any frozen protocol. All changes below are formatting, numbering, rounding and wording.

## 1. Consistency with the constructive reference (step 1)
Unchanged judgement. The lowest deviation (RPE1, linear-PCA, T5: −2.8458 pp against the simulated reference 0.508) is about −1.2 SD of the 21-split pooled spread (simulated single-split SD 10.77 pp, independence approximation, conservative), and the confirmatory 95% CIs of the primary-model T5 coverage contain the simulated reference in all three datasets. The supporting sentence was already in Results v0.10 and is retained.

## 2. Rounding and precision (step 2)
- Rule: round half away from zero, directly from unrounded archived values. Proportions (coverage, power, R, CV) 3 decimals; pp, percentages and their CIs 1 decimal; ratios, correlations, AUROC, log-odds and the CIs of proportions/ratios 2 decimals; protocol parameters as specified. Stated in Methods ("Rounding").
- Interpretation recorded: the default "CIs to 2 decimals" was applied to CIs of proportions, ratios and correlations; CIs of pp differences follow the pp rule (1 decimal); R and CV were treated as proportions (3 decimals).
- Changes: 18 explicit main-text edits from raw values (e.g. Abstract 0.48 → 0.480, 0.52 → 0.518; development marginal 0.90 → 0.904; power 0.96 → 0.960, 0.81 → 0.809; SE 0.57 pp → 0.6 pp; noise-floor T5 0.55/0.49/0.52/0.66 → 0.552/0.492/0.518/0.664; grid M0 T5 0.49 → 0.491; R 0.83 → 0.830; ρ 0.110/0.119 → 0.11/0.12; reliability "≈0.99; 0.992" → "0.99 in all three datasets"); 24 CIs with 3 decimals rounded to 2 (5 of them ending in 5 were re-rounded from raw values: [0.95, 0.96], [0.68, 0.74], [0.68, 0.80], [0.80, 0.85], [0.88, 0.90]); 5 Supplementary Note 1 coverages raised to 3 decimals from raw values; "~15 pp" → "14.6 pp". Sentence-level diff of every changed number: `reports/review/numeric_sentence_diff_round5.txt` (458 lines).
- CI notation unified to "(95% CI a to b)" in running text; "[a, b]" retained inside tables. Minus signs before numbers set as "−".
- Not harmonised (listed for transparency): correlations and other statistics reported to 3 decimals inside the tables of Supplementary Note 1 §1–§6 were left as computed; supplementary tables S1–S10 are generated directly from their source files with the precision rule applied per column type.

## 3. Numbering (step 3): order of first citation
Supplementary tables (old → new):

| old | new | content |
|---|---|---|
| S3 | S1 | Freeze timeline and deviations |
| S1 | S2 | Reference-method per-split width matching |
| S2 | S3 | Confirmatory decision tables (now cited in Methods) |
| S8 | S4 | Constructive reference by model (v2, with archived Monte Carlo) |
| S4 | S5 | Selection of evaluation units (<30 cells) |
| S6 | S6 | Score–effect-strength relation |
| S5 | S7 | Number-to-source map (regenerated on final text; now cited in Methods) |
| S10 | S8 | GEARS coupling test (new; previously proposed as S10) |
| S7 | S9 | Prospectively identifiable subsets |
| S9 | S10 | CORUM subsets |

Supplementary figures (old → new):

| old | new | content |
|---|---|---|
| Fig. S2 | Fig. S1 | Score–effect-strength relation (figS_scale) |
| Fig. S1 | Fig. S2 | Noise-floor correction (figS_noise_floor) |
| Fig. S3 / 'S_explore' | Fig. S3 | Exploratory modularity (figS_explore_modularity) |

Main figure: Fig. 1 (unchanged); a citation of Fig. 1c was added in Results. No main tables. Citations added where tables were previously uncited: decision tables (Methods, confirmatory analysis), number-to-source map (Methods, reproducibility), GEARS table (Results). The GEARS table is S8 under first-citation order (the round-4 proposal had S10; the user allowed either).

## 4. Terminology, references, structure
- "non-essential genes" → "non-library targets" (2 places, consistent with the definition in Methods); "split-half oracle" → "split-half reference"; "top-tier" → "T5"; "split conformal prediction" → "split-conformal prediction".
- Citations unified to one numbered list in order of first citation (12 references); per-section lists removed; inline "(Li et al., Sci. Adv. 2026, doi:…)" replaced by a numbered citation. Reference years follow the CrossRef print record: GEARS 2024 (online 2023), Systema 2026 (online 2025); journal names abbreviated; page ranges from CrossRef.
- Structure: Title → Abstract → Introduction → Methods → Results → Discussion → Limitations (moved out of Discussion) → Data and code availability → References → Supplementary information (separate file).
- Internal draft notes removed. The three unconfirmable statements about the GEARS run are kept with unchanged content, rendered in English as "[not confirmable from the records …]" because the manuscript is in English.
- Temporal uses of "first" rephrased (4 places) so the banned-word check is clean; none was a priority claim.
- Supplementary: "FIRST READ" → "INITIAL READ" in the timeline; enumerations "First/Second" → "(i)/(ii)"; the stale sentence describing the GEARS run as "planned" now states that it was run and detected no reduction (Supplementary Table S8); "GWPS" → "non-library target set" in the supplement (11 places).
- Main text: the calibration-methods table in Methods was given a caption and is cited as Table 1 (the only main table).
