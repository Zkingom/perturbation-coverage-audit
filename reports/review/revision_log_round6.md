# Revision log, round 6 (2026-11-01): manuscript v1.1 working copy

Changes without new analysis (B2–B5 and layout):
- B2 Methods server section -> 3 sentences
- B2 Limitations item 9 compressed
- B3 Fig1 caption placeholder
- list M0
- list controls
- bold sigma_boot
- B3 Supp §2 contradiction
- Supp Note 2 added
- B3 S1 timeline +3 freeze rows
- B3 S2 dev columns filled
- B3 S7 summary v2
- Methods S7 wording
- B5 AI disclosure
- Supplementary Figs. S2 and S3 re-rendered from their archived plotting code with label changes only: 'K562 non-essential (GWPS)' → 'K562 non-library targets', 'GWPS' → 'non-library', 'oracle' → 'reference method'; coverage labels in Fig. S2a to three decimals (rounded half away from zero).
- 'on non-library target set' → 'on the non-library target set' (grammar).
- LaTeX headers: \XeTeXgenerateactualtext=1 (guarded) so the PDF text layer carries the original Unicode (hyphens and hatted symbols were mis-mapped in v1.0).
- Supplementary Table S7 regenerated with trace levels (see summary table); S1 timeline and S2 sources: reports/review/freeze_timeline_v2.csv, oracle_width_matching_by_split_v2.csv.

Round 6b (2026-11-03), after the Norman analysis:
- B1: new title, abstract and introduction (diagnostic storyline); new Methods section 'Combinatorial perturbations (protocol combo_crosscell v1.2)', Results section 'The diagnostic departs from zero on combinatorial perturbations' with Fig. 2, Discussion paragraph 'When stratified coverage is informative', limitation 11, data-availability entry.
- References renumbered by first citation; Norman 2019 and scPerturb added (now 14 entries).
- Supplementary Table S9 (Norman) inserted; former S9 → S10, S10 → S11.
- S1 timeline: wording of round 5 restored after the table was rebuilt ('GWPS' → 'non-library target set', 'FIRST READ' → 'INITIAL READ').
- Main text: 'a first wave' → 'an initial wave'.
- S7: 167 rows added for the numbers in the new sections.

Round 6c (2026-11-04): protocols v1.2.1 (Nadig source -> scPerturb, QC Q1-Q5, GEO fallback on rented server) and v1.3 (Norman unseen-component ablation; gene-level resampling) frozen at 2026-11-04T09:03:48Z; v1.3 computed. Manuscript v1.1: abstract and introduction softened (coverage differences are point-estimate directions; fully unseen pairs lose accuracy gain); Results paragraph on two-way resampling replaces pair-level dT5 intervals (pre-written rule B.4); new Results paragraph on the ablation; Methods section for v1.3; Discussion cites GEARS's Norman generalization classes (GEARS not run); limitation 11 updated; Fig. 2 Norman intervals now two-way; Table S9 extended plus ablation table; references 14-15 added (Davezies 2021; Owen 2007), Nadig renumbered 16; S7 +101 rows. v1.4 draft (second combinatorial dataset) written, not frozen.

- 2026-11-04 incident: code/download_resumable.py (hashed in freeze_v0.6/v0.6.1/v0.6.2 manifests) was accidentally overwritten when writing a new Zenodo downloader under the same name. It was reconstructed byte-for-byte by replaying the archived edits (edit_file calls of 2026-10-05 plus two Python edits, Windows CRLF line endings); the restored file matches the frozen SHA-256 315dbe2e…8cf9. The new downloader is code/download_zenodo_resumable.py.
- 2026-11-04 (6d): protocol v1.4 frozen at 2026-11-04T10:15:10Z (option 1, Wessels only; section 2.4 added before freezing: timing of the 'only one passes' rule, invariance to gates G3/G5, hard narrowing to Norman if G1 < 60 or G2 < 90%, two-way resampling and unchanged thresholds). Results: mechanistic sentence added (with one gene held out, the measured single-gene component alone reduced error more than linear-PCA, point estimates 0.267 vs 0.188, paired difference untested; elevation of R depends on measured single-gene components). Discussion: sentence on the second combinatorial screen (pending). Supplementary Table S13 (gate table) added; former S10/S11 → S11/S12 and the ablation table is S10. References 16 Wessels, 17 Joung added (Crossref-verified); Nadig → 18. S7 +6 rows.


## 2026-11-04 (6e): wording, file protection, version control, data provenance

- Wording weakened (user request): Results ablation sentence now 'is consistent with the elevation of R depending on …; the linear model showed no measurable gain beyond the measured component', with the two-unseen R verdict 'uncertain' stated; 'was no more accurate' / 'lost its accuracy advantage' → 'showed no measurable accuracy gain over the mean-response model' in Abstract, Introduction, Results and Discussion; 'The ablation locates the regime' → 'The ablation is consistent with …'. Presubmission draft aligned (gene-level CIs; no 'only there did coverage differ').
- File protection: every file whose current SHA-256 equals a hash recorded in any manifest (113 files, plus the large raw single-cell file on G: matched by size) and all manifests set read-only (60 were still writable). All tracked code/protocol files match a recorded hash; the files that differ from an earlier manifest are later frozen versions or recorded post-freeze engineering changes (v0.7 post_freeze_execution entries).
- Create-only rule: new files are written with code/safe_io.create_new (raises if the path exists); see project rule in memory.
- Version control: no git binary on this machine; a bare git repository _vcs_repo.git was created with dulwich (bare, so no git operation can write into the working directory). Retroactive annotated tags, one per freeze manifest (names: v0.5 … v1.4; the GEARS protocols v1.0/v1.1 are tagged gears_coupling-v1.0/-v1.1 to avoid ambiguity). Each tagged commit holds the manifest and the files whose current hash is recorded in it; commit and tag times are the actual creation time, the freeze time is in the tag message. Not pushed to any remote (user decides).

| tag | manifest | frozen (UTC) | tag created (UTC) | files in tagged commit | recorded but not in commit |
|---|---|---|---|---|---|
| v0.5 | freeze_v0.5.json | 2026-10-04T06:23:27Z | 2026-11-04T10:46:15Z | 10 + manifest | 9 (auditB_clusters_original_r8.csv, config.toml, go_ase100.npz, p0v3_units_all_targets.csv, pseudobulk_genes.csv, …) |
| v0.5.1 | freeze_v0.5.1.json | 2026-10-04T07:55:54Z | 2026-11-04T10:46:15Z | 10 + manifest | 9 (auditB_clusters_original_r8.csv, config.toml, go_ase100.npz, p0v3_units_all_targets.csv, pseudobulk_genes.csv, …) |
| v0.6 | freeze_v0.6.json | 2026-10-06T10:38:40Z | 2026-11-04T10:46:15Z | 9 + manifest | 13 (K562_gwps_raw_bulk_01.h5ad, K562_gwps_raw_singlecell_01.h5ad, go_ase100.npz, p0v3_units_all_targets.csv, pseudobulk_genes.csv, …) |
| v0.6.1 | freeze_v0.6.1.json | 2026-10-06T11:52:14Z | 2026-11-04T10:46:15Z | 15 + manifest | 19 (K562_gwps_raw_bulk_01.h5ad, K562_gwps_raw_singlecell_01.h5ad, analysis_v06.py, eval_v06.py, go_ase100.npz, …) |
| v0.6.2 | freeze_v0.6.2.json | 2026-10-08T11:53:59Z | 2026-11-04T10:46:15Z | 22 + manifest | 18 (K562_gwps_raw_bulk_01.h5ad, K562_gwps_raw_singlecell_01.h5ad, d1_run_v061.py, d1_run_v061.py_sha256_at_freeze, go_ase100.npz, …) |
| v0.7 | freeze_v0.7.json | 2026-10-10T09:05:55Z | 2026-11-04T10:46:16Z | 55 + manifest | 28 (K562_gwps_raw_bulk_01.h5ad, K562_gwps_raw_singlecell_01.h5ad, ReactomePathways.gmt.zip, confirm_run_v07.py, corum_human_current.txt, …) |
| gears_coupling-v1.0 | freeze_gears_coupling_v1.0.json | 2026-10-19T11:50:36Z | 2026-11-04T10:46:16Z | 11 + manifest | 2 (gears_bundle_v07.tar.gz, sha256) |
| gears_coupling-v1.1 | freeze_gears_coupling_v1.1.json | 2026-10-28T00:17:18Z | 2026-11-04T10:46:16Z | 3 + manifest | 0 () |
| v1.2 | freeze_v1.2.json | 2026-11-01T08:39:44Z | 2026-11-04T10:46:16Z | 12 + manifest | 0 () |
| v1.2.1 | freeze_v1.2.1.json | 2026-11-04T09:03:48Z | 2026-11-04T10:46:16Z | 2 + manifest | 0 () |
| v1.3 | freeze_v1.3.json | 2026-11-04T09:03:48Z | 2026-11-04T10:46:16Z | 5 + manifest | 2 (v12_result_sha256, v12_test_rows_sha256) |
| v1.4 | freeze_v1.4.json | 2026-11-04T10:15:10Z | 2026-11-04T10:46:16Z | 7 + manifest | 0 () |

Snapshot commit of the tracked tree (297 files): 0e11e9f1f14f. Future freezes: `python code/vcs_freeze_tags.py tag <name> <manifest>`.

### Data files: who downloaded them, when, and checksum result

| file | source | downloaded by | when (UTC) | publisher checksum | check |
|---|---|---|---|---|---|
| K562_essential_raw_bulk_01.h5ad | https://ndownloader.figshare.com/files/35773070 | agent (scripted) 2026-10-02; re-downloaded by the user (browser) 2026-10-03 | 2026-10-02T05:21:54Z | figshare md5 8321d5d3ffc99db2a5c71edca4189735 | match |
| gene2go_all.pkl | https://dataverse.harvard.edu/api/access/datafile/6153417 | agent (scripted, sandbox) | 2026-10-02T05:21:59Z | dataverse md5 77c9af0c61c30ea4d7a85680f4d122dc | match |
| essential_all_data_pert_genes.pkl | https://dataverse.harvard.edu/api/access/datafile/6934320 | agent (scripted, sandbox) | 2026-10-02T05:22:03Z | dataverse md5 b7bc2a91ca513b86f27f090d963711a6 | match |
| go_essential_all.tar.gz | https://dataverse.harvard.edu/api/access/datafile/6934319 | agent (scripted, sandbox) | 2026-10-02T05:22:14Z | dataverse md5 b8bff0d53407f26648330d264df7fe16 | match |
| 9606.protein.info.v12.0.txt.gz | https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz | agent (scripted, sandbox) | 2026-10-02T05:22:58Z | no publisher checksum (STRING provides none) | n/a; SHA-256 matches the 2026-10-01 record |
| 9606.protein.links.detailed.v12.0.txt.gz | https://stringdb-downloads.org/download/protein.links.detailed.v12.0/9606.protein.links.detailed.v12.0.txt.gz | agent (scripted, sandbox) | 2026-10-02T05:23:12Z | no publisher checksum (STRING provides none) | n/a; SHA-256 matches the 2026-10-01 record |
| linear_perturbation_prediction-Paper-publication.zip | https://zenodo.org/api/records/16092690/files/const-ae/linear_perturbation_prediction-Paper-publication.zip/content | agent (scripted, sandbox) | 2026-10-02T05:23:27Z | zenodo md5 e352192c22e6577edba3dd788677caf2 | match |
| single_perturbation_jobs_stats.tsv | https://zenodo.org/api/records/16092690/files/single_perturbation_jobs_stats.tsv/content | agent (scripted, sandbox) | 2026-10-02T05:23:28Z | zenodo md5 2b1bdd900d4ffe219ced2cce2fa64001 | match |
| cell-gears-0.1.2.tar.gz | https://files.pythonhosted.org/packages/37/34/bb7c4a418fbb3c1f7216bc1d6396b328ef9c338631c616295da7d444d6b3/cell-gears-0.1.2.tar.gz | agent (scripted, sandbox) | 2026-10-02T05:23:28Z | pypi sha256 11b8a180b0af7fa797999c0272a5f1c39f5e89bb6d09b784df5cd2b892507959 | match |
| K562_essential_raw_singlecell_01.h5ad | https://ndownloader.figshare.com/files/35773219 | user (browser download) | 2026-10-02 (download); verified 2026-10-03T07:00:56Z | figshare md5 4f1122ce1c7f13299a68df6459a266d3 | match |
| K562_gwps_raw_singlecell_01.h5ad | Figshare+ article 20029387 v1, file id 35775507, https://ndownloader.figshare.com/files/35 | user (user-run launcher outside the sandbox, code/download_resumable.py) | 2026-10-06T01:38:59Z | MD5 887e3e6a8c8df6eadf7a3030a53c9546 | match |
| rpe1_raw_bulk_01.h5ad | Figshare+ article 20029387 v1, file id 35775581, https://ndownloader.figshare.com/files/35 | agent (code/download_resumable.py) | 2026-10-06T03:15:58Z | MD5 74765fa87635467a869ea972356ae0e7 | match |
| K562_gwps_raw_bulk_01.h5ad | Figshare+ article 20029387 v1, file id 35774443, https://ndownloader.figshare.com/files/35 | agent (code/download_resumable.py) | 2026-10-06T03:17:44Z | MD5 4570b53c9d62ff6df281e622f0350060 | match |
| rpe1_raw_singlecell_01.h5ad | Figshare+ article 20029387 v1, file id 35775606, https://ndownloader.figshare.com/files/35 | agent (code/download_resumable.py) | 2026-10-06T03:33:06Z | MD5 6a2a9d0d2bf4ec147f4d1104043b268c | match |
| NormanWeissman2019_filtered.h5ad | Zenodo 13350497 (scPerturb v1.4) | user (own network; placed in G:\biology\geo_GSE264667\ instead of the planned folder) | 2026-11-01T09:29:07Z (file modification time; download time not logged) | MD5 c870e6967d91c017d9da827bab183cd6 (Zenodo API) | match (verified by agent 2026-11-03) |
| NadigOConner2024_hepg2.h5ad | Zenodo 13350497 (scPerturb v1.4) | agent (code/download_zenodo_resumable.py, background, at the user's request 2026-11-04) | started 2026-11-04T10:11Z; in progress | MD5 af2be47f7477cf32fa6e4bec1c6a4868 (Zenodo API) | pending |
| NadigOConner2024_jurkat.h5ad | Zenodo 13350497 (scPerturb v1.4) | agent (code/download_zenodo_resumable.py, background, at the user's request 2026-11-04) | started 2026-11-04T10:11Z; in progress | MD5 d8b05d00bfbd686d37ffdd4293bc6c8c (Zenodo API) | pending |
| WesselsSatija2023.h5ad | Zenodo 13350497 (scPerturb v1.4) | user (requested 2026-11-04 after the v1.4 freeze) | not yet downloaded | MD5 6897bfdcda928a678208fecf4eeb282e (Zenodo API) | pending |
| go-basic.obo | https://release.geneontology.org/2026-05-19/ontology/go-basic.obo | agent (scripted, sandbox) | 2026-10-10 (before the v0.7 freeze; exact time not logged) | no publisher checksum; SHA-256 recorded in freeze_v0.7.json | current file matches recorded SHA-256 |
| goa_human.gaf.gz | https://release.geneontology.org/2026-05-19/annotations/goa_human.gaf.gz | agent (scripted, sandbox) | 2026-10-10 (before the v0.7 freeze; exact time not logged) | no publisher checksum; SHA-256 recorded in freeze_v0.7.json | current file matches recorded SHA-256 |
| go_release-date.json | https://release.geneontology.org/2026-05-19/metadata/release-date.json | agent (scripted, sandbox) | 2026-10-10 (before the v0.7 freeze; exact time not logged) | no publisher checksum; SHA-256 recorded in freeze_v0.7.json | current file matches recorded SHA-256 |
| ReactomePathways.gmt.zip | https://reactome.org/download/current/ReactomePathways.gmt.zip | agent (scripted, sandbox) | 2026-10-10 (before the v0.7 freeze; exact time not logged) | no publisher checksum; SHA-256 recorded in freeze_v0.7.json | current file matches recorded SHA-256 |
| corum_human_current.txt | https://mips.helmholtz-muenchen.de/fastapi-corum/public/file/download_current_file?file_id=human&file_format=txt | agent (scripted, sandbox) | 2026-10-10 (before the v0.7 freeze; exact time not logged) | no publisher checksum; SHA-256 recorded in freeze_v0.7.json | current file matches recorded SHA-256 |

- Update 2026-11-04 18:47 (UTC+8): NadigOConner2024_hepg2.h5ad complete (850,590,740 bytes; MD5 af2be47f7477cf32fa6e4bec1c6a4868 verified by the download script; 8 connection timeouts recovered by resume); set read-only. Jurkat download in progress.
- Update 2026-11-05 20:25 (UTC+8): NadigOConner2024_jurkat.h5ad completed 2026-11-04 19:26 (UTC+8) by the agent's background download (1,293,665,804 bytes; MD5 d8b05d00bfbd686d37ffdd4293bc6c8c verified; 10 timeouts recovered by resume); set read-only. WesselsSatija2023.h5ad downloaded by the agent at the user's request on 2026-11-05 20:12–20:22 (UTC+8) after the v1.4 freeze (219,393,529 bytes; MD5 6897bfdcda928a678208fecf4eeb282e verified; a first attempt failed at the initial request with an SSL EOF error and was restarted; 3 interruptions recovered by resume); set read-only. No content of either file has been read.
- 2026-11-05 (6f): QC run (code/qc_nadig_wessels_v121_v14.py; no scores). Nadig HepG2 and Jurkat (v1.2.1 Q1–Q5): all pass; scPerturb files have the same shape as the GEO originals, raw integer counts, 2,393 targets, 4,976 / 12,013 control cells; ncounts equals the X row sum exactly (used as library size per v1.2.1 §5); UMI_count correlates 0.99999 with the row sum but differs by up to 2,080 / 4,178. Results in results/v1.2/{hepg2,jurkat}/qc_v121.json. A5 awaits user confirmation. Wessels (v1.4 §3): identity confirmed (30,707 cells); G1 136 evaluable pairs (pass); G2 136/157 = 86.6% (FAIL; MYB single 15 cells, MED13L single 27 cells, affecting 21 pairs); G3 424 control cells (FAIL); G4 pass; G5 24 genes (pass). Per v1.4 §2.4 item 4 the conclusions are restricted to Norman (no loosening); Wessels is not analysed. Manuscript: Discussion sentence 'its analysis … is pending' replaced with the gate outcome; Table S13 gains the data-level row (to-do item 4 resolved for the Discussion; presubmission placeholder still open).

## 2026-11-05 (6g): A5 — HepG2 and Jurkat (protocol v1.2 §3.4 + v1.2.1)

- Code: code/nadig_v12.py (new file, create-only). Pseudobulk in one backed pass (row blocks of 2,000 cells; library size obs ncounts per v1.2.1 §5; half seed master_seed+15000+{1,2}); unit table by data_v07.build_units rules with obs-aggregated metadata (the bulk reader replaced by a function returning the single-cell obs/readout names; no rule changed); splits repeated_group_kfold, seed base master_seed+16000+{1,2}, written read-only; per-split linear-PCA as pipeline_v05/eval_v06 (K from {5,10,20,30,50} on validation, lambda 0.1); statistics, controls and judge() imported from norman_v12. Clusters for the bootstrap: TSS groups. Peak process memory 1.96 GB (C0), 1.9 GB (run).
- C0: HepG2 E 1632 (1569 TSS groups), Jurkat E 1908 (1835); half reconstruction error < 6e-15.
- Result (pre-specified prediction: R not elevated): HepG2 R 0.086 (0.07 to 0.11), Jurkat R 0.038 (0.03 to 0.05); both 'R not elevated' — prediction holds. Both negative controls pass; split-half control R 0.943 / 0.702 (threshold reachable). T5 coverage linear-PCA minus mean: 0.1 pp (−3.9 to 4.2) / −0.6 pp (−3.2 to 2.0).
- Manuscript: new Results subsection, Methods paragraph, Abstract clause, Discussion sentence (replacing 'had not been processed'), limitation 12 (same lab/platform; ncounts vs UMI_count library size, median ratio 1.01; K at grid top), new Supplementary Table S11 (former S11–S13 renumbered S12–S14 by first citation), Fig. 2 redrawn with HepG2/Jurkat and caption updated.

## 2026-11-08: manuscript v1.2

- Gemini language edit of v1.1 reviewed (polish_review_2026-11-06.md); v1.2 = Gemini text with format restored and 59 passages reverted to evidence-matched wording (manuscript_v1.2_changes.md). Author line and AI-use paragraph updated. No number, citation or reference changed. Supplementary file unchanged (still v1.1).

## 2026-11-08: manuscript v1.4 (readability edit)

- Style-only rewrite of 27 prose paragraphs in Background, Results, Discussion and Conclusions (long sentences split, fewer em-dashes and stock connectors); Abstract, Methods, Declarations and References unchanged. Each paragraph passed automatic checks (numbers, citations, figure/table references, code spans, bold runs, qualifier and negation counts, no new banned words) and two meaning audits; one paragraph kept unchanged after an audit flag, and 9 sentences corrected by hand for hedge scope. Files: manuscript_v1.4_GenomeBiology.md, manuscript_v1.4.md.
