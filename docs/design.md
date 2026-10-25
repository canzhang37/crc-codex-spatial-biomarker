# Design: how robust is a spatial immune biomarker in colorectal cancer?

This plan was written before the single-cell table or the patient table was downloaded. Only the paper (Schürch et al., Cell 2020; PMC7479520), its STAR Methods and the file listing of the Mendeley record were read. Changes made after the files are opened go in the Amendments section, with a reason.

## Background

Schürch et al. imaged the tumor invasive front of 35 advanced colorectal cancer patients with CODEX (56 protein markers, 140 tissue-microarray spots, 4 per patient, 258,385 cells in 29 cell clusters, one of them imaging artifacts). The patients were chosen as 17 with a Crohn's-like reaction (CLR, many tertiary lymphoid structures) and 18 with diffuse inflammatory infiltration (DII, none). Each cell's window of its 10 nearest neighbors was described by its 29 cell-type frequencies, and the windows were clustered with MiniBatchKMeans (k = 10) into cellular neighborhoods (CNs); one CN of artifacts was dropped, leaving nine.

The paper's survival claim: among the 18 DII patients (13 deaths), the log frequency (pseudocount 0.001) of PD-1+ CD4+ T cells within CN-9, the granulocyte-enriched CN, was associated with longer overall survival (Cox model, p = 0.006). Neither the overall frequency of these cells nor the amount of CN-9 alone was. The feature was chosen from CN-specific frequencies ranked by a classifier of CLR against DII, and each candidate was tested in its own Cox model. The Kaplan-Meier figure splits the patients at the best cut point.

The claim rests on 18 patients and 13 deaths, a feature chosen among many, patient frequencies estimated from four spots each, and a best-cut-point curve. This project checks how much of it holds once each of these is taken into account. It does not question the imaging or the cell calls.

## Data

- Single-cell table: Schürch C. (2020), Mendeley Data, doi:10.17632/mpjzbtfgfr.1, CC BY 4.0, file `CRC_clusters_neighborhoods_markers.csv` (222,999,273 bytes, SHA-256 416cc3926a7a900ce3b22a33be699535ca35ae85fca4585a8ba5dad6d7f3c677). One row per cell with coordinates, cell cluster, the authors' CN labels, marker intensities, patient, spot and group.
- Patient table: the paper's supplementary spreadsheet `mmc2.xlsx` (Table S1), from the Europe PMC supplementary files of PMC7479520, expected to hold overall survival, its censoring indicator and disease-free survival.
- Column names are taken from a public tutorial that uses the same files and are checked against the files on first opening (check 1 below).

## Checks before any analysis

1. Columns, types and missing values in both files; number of cells, spots and patients; cells per spot; patients per group.
2. Agreement with the paper: 258,385 cells, 35 patients (17 CLR, 18 DII), 140 spots, 29 clusters, ten CN labels with one of artifacts, 13 deaths among the DII patients.
3. How T-cell functional positivity (ICOS, Ki-67, PD-1) is recorded. If the table holds the authors' positive or negative calls, they are used. If it holds only intensities, positivity is set per marker by Otsu's threshold on arcsinh-transformed intensities within CD4+ T cells, CD8+ T cells and Tregs pooled, and the share of PD-1+ CD4+ T cells in CN-9 is compared with the paper's figures before going on.

These checks look at counts, codes and labels only, not at survival.

## Part A: replicate

1. Rebuild the windows and CNs as described (10 nearest neighbors within a spot, including the cell itself; 29 frequencies; MiniBatchKMeans with k = 10) with 20 random seeds. Agreement with the authors' labels by adjusted Rand index, and the best one-to-one matching of clusters to CNs by cell counts.
2. With the authors' CN labels, recompute the patient feature: PD-1+ CD4+ T cells in CN-9 as a proportion of all cells in CN-9 for that patient, pooled over the patient's spots. (To be checked when the data are opened: the denominator the paper used, all CN-9 cells or all CD4+ T cells in CN-9.) Fit the Cox model of overall survival on log(0.001 + frequency) in the 18 DII patients and report the hazard ratio, Wald and likelihood-ratio p values.

Replication target: the same direction and a likelihood-ratio p within a factor of two of 0.006. If it is not met, the reason is found and written up before Part B starts.

## Part B: robustness of the survival association

All in the 18 DII patients, with the authors' CN labels.

1. Multiplicity. The family is every CN-specific frequency of the nine functional T-cell subsets (ICOS+, Ki-67+ and PD-1+ within CD4+ T cells, CD8+ T cells and Tregs) in the eight CNs that are neither artifacts nor the follicle CN: 72 Cox models. Adjusted p values by Westfall and Young's max-T permutation, permuting the (time, status) pairs across patients 10,000 times, and by Holm's method.
2. Spot sampling. Each patient's frequency comes from four spots. The spots are resampled with replacement within each patient (2,000 times), the feature is recomputed and the Cox model refitted. Reported: the share of resamples with p < 0.05 and the 2.5th and 97.5th percentiles of the log hazard ratio. The intraclass correlation of the spot-level log frequency gives its reliability.
3. Influence. The Cox model is refitted 18 times, each time without one patient.
4. Pseudocount and scale. Pseudocounts of 0.0001 and 0.01, and the rank of the frequency instead of its log.
5. Cut point. The best-split log-rank test is corrected with the maximally selected rank statistic (Lausen and Schumacher), the usual correction for searching over cut points.
6. Proportional hazards. Scaled Schoenfeld residuals against time.

Predictions, fixed now: Part A replicates; the Westfall-Young adjusted p value is above 0.05; fewer than 80% of spot resamples give p < 0.05; at least one leave-one-out fit has p above 0.05. Each prediction is reported as met or not met. A failed prediction is a result, not a reason to change the analysis.

## Part C: what the spatial features add beyond cell composition

Question: does knowing where cells sit, beyond how many there are, separate CLR from DII patients? The groups were selected by the presence of tertiary lymphoid structures, so the follicle CN, B cells and follicular dendritic cells are left out of every feature set, or the task would be solved by construction.

- Set 1, composition: the patient's proportions of the remaining cell types, centered log-ratio transformed (pseudocount 0.5 cells).
- Set 2, set 1 plus CN-specific frequencies: each cell type's proportion within each remaining CN.
- Set 3, set 1 plus spatial contacts: for each pair of cell types, the log ratio of observed to expected contacts on the Delaunay graph of each spot, with the expectation from shuffling labels within the spot (200 shuffles), averaged over the patient's spots.
- Model: L1-penalized logistic regression on standardized features, penalty chosen by inner stratified five-fold cross-validation (log loss), scored by outer leave-one-patient-out cross-validation (35 fits).
- Measures: area under the ROC curve with a bootstrap interval over patients, and a label-permutation p value (1,000 permutations, each rerunning the whole nested cross-validation) for each set. If one set's permutations take more than 30 minutes on the laptop, they run on MSI. Set 2 or 3 is said to add to set 1 only if its AUC exceeds set 1's by at least 0.05 and the paired bootstrap interval of the difference lies above zero. With 35 patients the intervals will be wide, and the readout says so.

## Optional Part D: a graph neural network

If time allows after Parts A to C: a two-layer graph convolutional network (PyTorch Geometric) on each spot's cell graph, nodes carrying one-hot cell types, pooled to the spot and averaged to the patient, scored by the same leave-one-patient-out scheme as Part C. It is exploratory and reported beside Part C, not instead of it.

## Outputs

`results/checks/`, `results/replication.csv`, `results/robustness_*.csv`, `results/composition_vs_spatial.csv`, figures in `figures/`, and a readout in `docs/readout.md`. Figures show patients only through model outputs and spot images of cell positions, never clinical details beyond group and survival.

## Amendments

### Amendment 1 (after the checks, before any survival or classification result)

The checks (`results/checks/`) match the paper: 258,385 cells, 35 patients (17 CLR coded 1, 18 DII coded 2), 140 spots with four per patient (84 to 3,887 cells per spot), 29 cell clusters (one called "dirt"), ten CN labels in `neighborhood10`. `neighborhood number final` holds the paper's CN-1 to CN-9 (CN-5 follicle, CN-9 granulocyte enriched) and is empty for the 4,566 cells of the artifact CN. In Table S1 (sheet A) `OS` is in months and `OS_Censor` = 1 marks a death: 13 deaths among the 18 DII patients and 4 among the 17 CLR patients.

Decisions that follow, none of them informed by a survival result:

1. Positivity (check 3). The table holds the authors' 0/1 calls for the nine T-cell subsets (`CD4+ICOS+`, `CD4+Ki67+`, `CD4+PD-1+`, `CD8+...`, `Treg-...`) and six macrophage subsets. They are used as given. The calls are marker gates, not cluster restricted: `CD4+PD-1+` is 1 for 2,485 cells, 1,112 of them in the cluster "CD4+ T cells CD45RO+" and others in macrophage, B-cell and granulocyte clusters. As a sensitivity analysis in Part B, the call is counted only in the CD4+ T-cell clusters ("CD4+ T cells", "CD4+ T cells CD45RO+", "CD4+ T cells GATA3+").
2. Denominator (Part A2). The authors' code (nolanlab/NeighborhoodCoordination, `cca_cleaned_up.ipynb`) computes `np.log(1e-3 + cells.groupby(['patients','neighborhood10'])[subsets].mean())`, so the frequency is the share of all of the patient's cells in that CN carrying the call, with a natural log. This is the main definition. The share among the patient's CD4+ T-cell-cluster cells in CN-9 is a further sensitivity analysis.
3. Zero counts. All 18 DII patients have cells in CN-9 (60 to several hundred); 4 of them have no `CD4+PD-1+` cell there, so their value is log(0.001). A patient without any cell in a CN gets frequency 0 for that CN in the 72-feature family, as in a pseudocount analysis.
4. Part C cell types. The clusters "dirt", "undefined", "tumor cells / immune cells" and "immune cells / vasculature" are left out with "B cells", as the authors also removed undefined and mixed clusters before their contact analysis. Cells in CN-5 and in the artifact CN are left out of every Part C feature.

### Amendment 2 (post hoc, after Parts A and B)

Part A1 showed that the rebuilt CNs agree with the published labels only moderately and vary with the k-means seed (adjusted Rand index 0.52 to 0.80 over 20 seeds, median 0.60). This was not foreseen, so the following check is post hoc and reported as such: for each seed, the rebuilt cluster matched to the published CN-9 replaces it, and the Part A2 Cox model is refitted. It shows how much the survival result depends on one run of the clustering.

A second post hoc check comes from the spot figure of Part B: several DII patients with the lowest values died within weeks of surgery, which may reflect postoperative deaths rather than tumor progression. The Part A2 model is refitted without patients whose overall survival is under 3 months. Neither check changes the pre-registered verdicts.
