# BGI Amplicon Pipeline Refinement Log & Edge Cases
*Date: April 25, 2026*

## Implemented Modifications

1. **`10_venn_flower.R` Integration**:
   - Deployed a completely decoupled base R layout handler implementing purely parametric custom Flower Plots matching strict BGI geometries without secondary libraries.
   - Handled bifurcation dynamically to force edge layouts onto `N=2` comparisons protecting parallel dual box topologies natively bridging away from radical ellipses.
   - Forced algorithm dependencies to govern aspect ratio (`semi-minor = max(0.25, 2.2 / n)`) to keep intersections strictly uncluttered during massive 17-group analysis overlays.
   
2. **`12_plsda.R` Structural Integrities**:
   - Scrubbed runtime `BiocManager::install("mixOmics")` dependencies strictly adhering to pre-installs under `install_packages.R` limiting internet constraints on high-performance node blocks.
   - Validated standard `mixOmics::perf` models to compute cross-validation explicitly via M-Folds capturing error validations (`BER`) output. Implemented dynamic fold logic guaranteeing iterations bypass crashing boundaries automatically when comparative structures fail minimum thresholds (`< 2`).
   
3. **`15_multilevel_taxa.R` Namespace Safeties**:
   - Overhauled array mapping inside `aggregate()` enforcing index constraints by natively verifying name mappings explicitly against target column variants natively dodging blind indexed extractions `[, -1]`.
   
4. **`03` & `15` Outputs File Tracking**:
   - Rebound native `ggsave()` and `pheatmap()` commands internally generating `.comp_prefix` bounds automatically locking files individually rather than letting independent runs overwrite base names during standalone multi-group sweeps matching native BGI nested architecture.

## Encountered Edge-Cases

* **Taxonomic Disambiguation Crash (`15_multilevel_taxa.R`)**: During analysis over sparse relative subsets, slicing variables by final nested taxa groupings (`strip_taxa_prefix()`) consolidated multiple distinctly isolated biological arrays strictly back to `"Other"` identically causing R to trigger fatal `duplicate 'row.names' are not allowed` during matrix constructions into `pheatmap`. 
  - *Fix Deployed*: Processed the `heatmap_labels` array natively through `make.unique()` guaranteeing arbitrary unique assignments (e.g. `Other.1`, `Other.2`) resolving matrix constraints correctly without losing grouping boundaries. (Future iteration consideration: Update this to a fully parent-qualified disambiguation function (e.g. `Other (Euryarchaeota)`) for further scientific clarity natively avoiding generic enumeration labels).
