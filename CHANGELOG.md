# Changelog

All notable changes to the **gms-nallo-methylation-dashboard** project will be documented in this file.

## [0.3.0] - 2026-09-11

### Added
- **Cohort Reference Benchmark Visualization**:
  - Integrated interactive horizontal box-whisker-scatter plot powered by Chart.js comparing patient sample bulk methylation against population healthy controls (Epigenome of Sweden / GoS, 288 controls).
  - Highlights normative imprinting reference envelopes: 95% normal reference bar (`#D1E8FF`), Interquartile Range (IQR Q1–Q3 bar, `#4A90E2`), cohort median marker (`#003366`), and patient sample methylation floating point (`#FFCC00` Safety Gold).
  - Built-in visual boundary shading for normative ASM window (`[0.35, 0.65]`), hypomethylation threshold (`<0.35`), and hypermethylation threshold (`>0.65`).
  - Added "All Loci" vs. "Selected Only" toggle filter to isolate the specific locus selected in the clinical triage grid.
  - Interactive cross-linking: clicking any bar in the benchmark chart selects and scrolls to the corresponding row in the clinical grid.
  - One-click **Export Plot** to high-resolution PNG.
- **Collapsible Reference Statistics Table**:
  - Detailed locus metrics table reporting control locus tag, cohort sample size ($N$), median, IQR ($Q_1 - Q_3$), 95% normal interval ($[\text{CI}_{\text{lo}}, \text{CI}_{\text{hi}}]$), sample score, and clinical status badge.
  - Clicking any table row synchronizes and highlights the row in the main Ag-Grid.
- **Dynamic Reference Dataset Upload (Zero-Backend Offline Support)**:
  - Added **"Upload Reference"** button supporting custom reference files (`.tsv`, `.csv`, `.txt`, `.json`) parsed entirely in-browser via the FileReader API.
  - Instant re-rendering of benchmark plot and reference table upon file selection without network dependencies.
  - Added **"Download Template"** button to export reference dataset and column headers as a ready-to-edit TSV template.
  - Added **"Reset"** button to revert back to built-in GoS (288 controls) baseline with one click.
  - Dynamic status badge reflecting active reference source and locus count.
- **Bi-directional CRUD Reactivity**:
  - Deleting, adding, or modifying loci in the Ag-Grid live-updates the GoS benchmark chart, sample distribution charts, global ASM pie chart, and reference table in real time.

### Changed
- **ASM Classification Scale Adjustment**:
  - Updated Allele-Specific Methylation (ASM) normative window to **35% – 65%** (hypomethylated $<35\%$, hypermethylated $>65\%$) based on empirical distribution analysis of the 288 healthy GoS control cohort.
- **Clinical Color Scheme**:
  - Transitioned chart elements to high-contrast clinical palette with focused prominence on user-selected locus and dimmed opacity on unselected loci.
  - Increased Y-axis locus font size to 13px bold for enhanced legibility on high-resolution clinical monitors.

## [0.2.4] - 2026-06-25

### Added
- **Multi-profile Scanning Support**:
  - The `nallo_methylation_manager.py` scanning mechanism has been upgraded to discover and process *all* matching `.csv`/`.tsv` files (containing `'methbat'`, `'profile'`, or `'vs'`) within each sample directory instead of stopping at the first candidate.
  - Outputs separate, descriptively-named report files based on the input file's stem (e.g., `BWS-1_methbat_report_GOS_methylation_report.html`), allowing easy comparison of the same sample against different reference datasets.

## [0.2.3] - 2026-06-24

### Added
- **Dedicated External Links Column**:
  - Integrated quick-access badges for genomic browsers:
    - **G**: Link directly to the region on **gnomAD v4 (GRCh38)**.
    - **N**: Link directly to the region on the **NCBI Genome Data Viewer (GRCh38)**.
  - Smart string parsing to strip `"chr"` prefixes for compatibility with NCBI coordinate URLs.
- **Dynamic Columns Toggle Dropdown**:
  - Added a "Columns" selector button to show/hide any data column dynamically in the Ag-Grid table.
- **Drag-and-Drop Column Re-ordering**:
  - Enabled header-dragging to allow users to completely re-arrange the table layout to their preference.
- **Collapsible IGV-Web Instructions**:
  - Moved the step-by-step instructions for loading phased BAMs and setting up 5mC color modes into a collapsible `<details>` panel directly above the viewer, preserving vertical screen space while remaining accessible when the sidebar is hidden.
- **Interactive IGV Panel Trigger**:
  - Hidden the IGV-Web viewer by default to boost initial page performance, replacing it with an on-demand **"Launch IGV"** toggle button.

## [0.2.2] - 2026-06-18

### Added
- **Interactive Table Editing (CRUD)**:
  - Enabled editing cells directly within the Ag-Grid interface.
  - Added **"Add Locus"** button to dynamically insert new genomic regions.
  - Added a **Trash / Delete** button next to each row for quick removal of loci.
  - Recalculates the genomic `Interval` automatically when `Start` or `End` coordinates are edited.
- **Separate HTML Reports**:
  - Refactored `nallo_methylation_manager.py` to output a unique, self-contained HTML report per sample (e.g., `{sample_id}_methylation_report.html`) instead of a single multi-sample document.
- **New Columns**:
  - Added `End Position` and `Expected Category` to the main report table.
- **Chart.js Visualizations**:
  - Integrated interactive vertical bar charts showing distribution tallies (ASM vs. Methylated vs. Unmethylated) per selected sample.
  - Added a project-wide pie chart summarizing total ASM counts in the sidebar.
