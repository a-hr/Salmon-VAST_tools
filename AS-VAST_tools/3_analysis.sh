#!/usr/bin/env bash -l

# ---- INPUTS ----

input_table_filt="/Users/varo/Library/Mobile Documents/com~apple~CloudDocs/Biodonostia/GB/AS_processing/data/coverage_variance_filtered_INCLUSION_LEVELS_CGGA.tab.gz"
input_table_full="/Users/varo/Library/Mobile Documents/com~apple~CloudDocs/Biodonostia/GB/AS_processing/data/INCLUSION_LEVELS_FULL-hg38-121_CGGA.tab.gz"
design_tab="/Users/varo/Library/Mobile Documents/com~apple~CloudDocs/Biodonostia/GB/AS_processing/designs/design_CGGA.tab"

group_A=Oligodendroglioma_G2
group_B=Glioblastoma

output_dir="/Users/varo/Library/Mobile Documents/com~apple~CloudDocs/Biodonostia/GB/AS_processing/perl_R_python"

min_comparison_pct=50
min_dpsi=10

# ---- SETUP ----

# to create AS-env environment:
# conda create -n AS-env
# conda install -c conda-forge pandas numpy plotly python-kaleido matplotlib scipy

conda activate AS-env

# ---- RUN ----
python AS_analysis.py "${input_table_filt}" "${input_table_full}" "${design_tab}" "${group_A}" "${group_B}" "${output_dir}" "${min_comparison_pct}" "${min_dpsi}"

# ---- COMPRESS OUTPUT ----
for tab in `ls "${output_dir}/${group_B}_vs_${group_A}/tables/"`; do
  gzip "${output_dir}/${group_B}_vs_${group_A}/tables/${tab}"
done
