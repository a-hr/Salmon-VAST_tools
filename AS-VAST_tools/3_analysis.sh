#!/usr/bin/env bash -l

project_dir=/scratch/heral/Salmon-VAST_tools/AS-VAST_tools

# ---- INPUTS ----
cd $project_dir

input_table_filt="$project_dir/data/coverage_filtered_INCLUSION_LEVELS_BIOD.tab.gz"
input_table_full="$project_dir/data/INCLUSION_LEVELS_FULL-hg38-78_BIOD.tab.gz"
design_tab="$project_dir/designs/design_BIOD.tab"
group_A=Normal_brain
group_B=GB_tumor

output_dir="$project_dir/perl_python"

min_comparison_pct=50
min_dpsi=10

# ---- SETUP ----

# to create AS-env environment:
# conda create -n AS-env
# conda install -c conda-forge pandas numpy plotly python-kaleido matplotlib scipy

# conda activate /scratch/$USER/envs/AS_env

# ---- RUN ----
python AS_analysis.py "${input_table_filt}" "${input_table_full}" "${design_tab}" "${group_A}" "${group_B}" "${output_dir}" "${min_comparison_pct}" "${min_dpsi}"

# ---- COMPRESS OUTPUT ----
for tab in `ls "${output_dir}/${group_B}_vs_${group_A}/tables/"`; do
  gzip "${output_dir}/${group_B}_vs_${group_A}/tables/${tab}"
done
