# Summary of Salmon and VAST-TOOLS based workflows

> Álvaro Herrero Reiriz
> 2024-05-10

This document provides a summary on how to use the Salmon and VAST-TOOLS based workflows. The folders `salmon` and `vast-tools` contain the scripts and instructions to run the workflows.

- [Summary of Salmon and VAST-TOOLS based workflows](#summary-of-salmon-and-vast-tools-based-workflows)
  - [Environment setup](#environment-setup)
  - [Differential expression analysis with Salmon](#differential-expression-analysis-with-salmon)
  - [Alternative splicing analysis with VAST-TOOLS](#alternative-splicing-analysis-with-vast-tools)
    - [Previous steps](#previous-steps)
    - [Table processing](#table-processing)

## Environment setup

Create two conda environments with the following commands and clone the environments from the provided yml files:

```bash
conda env create -p /scratch/heral/envs/salmon_env -f salmon_env.yml
conda env create -p /scratch/heral/envs/vast_env -f vast_env.yml
```

## Differential expression analysis with Salmon

1. Once the output from Salmon is obtained (quant.sf files), the first thing to do is to generate the TPM matrices. This can be done with the `1_create_table.sh` script. We will generate a TPM matrix where each of the columns corresponds to a sample and each of the rows corresponds to a gene.

    The requited input arguments are:
    - `outputDir`: Directory containing the output from Salmon. 
    - `metric`: Metric used to quantify the expression. Default is lengthScaledTPM. For more information, see the `tximport` package documentation.
    - `designFile`: File containing the design of the experiment (i.e. the samples and the groups they belong to).
    - `gtf`: Path to the GTF file, used to assign feature names.
    - `outname`: Name prefix for the output file.

    After modifying the script with the required arguments, run it with the following command:

    ```bash
    sbatch 1_create_table.sh
    ```

2. Once the TPM matrix is generated, we can perform the differential expression analysis with the `2_DiffExpressionAnalysis.R` script. **Recommended to run in local.** The script requires the following arguments, which can be modified in the script itself:
   
    - `project_dir`: Directory where the output files will be saved.
    - `gtf_file`: Path to the GTF file.
    - `design_tab`: Path to the design file.
    - `tpms_path`: Path to the TPM matrix. It can be gzipped.

    - `group_a`: Reference group for the analysis.
    - `group_b`: Group to compare with the reference group.

    - `nreads`: Minimum number of reads to consider a gene for the analysis. Default is 10.
    - `fold_change`: Minimum log2 fold change to consider a gene as differentially expressed. Default is 2.

    After modifying the script with the required arguments, it can be run by clicking the `Source` button in RStudio.

    > Make sure to provide all of the paths as ABSOLUTE paths.
    
3. All the output files will be saved in the `project_dir/$metric_$group_b_vs_$group_a` directory.

## Alternative splicing analysis with VAST-TOOLS

### Previous steps

Starting from the `FASTQ` files, run `vast_align.sh` and then `vast_combine.sh` as explained in the VAST-TOOLS documentation.

> The output of these scripts will be a `INCLUSION_LEVELS_FULL.tab` file, the starting point for the following steps.

### Table processing

1. Coverage filter: first off, we need to filter out events with low coverage. This can be done with the `1_filter_coverage.sh` script. The script requires the following arguments:

    - `vast_tools_full_table`: Path to the full output table.
    - `min_samples`: Integer, equals the 70% of the sample count. Sets the minimum number of samples where an event must be present with the VLOW coverage.

    After modifying the script with the required arguments, run it with the following command:

    ```bash
    bash 1_filter_coverage.sh
    ```

2. Variance filter (optional): if desired, we can filter out events with low variance. This can be done with the `2_filter_variance.R` script. The script requires the following arguments:

    - `dir`: Path to the project directory.
    - `coverage_filtered_table`: Path to the output table from the previous step.
    - `full_table`: Path to the full VAST-TOOLS output table.
    - `output_path`: Path to store the variance filtered table.

    After modifying the script with the required arguments, run it from RStudio by clicking the `Source` button.

    > Make sure to activate the `salmon_env` environment in case you run this on the cluster. In local, just make sure you've installed `tidyverse` and `matrixStats` packages.

3. Table and plot generation: finally, we can generate the table and plots with the `3_analysis.sh` script. The script requires the following arguments:

    - `input_table_filt`: Path to the coverage (and variance if done) filtered table.
    - `input_table_full`: Path to the full VAST-TOOLS output table.
    - `design_tab`: Path to the file containing sample and group names. Expected columns are: `SampleName` and `GroupName`.
    - `group_A`: Exact name of group A (reference group) in `design.tab`.
    - `group_B`: Exact name of group B (contrast group) in `design.tab`.
    - `min_comparison_pct`: Minimum percentage of sample comparisons where minimum filtering conditions must be met in order to include an event.
    - `min_dpsi`: Minimum delta PSI to consider a comparison relevant.

    After modifying the script with the required arguments, run it with the following command:

    ```bash
    bash 3_analysis.sh
    ```

    > Make sure to activate the `vast_env` environment before running the script.

