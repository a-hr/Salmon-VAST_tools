#!/usr/bin/env bash

vast_tools_full_table=$1
min_samples=$2  # 70% muestras

perl Get_Event_Stats.pl $vast_tools_full_table \
	--N $min_samples \
	--VLOW \
	--outfile coverage_filtered_$(basename $vast_tools_full_table)
