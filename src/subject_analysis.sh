#!/bin/bash

Rscript mage_processing.R ../configs/mage_processing.yaml

Rscript prepare_counts_for_subject_analysis.R ../configs/prepare_counts_for_subject_analysis.yaml

Rscript export_subject_genotypes.R ../configs/export_subject_genotypes.yaml

Rscript subject_analysis.R ../configs/subject_analysis.yaml
