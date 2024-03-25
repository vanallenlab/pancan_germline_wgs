#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Shell subroutine to stage sample lists and manifests in an AoU workspace

# Note that this code is designed to be run inside the AoU Researcher Workbench
# It also depends on variables set in aou_gatk.sh

# Copy sample info to local disk
gsutil -m cp -r $MAIN_WORKSPACE_BUCKET/data/sample_info/sample_lists ./

# Copy manifest of CRAM/CRAI paths
gsutil -u $GPROJECT -m cp \
  gs://fc-aou-datasets-controlled/v7/wgs/cram/manifest.csv \
  ./cram_manifest.csv

# Write list of cancers ready to be processed
cat << EOF > cancers.list
pancreas
esophagus
stomach
lung
liver
colorectal
kidney
melanoma
prostate
ovary
breast
uterus
cns
bladder
oral
other
male_controls.p1
male_controls.p2
female_controls.p1
female_controls.p2
multiple_sites
pediatric
ufc_cases
ufc_controls
EOF
# Note: as of 1/8/24, single-sample processing for cancers were divided among 
# multiple workspaces as follows:
# Main workspace: pancreas, esophagus, stomach, lung, liver, kidney, ovary, cns, male_controls.p1, pediatric
# Second workspace: colorectal, bladder, male_controls.p2, multiple_sites, ufc_cases
# Third workspace: melanoma, uterus, female_controls.p1, ufc_controls
# Fourth workspace: prostate, oral, female_controls.p2
# Fifth workspace: breast, other (TGCT)

# Make .tsv mapping person_id, cram path, and crai path for each cancer type
while read cancer; do
  fgrep -wf sample_lists/$cancer.samples.list cram_manifest.csv \
  | sed 's/,/\t/g' | sort -Vk1,1 \
  > data/cram_paths/$cancer.cram_paths.tsv
done < cancers.list
