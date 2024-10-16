#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Shell code to curate selected sample phenotypes for DFCI G2C phase 1 cohort

# Note that this code is designed to be run on RLC's local machine
# and references local paths that will not generalize to other systems

# This code is intended for record-keeping/reproducibility purposes
# rather than reapplication by other users


#########
# SETUP #
#########

# Set up local environment
export GPROJECT="vanallen-pancan-germline-wgs"
export WRKDIR=`mktemp -d`
export BASEDIR="/Users/ryan/Desktop/Collins/VanAllen/pancancer_wgs"
export CODEDIR="$BASEDIR/pancan_germline_wgs"


###########
# NCI GDC #
###########

# Apollo
export APOLLODIR="$BASEDIR/data_and_cohorts/apollo/clinical.cases_selection.2024-03-26"
$CODEDIR/scripts/sample_info/phenotypes/curate_nci_gdc_phenotypes.R \
  --clinical-tsv $APOLLODIR/clinical.tsv \
  --exposure-tsv $APOLLODIR/exposure.tsv \
  --biospecimen-tsv $APOLLODIR/../biospecimen.cases_selection.2024-03-26/sample.tsv \
  --cohort apollo \
  --out-tsv $WRKDIR/apollo.phenos.tsv

# CPTAC
export CPTACDIR="$BASEDIR/data_and_cohorts/cptac/clinical.cases_selection.2024-03-27"
$CODEDIR/scripts/sample_info/phenotypes/curate_nci_gdc_phenotypes.R \
  --clinical-tsv $CPTACDIR/clinical.tsv \
  --exposure-tsv $CPTACDIR/exposure.tsv \
  --biospecimen-tsv $CPTACDIR/../biospecimen.project-cptac-3.2024-08-20/sample.tsv \
  --cohort cptac \
  --out-tsv $WRKDIR/cptac.phenos.tsv

# Eagle
export EAGLEDIR="$BASEDIR/data_and_cohorts/eagle/clinical.cases_selection.2024-03-29"
$CODEDIR/scripts/sample_info/phenotypes/curate_nci_gdc_phenotypes.R \
  --clinical-tsv $EAGLEDIR/clinical.tsv \
  --exposure-tsv $EAGLEDIR/exposure.tsv \
  --biospecimen-tsv $EAGLEDIR/../biospecimen.cases_selection.2024-03-29/sample.tsv \
  --cohort eagle \
  --out-tsv $WRKDIR/eagle.phenos.tsv

# HCMI
export HCMIDIR="$BASEDIR/data_and_cohorts/hcmi/clinical.cases_selection.2024-03-26"
$CODEDIR/scripts/sample_info/phenotypes/curate_nci_gdc_phenotypes.R \
  --clinical-tsv $HCMIDIR/clinical.tsv \
  --exposure-tsv $HCMIDIR/exposure.tsv \
  --biospecimen-tsv $HCMIDIR/../biospecimen.cases_selection.2024-03-26/sample.tsv \
  --cohort hcmi \
  --out-tsv $WRKDIR/hcmi.phenos.tsv

# WCDT
export WCDTDIR="$BASEDIR/data_and_cohorts/wcdt/clinical.cases_selection.2024-03-26"
$CODEDIR/scripts/sample_info/phenotypes/curate_nci_gdc_phenotypes.R \
  --clinical-tsv $WCDTDIR/clinical.tsv \
  --exposure-tsv $WCDTDIR/exposure.tsv \
  --biospecimen-tsv $WCDTDIR/..//biospecimen.cases_selection.2024-03-26/sample.tsv \
  --cohort wcdt \
  --out-tsv $WRKDIR/wcdt.phenos.tsv



#########
# dbGaP #
#########

# BioMe
$CODEDIR/scripts/sample_info/phenotypes/curate_biome_phenotypes.R \
  --in-tsv $BASEDIR/data_and_cohorts/topmed/biome/phs001644.v2.pht009946.v2.p2.c1.TOPMed_CCDG_BioME_Subject_Phenotypes.HMB-NPU.txt \
  --out-tsv $WRKDIR/biome.phenos.tsv

# COPDGene
# Note: working on getting complete phenotype data access via dbGaP
# For now, just using the minimal data available from TOPMed accession
# We can udpate this later if/when we get access to more phenotype data
$CODEDIR/scripts/sample_info/phenotypes/curate_copdgene_phenotypes.R \
  --sample-attributes $BASEDIR/data_and_cohorts/copdgene/dbgap_downloads/phs000951.v5.pht005052.v5.p5.c1.TOPMed_WGS_COPDGene_Sample_Attributes.HMB.txt.gz \
  --topmed-subject-tsv $BASEDIR/data_and_cohorts/copdgene/dbgap_downloads/phs000951.v5.pht005050.v5.p5.TOPMed_WGS_COPDGene_Subject.MULTI.txt.gz \
  --id-link-tsv $BASEDIR/data_and_cohorts/copdgene/dbgap_downloads/phs000951.v5.pht005051.v5.p5.TOPMed_WGS_COPDGene_Sample.MULTI.txt.gz \
  --out-tsv $WRKDIR/copdgene.phenos.tsv

# CEPH
# Note: working on getting complete phenotype data access for CEPH
# For now, just using the minimal data available online
# We can udpate this later if/when we get access to detailed health history info
export CEPHDIR=$BASEDIR/data_and_cohorts/ceph/phenotypes
gsutil -m cp gs://dfci-g2c-inputs/sample-lists/ceph.samples.list $WRKDIR/
$CODEDIR/scripts/sample_info/phenotypes/curate_ceph_phenotypes.R \
  --sample-ids $WRKDIR/ceph.samples.list \
  --pedigree-tsv $CEPHDIR/phs001872.v1.pht009364.v1.p1.CEPH_Utah_Pedigree.MULTI.txt.gz \
  --subject-info-tsv $CEPHDIR/phs001872.v1.pht009363.v1.p1.CEPH_Utah_Subject.MULTI.txt.gz \
  --out-tsv $WRKDIR/ceph.phenos.tsv

# GTEx
export GTEXDIR=$BASEDIR/data_and_cohorts/gtex
$CODEDIR/scripts/sample_info/phenotypes/curate_gtex_phenotypes.R \
  --phenotypes-tsv $GTEXDIR/dbgap_phenotypes/phs000424.v9.pht002742.v9.p2.c1.GTEx_Subject_Phenotypes.GRU.txt.gz \
  --sample-attributes-tsv $GTEXDIR/GTEx_Analysis_2021-02-11_v9_Annotations_SampleAttributesDS.txt \
  --out-tsv $WRKDIR/gtex.phenos.tsv

# MESA
export MESADIR=$BASEDIR/data_and_cohorts/mesa/mesa_phenotypes/MESA_tsvs
$CODEDIR/scripts/sample_info/phenotypes/curate_mesa_phenotypes.R \
  --first-exam-tsv $MESADIR/phs000209.v13.pht001116.v10.p3.c1.MESA_Exam1Main.HMB.txt.gz \
  --followup-tsv $MESADIR/phs000209.v13.pht001118.v8.p3.c1.MESA_Exam2Main.HMB.txt.gz \
  --followup-tsv $MESADIR/phs000209.v13.pht001119.v8.p3.c1.MESA_Exam3Main.HMB.txt.gz \
  --followup-tsv $MESADIR/phs000209.v13.pht001120.v10.p3.c1.MESA_Exam4Main.HMB.txt.gz \
  --followup-tsv $MESADIR/phs000209.v13.pht003091.v3.p3.c1.MESA_Exam5Main.HMB.txt.gz \
  --family-exam-tsv $MESADIR/../../PhenotypeFiles/phs000209.v13.phenotypes/phs000209.v13.pht001121.v3.p3.c1.MESA_FamilyExamMain.HMB.txt.gz \
  --out-tsv $WRKDIR/mesa.phenos.tsv

# LCINS
export LCINSDIR=$BASEDIR/data_and_cohorts/lcins/NatGenet_LCINS_phenotypes
$CODEDIR/scripts/sample_info/phenotypes/curate_lcins_phenotypes.R \
  --supp-table-1-tsv $LCINSDIR/LCINS.supp_table1.tsv \
  --subject-phenotypes-tsv $LCINSDIR/phs001697.v1.pht010578.v1.p1.c3.EAGLE_Never_Smokers_Subject_Phenotypes.GRU.txt.gz \
  --out-tsv $WRKDIR/lcins.phenos.tsv



#################
# OTHER|BESPOKE #
#################

# HGSVC
export HGSVCDIR=$BASEDIR/data_and_cohorts/hgsvc
$CODEDIR/scripts/sample_info/phenotypes/curate_hgsvc_phenotypes.R \
  --metadata-tsv $HGSVCDIR/HGSV_metadata.DFCI_G2C_selected.tsv \
  --out-tsv $WRKDIR/hgsvc.phenos.tsv

# HMF
export HMFDIR=$BASEDIR/data_and_cohorts/hmf/hmf-dr-355
$CODEDIR/scripts/sample_info/phenotypes/curate_hmf_phenotypes.R \
  --preprocessed-tsv $HMFDIR/HMF.combined_metadata_and_manifest.unique_patients.tsv.gz \
  --cancer-classification-tsv $HMFDIR/hmf-dr-355.g2c_classification_scheme.tsv \
  --out-tsv $WRKDIR/hmf.phenos.tsv

# ICGC
export ICGCDIR=$BASEDIR/data_and_cohorts/icgc/icgc_wgs_download_may26_2023/icgc_donor_metadata_release_28_may_2023
export TCGADIR=$BASEDIR/data_and_cohorts/icgc/tcga_sarcoma_thyroid_oct2024/
ICGC_PROCDIR=`mktemp -d`
$CODEDIR/scripts/sample_info/phenotypes/curate_icgc_phenotypes.R \
  --donors-tsv $ICGCDIR/donor.all_projects.tsv.gz \
  --repository-tsv $ICGCDIR/../repository_1685119483.tsv \
  --project-cancer-map $BASEDIR/data_and_cohorts/icgc/icgc_project_to_g2c_map.tsv \
  --specimen-tsv $ICGCDIR/specimen.all_projects.tsv.gz \
  --exposure-tsv $ICGCDIR/donor_exposure.all_projects.tsv.gz \
  --out-tsv $ICGC_PROCDIR/icgc.phenos.tsv
for cancer in thyroid sarcoma; do
  $CODEDIR/scripts/sample_info/phenotypes/curate_nci_gdc_phenotypes.R \
    --clinical-tsv $TCGADIR/$cancer/clinical.tsv \
    --exposure-tsv $TCGADIR/$cancer/exposure.tsv \
    --biospecimen-tsv $TCGADIR/$cancer/sample.tsv \
    --cohort icgc \
    --out-tsv $ICGC_PROCDIR/tcga.$cancer.phenos.tsv
done
cat \
  $ICGC_PROCDIR/icgc.phenos.tsv \
  <( sed '1d' $ICGC_PROCDIR/tcga.thyroid.phenos.tsv ) \
  <( sed '1d' $ICGC_PROCDIR/tcga.sarcoma.phenos.tsv ) \
> $WRKDIR/icgc.phenos.tsv
rm -rf $ICGC_PROCDIR

# Proactive
gsutil -m cat gs://dfci-g2c-inputs/sample-lists/proactive-*.samples.list \
| sort | uniq > $WRKDIR/proactive.samples.list
export PROCADIR="~/Partners HealthCare Dropbox/Ryan Collins/PROCA_DATA"
export PROCATERRA="$BASEDIR/data_and_cohorts/proactive/GARBER-PROACTIVE_WGS_terra_metadata_sept24_2024"
$CODEDIR/scripts/sample_info/phenotypes/curate_proactive_phenotypes.R \
  --meta-xlsx "$PROCADIR/PROCA IDs for PROACTIVE Data_15Mar24.xlsx" \
  --meta-xlsx "$PROCADIR/TOP PROACTIVE Cohort_PROCA IDs_4June24.xlsx" \
  --meta-xlsx "$PROCADIR/PROCA IDs for PROACTIVE Data_8Feb24.xlsx" \
  --meta-xlsx "$PROCADIR/PROCA IDs for PROACTIVE Data_7Nov23.xlsx" \
  --meta-xlsx "$PROCADIR/PROCA IDs for PROACTIVE Data_21Nov23.xlsx" \
  --meta-xlsx "$PROCADIR/PROACTIVE All Consented 15Aug23_for Jackie.xlsx" \
  --participant-tsv "$PROCATERRA/participant.tsv" \
  --sample-tsv "$PROCATERRA/sample.tsv" \
  --processed-samples $WRKDIR/proactive.samples.list \
  --young-lung-meta $BASEDIR/../jackie_younglung/younglung_metadata/combined_data_ALL.tsv \
  --out-dir $WRKDIR

# UFC
# UFC phenotypes were curated via manual chart review and are stored on RLCs laptop
# Thus, they are directly copied over here
cp $BASEDIR/data_and_cohorts/ufc/ufc.phenos.tsv $WRKDIR/

###################
# STAGE & CLEANUP #
###################

# Pool all phenotype data
head -n1 $( find $WRKDIR/ -name "*.phenos.tsv" | head -n1 ) \
> $WRKDIR/dfci-g2c.non_aou.phenos.tsv
find $WRKDIR/ -name "*.phenos.tsv" \
| xargs -I {} sed '1d' {} \
| sort -Vk2,2 -k1,1V \
>> $WRKDIR/dfci-g2c.non_aou.phenos.tsv

# Compress all phenotype data
find $WRKDIR/ -name "*.phenos.tsv" | xargs -I {} gzip -f {}

# Copy all per-cohort phenotype files to main (non-AoU) project bucket for storage
gsutil -m cp \
  $WRKDIR/*.phenos.tsv.gz \
  gs://dfci-g2c-inputs/phenotypes/

# Clean up working directory
rm -rf $WRKDIR

