# WDL to run DeepVariant on a single CRAM to generate single-sample gVCF
# Contact: Ryan_Collins@dfci.harvard.edu

# Based loosely on DNANexus' original (archived) workflow. See:
# https://github.com/dnanexus-rnd/DeepVariant-GLnexus-WDL/blob/master/wdl/DeepVariant.wdl

version 1.0

workflow DeepVariant {
  File cram
  File ref_fa
  File deepvariant_model_tar
  String deepvariant_docker

  # DV gVCF optimization parameter to reduce output file size
  # See: https://github.com/google/deepvariant/blob/r0.5/docs/deepvariant-gvcf-support.md
  Int? gq_binsize = 5
}
