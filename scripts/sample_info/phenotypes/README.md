# DFCI G2C Database

Copyright (c) 2023-Present, [Ryan L. Collins](mailto:Ryan_Collins@dfci.harvard.edu), [Noah Fields](mailto:noah_fields@dfci.harvard.edu), and the Van Allen, Gusev, and Haigis laboratories at Dana-Farber Cancer Institute.  
Distributed under terms of the [GNU GPL v2.0 License](/LICENSE) (see `LICENSE`).  


---  

This subdirectly contains code to curate participant- and sample-level metadata for all cohorts included in G2C.  

Given the highly heterogeneous nature of the cohorts included in G2C, we have mapped all available phenotype information onto a simplified data dictionary as follows:  

1. `reported_sex` : reported sex of participant, if known. Eligible values are `male`, `female`, or `other`. Note that this is not necessarily consistent with inferred sex from WGS ploidy.  

2. `reported_race_or_ethnicity` : participant race and/or ethnicity, if either are reported. Note that this will frequently differ from WGS-inferred genetic ancestry.  

3. `age` : individual age in years, if known. Ideally this reflects age at diagnosis, if reported separately, otherwise should reflect age at recruitment.  

4. `birth_year` : year of birth, if known.  

5. `vital_status` : one-hot indicator for participant vital status, if known. Eligible values are `0` if individual is known to be deceased, `1` if individual is known to be alive at last contact, or `NA` if vital status is unknown or not reported.  

6. `age_at_last_contact` : age in years at time of last contact or death. If no follow-up was conducted or reported, this value will match `age`.  

7. `days_from_dx_to_last_contact` : number of days from diagnosis to last contact or death, if reported; `NA` if unknown, not reported, or not relevant (e.g., for cancer-free controls).  

8. `height` : participant height in centimeters, if reported, presumably at time of recruitment.  

9. `weight` : participant height in kilograms, if reported, presumably at time of recruitment.  

10. `bmi` : participant BMI, if reported or if both `height` and `weight` are available, presumably at time of recruitment.  

11. `cancer` : primary cancer diagnos(es) for each patient. Metastatic/secondary sites are excluded from consideration in these fields where clinical records are sufficiently detailed to allow disambiguation. This field can have one or more values delimited with a semicolon; note that in most cases there will be a single reported primary cancer diagnosis for a patient. Diagnoses were collapsed into the following 16 domains: `prostate`, `breast`, `lung`, `colorectal`, `melanoma`, `uterus`, `kidney`, `bladder`, `oral_cavity`, `ovary`, `cns`, `pancreas`, `esophagus`, `liver`, `stomach`, `other`. Cancer-free controls have `control` in this field.  

12. `stage` : AJCC disease stage for first primary cancer, if reported. Ideally recorded at time of diagnosis. Pathology staging preferred, followed by clinical staging, followed by pathology-reported T category. Eligible values are `0`, `I`, `II`, `III`, `IV`, `unknown` if not reported, and `NA` for cancer-free controls.  

13. `metastatic` : one-hot indicator if patient is known to have developed metastatic disease. Eligible values are `0`, `1`, `unknown` if not reported, and `NA` for cancer-free controls. Note that this value may disagree with `stage` if the patient was staged at diagnosis then subsequently developed metastatic disease.  

14. `grade` : Grade of primary tumor, if reported. Ideally reported at time of diagnosis. Eligible values are `1`, `2`, `3`, `4`, `unknown` if not reported, and `NA` for cancer-free controls. Borderline grade tumors (`GB`) were treated as grade `1`.  

15. `smoking_history` : one-hot indicator for patient smoking history, if known. Eligible values include `0` for reported never smoker, `1` for positive history of smoking, and `NA` for unreported or missing data.  

16. `cancer_icd10` : semicolon-delimited ICD-10 codes for this patient's cancer diagnoses, if reported. `NA` if missing and for cancer-free controls.  

17. `original_dx` : free text field to record original specific cancer diagnosis, which is not standardized across cohorts. `NA` for missing values and cancer-free controls.  

18. `wgs_tissue` : DNA source (i.e., tissue or biospecimen) used for germline whole-genome sequencing. `unknown` if not reported.  

---  

Note that access to individual level phenotype data records are protected by agreements and consents from various entities (e.g., dbGaP, NIH All of Us, Hartwig Medical Foundation), so no blanket access can be provided. The above data descriptors and scripts contained in this subdirectory are sufficient to recreate the phenotype data used in all analyses, however.  

