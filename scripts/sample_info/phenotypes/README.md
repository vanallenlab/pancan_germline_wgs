# DFCI G2C Database
Copyright (c) 2024 Ryan Collins & the Dana-Farber Cancer Institute
Contact: [Ryan Collins](mailto:Ryan_Collins@dfci.harvard.edu)  

---  

This subdirectly contains code to curate participant- and sample-level metadata for all cohorts included in G2C.  

Given the highly heterogeneous nature of the cohorts included in G2C, we have mapped all available phenotype information onto a simplified data dictionary as follows:  

1. `reported_sex` : reported sex of participant, if known. Eligible values are `male`, `female`, or `other`. Note that this is not necessarily consistent with inferred sex from WGS ploidy.  

2. `reported_race_or_ethnicity` : participant race and/or ethnicity, if either are reported. Note that this will frequently differ from WGS-inferred genetic ancestry. Eligible values include: TBD.  

3. `age` : individual age in years, if known. Ideally this reflects age at diagnosis, if reported separately.  

4. `birth_year` : year of birth, if known.  

5. `vital_status` : one-hot indicator for participant vital status, if known. Eligible values are `0` if individual is known to be deceased, `1` if individual is known to be alive at last contact, or `NA` if vital status is unknown or not reported.  

6. `age_at_last_contact` : age in years at time of last contact or death. If no follow-up was conducted or reported, this value will match `age`.  

7. `time_from_dx_to_last_contact` : number of days from diagnosis to last contact or death, if reported; `NA` if unknown, not reported, or not relevant (e.g., for cancer-free controls).   

8. `height` : participant height in centimeters, if reported.  

9. `weight` : participant height in kilograms, if reported.  

10. `bmi` : participant BMI, if reported or if both `height` and `weight` are available.  

Columns 11 - 26 contain one-hot indicator encodings for the primary cancer diagnos(es) for each patient. Metastatic/secondary sites are excluded from consideration in these fields where clinical records are sufficiently detailed to allow disambiguation. Note that in many cases there will be a single reported primary cancer diagnosis for a patient. Diagnoses were collapsed into the following 16 domains: `prostate`, `breast`, `lung`, `colorectal`, `melanoma`, `uterus`, `kidney`, `bladder`, `oral_cavity`, `ovary`, `cns`, `pancreas`, `esophagus`, `liver`, `stomach`, `other`.  

27. `stage_at_dx` : AJCC disease stage for first primary cancer at time of diagnosis, if reported. Eligible values are `I`, `II`, `III`, `IV`, `unknown` if not reported, and `NA` for cancer-free controls.  

28. `metastatic` : one-hot indicator if patient is known to have metastatic disease. Eligible values are `0`, `1`, `unknown` if not reported, and `NA` for cancer-free controls.  

29. `grade_at_dx` : Gleason grade for first primary tumor at time of diagnosis, if reported. Eligible values are `low`, `intermediate`, `high`, `unknown` if not reported, and `NA` for cancer-free controls.  

30. `smoking_history` : patient smoking history, if known. Eligible values are `never` for no smoking history, `former` for positive smoking history but not current smoker, and `current` for current smoker.  

31. `cancer_icd10` : ICD-10 codes for this patient's cancer diagnoses, if reported. `NA` if missing and for cancer-free controls.  

32. `original_dx` : free text field to record original specific cancer diagnosis, which is not standardized across cohorts. `NA` for missing values and cancer-free controls.  

