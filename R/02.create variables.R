# Import Library
library(tidyverse)
library(lubridate)


################################################################################# I ### Load data
path_raw <- fs::path("", "Volumes", "Lab_Gillis",
                     "ORIEN_Avatar", "resources")
# load("~/Documents/GitHub/Gillis/ch_orien/Rmd/R01 Gillis-Teng 2023-2024/initial_files.RData")

clinical_data <- 
  read_rds(paste0(
    here::here(), "/data/processed data",
    "/ORIEN_data_with_ovarian_sample_20260127.rds")) %>% 
  janitor::clean_names()
old_histo <- 
  read_rds(paste0(
    here::here(), "/data/processed data/Oldversion",
    "/ORIEN_data_with_ovarian_sample_2025-09-04.rds")) %>% 
  janitor::clean_names() %>% 
  select(avatar_key, old_histo_bev = histology_behavior, old_histo = histology)

################################################################################# II ### Create new variables
gyn_data <- clinical_data %>% 
  mutate(raceeth = factor(raceeth, levels = c("White Non-Hispanic",
                                                "Black Non-Hispanic",
                                                "Hispanic any race"))) %>% 
  mutate(bmi_cat = case_when(
    bmi_cat == "Underweight" |
      bmi_cat == "Healthy"                         ~ "Underweight - Healthy",
    !is.na(bmi_cat)                                ~ bmi_cat
  ), bmi_cat = factor(bmi_cat, levels = c("Underweight - Healthy",
                                          "Overweight",
                                          "Obese"))) %>% 
  mutate(clin_group_stage = case_when(
    str_detect(clin_group_stage, "IV")                        ~ "4",
    # should not be any 5 but just in case
    str_detect(clin_group_stage, "V")                         ~ "5",
    str_detect(clin_group_stage, "III")                       ~ "3",
    str_detect(clin_group_stage, "II")                        ~ "2",
    str_detect(clin_group_stage, "I")                         ~ "1"
  )) %>% 
  mutate(path_group_stage = case_when(
    str_detect(path_group_stage, "IV")                        ~ "4",
    # should not be any 5 but just in case
    str_detect(path_group_stage, "V")                         ~ "5",
    str_detect(path_group_stage, "III")                       ~ "3",
    str_detect(path_group_stage, "II")                        ~ "2",
    str_detect(path_group_stage, "I")                         ~ "1"
  )) %>% 
  group_by(avatar_key) %>% 
  mutate(stage = max(clin_group_stage, path_group_stage, na.rm = TRUE)) %>%
  ungroup() %>% 
  mutate(stage_cat = case_when(
    stage == "1" |
      stage == "2"                                 ~ "1-2",
    stage == "3" |
      stage == "4"                                 ~ "3-4",
  ), stage_cat = factor(stage_cat, levels = c("1-2",
                                              "3-4"))) %>% 
  # Create new var
  left_join(., old_histo, by = "avatar_key") %>% 
  mutate(histology_behavior = str_match(histology_behavior, "([1-9]+ (.*))")[,3]) %>% 
  mutate(old_histo_bev = str_match(old_histo_bev, "([1-9]+ (.*))")[,3]) %>% 
  mutate(histotype = case_when(
    
    old_histo_bev == "Endometrioid adenocarcinoma, NOS" &
      old_histo == "Adenocarcinoma, NOS"                                 ~ "Endometrioid",
    old_histo_bev == "Endometrioid adenocarcinoma, NOS" &
      old_histo == "Mixed cell adenocarcinoma"                           ~ "Endometrioid",
    
    old_histo_bev == "Serous surface papillary tumor of borderline malignancy" &
      old_histo == "Serous cystadenoma, borderline malignancy"           ~ "Exclude",
    
    old_histo_bev == "Serous cystadenocarcinoma, NOS" &
      old_histo == "Neoplasm, malignant"                                 ~ "Serous",
    old_histo_bev == "Serous cystadenocarcinoma, NOS" &
      old_histo == "Neoplasm, malignant"                                 ~ "Serous",
    old_histo_bev == "Endometrioid adenocarcinoma, NOS" &
      old_histo == "Other"                                               ~ "Endometrioid",
    old_histo_bev == "Serous surface papillary carcinoma" &
      old_histo == "Serous endometrial/tubal intraepithelial carcinoma"  ~ "Serous",
    old_histo_bev == "Clear cell adenocarcinoma, NOS" &
      old_histo == "Clear cell adenocarcinofibroma"                      ~ "Clear cell",
    
    old_histo == "Serous cystadenocarcinoma, NOS"                        ~ "Serous",
    old_histo == "Endometrioid adenocarcinoma, NOS"                      ~ "Endometrioid",
    old_histo == "Papillary serous cystadenocarcinoma"                   ~ "Serous",
    old_histo == "Clear cell adenocarcinoma, NOS"                        ~ "Clear cell",
    old_histo == "Serous surface papillary carcinoma"                    ~ "Serous",
    old_histo == "Mixed cell adenocarcinoma"                             ~ "Other epithelial",
    old_histo == "Adenocarcinoma, NOS"                                   ~ "Other epithelial",
    old_histo == "Carcinosarcoma, NOS"                                   ~ "Other epithelial",
    old_histo == "Mucinous adenocarcinoma"                               ~ "Mucinous",
    old_histo == "Carcinoma, NOS"                                        ~ "Other epithelial",
    old_histo == "Mullerian mixed tumor"                                 ~ "Other epithelial",
    old_histo == "Carcinoma, metastatic, NOS"                            ~ "Other epithelial",
    old_histo == "Mucinous cystadenocarcinoma, NOS"                  ~ "Mucinous",
    old_histo == "Neoplasm, malignant"                                  ~ "Exclude",
    old_histo == "Papillary adenocarcinoma, NOS"                           ~ "Serous",
    old_histo == "Squamous cell carcinoma, NOS"                           ~ "Other epithelial",
    old_histo == "Transitional cell carcinoma, NOS"                  ~ "Serous",
    old_histo == "Adenocarcinoma, intestinal type" &
      old_histo_bev == "Pseudomyxoma peritonei"                     ~ "Mucinous",
    old_histo == "Carcinoma, anaplastic, NOS"                           ~ "Other epithelial",
    old_histo == "Clear cell adenocarcinofibroma"                           ~ "Exclude",
    old_histo == "Mesodermal mixed tumor"                                      ~ "Other epithelial",
    old_histo == "Serous endometrial/tubal intraepithelial carcinoma"         ~ "Exclude",
    old_histo == "Granulosa cell tumor, malignant"                        ~ "Exclude",
    old_histo == "Serous cystadenoma, borderline malignancy"         ~ "Serous but not invasive so probably exclude?",
    old_histo == "Sertoli-leydig cell tumor, poorly differentiated"         ~ "Exclude",
    old_histo == "Sertoli cell carcinoma"                           ~ "Exclude",
    old_histo == "Adenocarcinoma in situ, NOS"                           ~ "Other epithelial",
    old_histo == "Large cell carcinoma, NOS"                           ~ "Other epithelial",
    old_histo == "Small cell carcinoma, NOS"                           ~ "Other epithelial",
    old_histo == "Mucinous cystic tumor of borderline malignancy"         ~ "Exclude",
    old_histo == "Other"                                             ~ "Exclude",
    old_histo == "Papillary carcinoma, NOS"                           ~ "Serous",
    old_histo == "Sebaceous adenocarcinoma"                           ~ "Exclude",
    old_histo == "Thecoma, malignant"                                    ~ "Exclude",
    old_histo == "Dysgerminoma"                                          ~ "Exclude",
    old_histo == "Mixed germ cell tumor"                                 ~ "Exclude",
    old_histo == "Steroid cell tumor, malignant"                         ~ "Exclude"
  )) %>% 
  mutate(grade = case_when(
    histotype == "Serous" &
      grade_pathological %in% c(
        "Low grade", 
        "Well differentiated")                         ~ "Low",
    histotype == "Serous" &
      grade_pathological %in% c(
        "High grade", 
        "Moderately differentiated",
        "Poorly differentiated",
        "Undifferentiated")                            ~ "High",
    histotype == "Serous" &
      grade_pathological %in% c(
        "Grade cannot be assessed; Unknown",
        "Site-specific grade system category",
        "Unknown/Not Applicable"
      )                                                ~ "High",
    histotype == "Exclude"                             ~ "Exclude"#,
    # TRUE                                               ~ grade_pathological
  ), grade = factor(grade, levels = c("Low",
                                      "High"))) %>% 
  mutate(hgsoc_vs_others = case_when(
    grade == "High" #|
      # grade == "High - Unknown Serous"                 
                                                       ~ "High-grade serous/carcinosarcoma",
    histology == "Carcinosarcoma, NOS"                 ~ "High-grade serous/carcinosarcoma",
    histotype == "Exclude"                             ~ "Exclude",
    !is.na(grade) | 
      !is.na(histotype)                                ~ "Others"
  )) %>% 
  # Age at last contact
  mutate(is_agelastcontact_last_date = case_when(
    age_at_last_contact >= age_at_lab                  ~ "Yes",
    !is.na(age_at_last_contact) & 
      is.na(age_at_lab)                                ~ "Yes",
    age_at_last_contact < age_at_lab                   ~ "No"
  )) %>% 
  # mutate(age_at_last_contact = case_when(
  #   is_agelastcontact_last_date == "Yes"            ~ age_at_last_contact,
  #   is_agelastcontact_last_date == "No"             ~ age_at_lab
  # )) %>% 
  mutate(is_age_death_last_ = case_when(
    age_at_death == age_at_last_contact                ~ "Same",
    age_at_death > age_at_last_contact                 ~ "Yes",
    !is.na(age_at_last_contact) & 
      is.na(age_at_death)                              ~ "Yes",
    age_at_death < age_at_last_contact                 ~ "No"
  )) %>% 
  mutate(age_at_last_contact = case_when(
    is_age_death_last_ == "Yes"                        ~ age_at_last_contact,
    is_age_death_last_ == "Same"                       ~ age_at_last_contact,
    is_age_death_last_ == "No"                         ~ age_at_death
  )) %>%
  mutate(is_age_death_last = case_when(
    age_at_death == age_at_last_contact                ~ "Same",
    age_at_death > age_at_last_contact                 ~ "Yes",
    !is.na(age_at_last_contact) & 
      is.na(age_at_death)                              ~ "Yes",
    age_at_death < age_at_last_contact                 ~ "No"
  )) %>% 
  # Age at first treatment
  mutate(first_treatment = case_when(
    (age_at_first_radiation < age_at_first_surgery |
       (is.na(age_at_first_surgery) &
          !is.na(age_at_first_radiation))) &
      
      (age_at_first_radiation < age_at_med_start_1 |
         (is.na(age_at_med_start_1) &
            !is.na(age_at_first_radiation)))        ~ "Radiation",
    
    (age_at_first_surgery < age_at_first_radiation | 
       (is.na(age_at_first_radiation) &
          !is.na(age_at_first_surgery))) &
      
      (age_at_first_surgery < age_at_med_start_1 |
         (is.na(age_at_med_start_1) &
            !is.na(age_at_first_surgery)))          ~ "Surgery",
    
    (age_at_med_start_1 < age_at_first_radiation | 
       (is.na(age_at_first_radiation) &
          !is.na(age_at_med_start_1))) &
      
      (age_at_med_start_1 < age_at_first_surgery |
         (is.na(age_at_first_surgery) &
            !is.na(age_at_med_start_1)))             ~ "Drugs",
    
    age_at_med_start_1 == age_at_first_radiation     ~ "Drugs",
    age_at_med_start_1 == age_at_first_surgery       ~ "Surgery"
    
    # (age_at_first_sct < age_at_first_radiation | 
    #    (is.na(age_at_first_radiation) &
    #       !is.na(age_at_first_sct))) &
    #   
    #   (age_at_first_sct < age_at_first_surgery |
    #      (is.na(age_at_first_surgery) &
    #         !is.na(age_at_first_sct))) &
    #   
    #   (age_at_first_sct < age_at_med_start_1 |
    #      (is.na(age_at_med_start_1) &
    #         !is.na(age_at_first_sct)))             ~ "SCT"
    
  )) %>% 
  mutate(age_at_first_treatment = case_when(
    first_treatment == "Radiation"                  ~ age_at_first_radiation,
    first_treatment == "Surgery"                    ~ age_at_first_surgery,
    first_treatment == "Drugs"                      ~ age_at_med_start_1#,
    # first_treatment == "SCT"                        ~ age_at_first_sct
  )) %>% 
  # Drugs----
  mutate(received_carboplatin_any_regimen = case_when(
    if_any(any_of(contains("regimen_")), 
           ~ str_detect(., "carboplatin"))          ~ "Yes",
    has_medication_data == "Yes"                    ~ "No"
  )) %>%
  mutate(received_any_platin_any_regimen = case_when(
    if_any(any_of(contains("regimen_")), 
           ~ str_detect(., "platin"))               ~ "Yes",
    has_medication_data == "Yes"                    ~ "No"
  )) %>%
  mutate(received_any_taxel_any_regimen = case_when(
    if_any(any_of(contains("regimen_")), 
           ~ str_detect(., "taxel"))                ~ "Yes",
    has_medication_data == "Yes"                    ~ "No"
  )) %>%
  mutate(received_bev_any_regimen = case_when(
    if_any(any_of(contains("regimen_")), 
           ~ str_detect(., "bevacizumab"))          ~ "Yes",
    has_medication_data == "Yes"                    ~ "No"
  )) %>% 
  mutate(received_carboplatin_regimen1 = case_when(
    str_detect(regimen_drugs_1, "carboplatin")            ~ "Yes",
    has_medication_data == "Yes"                    ~ "No"
  )) %>%
  mutate(received_any_platin_regimen1 = case_when(
    str_detect(regimen_drugs_1, "platin")                 ~ "Yes",
    has_medication_data == "Yes"                    ~ "No"
  )) %>%
  mutate(received_any_taxel_regimen1 = case_when(
    str_detect(regimen_drugs_1, "taxel")                  ~ "Yes",
    has_medication_data == "Yes"                    ~ "No"
  )) %>%
  mutate(received_bev_regimen1 = case_when(
    str_detect(regimen_drugs_1, "bevacizumab")            ~ "Yes",
    has_medication_data == "Yes"                    ~ "No"
  )) %>% 
  # OS----
  mutate(os_event = case_when(
    vital_status == "Alive"                          ~ 0,
    vital_status == "Lost to follow-up"              ~ 0,
    vital_status == "Dead"                           ~ 1
  )) %>% 
  mutate(os_age = coalesce(age_at_death, age_at_last_contact)) %>% 
  mutate(os_time_from_dx_years = os_age - age_at_diagnosis) %>% 
  mutate(os_time_from_treatment_years = os_age - age_at_first_treatment) %>% 
  # RFS
  mutate(rfs_event = pfs_event) %>% 
  mutate(rfs_age = pfs_age) %>% 
  mutate(rfs_time_from_dx_years = rfs_age - age_at_diagnosis) %>% 
  mutate(rfs_time_from_treatment_years = rfs_age - age_at_first_treatment) %>% 
  # PFS - need to add death----
  # If no PFS data (event/age), I don't include their death and will be NA
  # as we don't know if there was a progression before
  mutate(pfs_age = case_when(
    pfs_event == 0 &
      os_event == 1                                 ~ os_age,
    
    pfs_event == 1                                  ~ pfs_age,
    
    pfs_event == 0                                  ~ pfs_age
  )) %>%
  mutate(pfs_event = case_when(
    pfs_event == 0 &
      os_event == 1                                 ~ 1,
    pfs_event == 1                                  ~ 1,
    pfs_event == 0                                  ~ 0
    # !is.na(pfs_event)                               ~ pfs_event,
  )) %>%
  mutate(pfs_time_from_dx_years = pfs_age - age_at_diagnosis) %>% 
  mutate(pfs_time_from_treatment_years = pfs_age - age_at_first_treatment) %>% 
  # Metastasis
  mutate(metastase_event = case_when(
    had_metastasis == "Yes"                         ~ 1,
    !is.na(has_metastasis_data)                     ~ 0
  )) %>% 
  mutate(met_age = coalesce(age_at_metastatic_site, age_at_last_contact)) %>% 
  mutate(met_time_from_dx_years = met_age - age_at_diagnosis) %>% 
  mutate(met_time_from_treatment_years = met_age - age_at_first_treatment) %>% 
  mutate(time_dx_to_first_treatment_days = (age_at_first_treatment - age_at_diagnosis) * 365)

# All patients who have outcome data but have missing pfs event have unknown data in Outcome raw file

# Save
write_csv(gyn_data,
          paste0(here::here(), 
                 "/data/processed data",
                 "/Ovary and CH data_",
                 str_remove_all(today(), "-"), ".csv"))
write_rds(gyn_data, 
          paste0(here::here(), 
                 "/data/processed data",
                 "/Ovary and CH data_",
                 str_remove_all(today(), "-"), ".rds"))

path_save <- fs::path("", "Volumes", "Gillis_Research",
                      "Christelle Colin-Leitzinger", "Ovary CH",
                      "ovary_ch")
write_csv(gyn_data, 
          paste0(path_save, 
                 "/data/processed data",
                 "/Ovary and CH data_",
                 str_remove_all(today(), "-"), ".csv"))

write_csv(gyn_data, 
          paste0(path_save, 
                 "/data/processed data",
                 "/Ovary and CH data_",
                 str_remove_all(today(), "-"), ".rds"))

path_save <- fs::path("", "Volumes", "Peres_Research",
                      "ORIEN analysis", "Ovary CH")

write_csv(gyn_data, 
          paste0(path_save, 
                 "/data/processed data",
                 "/Ovary and CH data_",
                 str_remove_all(today(), "-"), ".csv"))

# End create variables

