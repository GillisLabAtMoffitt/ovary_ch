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
    "/ORIEN_data_with_ovarian_sample_2025-09-04.rds")) %>% 
  janitor::clean_names()


################################################################################# II ### Create new variables
gyn_data <- clinical_data %>% 
  mutate(race_eth = factor(race_eth, levels = c("White Non-Hispanic",
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
  mutate(histology_behavior = str_match(histology_behavior, "([1-9]+ (.*))")[,3]) %>% 
  mutate(histotype = case_when(
    
    histology_behavior == "Endometrioid adenocarcinoma, NOS" &
      histology == "Adenocarcinoma, NOS"                                 ~ "Endometrioid",
    histology_behavior == "Endometrioid adenocarcinoma, NOS" &
      histology == "Mixed cell adenocarcinoma"                           ~ "Endometrioid",
    
    histology_behavior == "Serous surface papillary tumor of borderline malignancy" &
      histology == "Serous cystadenoma, borderline malignancy"           ~ "Exclude",
    
    histology_behavior == "Serous cystadenocarcinoma, NOS" &
      histology == "Neoplasm, malignant"                                 ~ "Serous",
    histology_behavior == "Serous cystadenocarcinoma, NOS" &
      histology == "Neoplasm, malignant"                                 ~ "Serous",
    histology_behavior == "Endometrioid adenocarcinoma, NOS" &
      histology == "Other"                                               ~ "Endometrioid",
    histology_behavior == "Serous surface papillary carcinoma" &
      histology == "Serous endometrial/tubal intraepithelial carcinoma"  ~ "Serous",
    histology_behavior == "Clear cell adenocarcinoma, NOS" &
      histology == "Clear cell adenocarcinofibroma"                      ~ "Clear cell",
    
    histology == "Serous cystadenocarcinoma, NOS"                        ~ "Serous",
    histology == "Endometrioid adenocarcinoma, NOS"                      ~ "Endometrioid",
    histology == "Papillary serous cystadenocarcinoma"                   ~ "Serous",
    histology == "Clear cell adenocarcinoma, NOS"                        ~ "Clear cell",
    histology == "Serous surface papillary carcinoma"                    ~ "Serous",
    histology == "Mixed cell adenocarcinoma"                             ~ "Other epithelial",
    histology == "Adenocarcinoma, NOS"                                   ~ "Other epithelial",
    histology == "Carcinosarcoma, NOS"                                   ~ "Other epithelial",
    histology == "Mucinous adenocarcinoma"                               ~ "Mucinous",
    histology == "Carcinoma, NOS"                                        ~ "Other epithelial",
    histology == "Mullerian mixed tumor"                                 ~ "Other epithelial",
    histology == "Carcinoma, metastatic, NOS"                            ~ "Other epithelial",
    histology == "Mucinous cystadenocarcinoma, NOS"                  ~ "Mucinous",
    histology == "Neoplasm, malignant"                                  ~ "Exclude",
    histology == "Papillary adenocarcinoma, NOS"                           ~ "Serous",
    histology == "Squamous cell carcinoma, NOS"                           ~ "Other epithelial",
    histology == "Transitional cell carcinoma, NOS"                  ~ "Serous",
    histology == "Adenocarcinoma, intestinal type" &
      histology_behavior == "Pseudomyxoma peritonei"                     ~ "Mucinous",
    histology == "Carcinoma, anaplastic, NOS"                           ~ "Other epithelial",
    histology == "Clear cell adenocarcinofibroma"                           ~ "Exclude",
    histology == "Mesodermal mixed tumor"                                      ~ "Other epithelial",
    histology == "Serous endometrial/tubal intraepithelial carcinoma"         ~ "Exclude",
    histology == "Granulosa cell tumor, malignant"                        ~ "Exclude",
    histology == "Serous cystadenoma, borderline malignancy"         ~ "Serous but not invasive so probably exclude?",
    histology == "Sertoli-leydig cell tumor, poorly differentiated"         ~ "Exclude",
    histology == "Sertoli cell carcinoma"                           ~ "Exclude",
    histology == "Adenocarcinoma in situ, NOS"                           ~ "Other epithelial",
    histology == "Large cell carcinoma, NOS"                           ~ "Other epithelial",
    histology == "Small cell carcinoma, NOS"                           ~ "Other epithelial",
    histology == "Mucinous cystic tumor of borderline malignancy"         ~ "Exclude",
    histology == "Other"                                             ~ "Exclude",
    histology == "Papillary carcinoma, NOS"                           ~ "Serous",
    histology == "Sebaceous adenocarcinoma"                           ~ "Exclude",
    histology == "Thecoma, malignant"                                    ~ "Exclude",
    histology == "Dysgerminoma"                                          ~ "Exclude",
    histology == "Mixed germ cell tumor"                                 ~ "Exclude",
    histology == "Steroid cell tumor, malignant"                         ~ "Exclude"
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
    str_detect(regimen_1, "carboplatin")            ~ "Yes",
    has_medication_data == "Yes"                    ~ "No"
  )) %>%
  mutate(received_any_platin_regimen1 = case_when(
    str_detect(regimen_1, "platin")                 ~ "Yes",
    has_medication_data == "Yes"                    ~ "No"
  )) %>%
  mutate(received_any_taxel_regimen1 = case_when(
    str_detect(regimen_1, "taxel")                  ~ "Yes",
    has_medication_data == "Yes"                    ~ "No"
  )) %>%
  mutate(received_bev_regimen1 = case_when(
    str_detect(regimen_1, "bevacizumab")            ~ "Yes",
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
                 today(), ".csv"))
write_rds(gyn_data, 
          paste0(here::here(), 
                 "/data/processed data",
                 "/Ovary and CH data_",
                 today(), ".rds"))

path_save <- fs::path("", "Volumes", "Gillis_Research",
                      "Christelle Colin-Leitzinger", "Ovary CH",
                      "ovary_ch")
write_csv(gyn_data, 
          paste0(path_save, 
                 "/data/processed data",
                 "/Ovary and CH data_",
                 today(), ".csv"))

write_csv(gyn_data, 
          paste0(path_save, 
                 "/data/processed data",
                 "/Ovary and CH data_",
                 today(), ".rds"))

path_save <- fs::path("", "Volumes", "Peres_Research",
                      "ORIEN analysis", "Ovary CH")

write_csv(gyn_data, 
          paste0(path_save, 
                 "/data/processed data",
                 "/Ovary and CH data_",
                 today(), ".csv"))

# End create variables

