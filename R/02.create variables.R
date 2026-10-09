# Import Library
library(tidyverse)
library(lubridate)


################################################################################# I ### Load data
path_orien <- fs::path("", "Volumes", "Lab_Gillis",
                     "ORIEN_Avatar", "resources")

clinical_data <- 
  read_rds(paste0(
    here::here(), "/data/processed_data",
    "/CHinOvary_CleanedClinicalData_20261009.rds")) %>% 
  janitor::clean_names()


################################################################################# II ### Create new variables
ovarian_data <- clinical_data %>% 
  # mutate(raceeth = factor(raceeth, levels = c("White Non-Hispanic",
  #                                               "Black Non-Hispanic",
  #                                               "Hispanic, any race"))) %>% 
  mutate(bmi_cat = case_when(
    bmi_cat == "Underweight" |
      bmi_cat == "Healthy"                         ~ "Underweight - Healthy",
    !is.na(bmi_cat)                                ~ bmi_cat
  ), bmi_cat = factor(bmi_cat, levels = c("Underweight - Healthy",
                                          "Overweight",
                                          "Obese"))) %>% 
  mutate(stage_cat = case_when(
    stage == "1" |
      stage == "2"                                 ~ "1-2",
    stage == "3" |
      stage == "4"                                 ~ "3-4",
  ), stage_cat = factor(stage_cat, levels = c("1-2",
                                              "3-4"))) %>% 
  # Create histo variables
  mutate(histology_behavior = str_match(histology_behavior, "([1-9]+ (.*))")[,3]) %>% 
  mutate(histotype = case_when(
    
    primary_diagnosis_site == "Ovary" &
      histology == "Adenocarcinoma, NOS" & 
      histology_behavior == "Serous cystadenocarcinoma, NOS"             ~ "Serous",
    primary_diagnosis_site == "Ovary" &
      histology == "Adenocarcinoma, NOS" & 
      histology_behavior == "Adenocarcinoma, NOS"                        ~ "Other epithelial",
    (primary_diagnosis_site == "Ovary" |
       primary_diagnosis_site == "Peritoneum, NOS" |
       primary_diagnosis_site == "Specified parts of peritoneum" |
       primary_diagnosis_site == "Fallopian tube") &
      histology == "Serous carcinoma, NOS"                               ~ "Serous",
    (primary_diagnosis_site == "Ovary" |
       primary_diagnosis_site == "Peritoneum, NOS") &
       histology == "Endometrioid adenocarcinoma, NOS"                   ~ "Endometrioid",
    
    (primary_diagnosis_site == "Ovary" |
       primary_diagnosis_site == "Peritoneum, NOS" |
       primary_diagnosis_site == "Fallopian tube") &
      histology == "Low-grade serous carcinoma"                          ~ "Serous", # LGSOC
    (primary_diagnosis_site == "Ovary" |
       primary_diagnosis_site == "Peritoneum, NOS" |
       primary_diagnosis_site == "Fallopian tube") &
      histology == "High-grade serous carcinoma"                         ~ "Serous", # HGSOC
    
    primary_diagnosis_site == "Ovary" &
      histology == "Clear cell adenocarcinoma, NOS"                      ~ "Clear cell",
    primary_diagnosis_site == "Ovary" &
      histology == "Clear cell adenocarcinofibroma"                      ~ "Clear cell",
    primary_diagnosis_site == "Ovary" &
      histology == "Mixed cell adenocarcinoma"                           ~ "Other epithelial",
    (primary_diagnosis_site == "Ovary" |
      primary_diagnosis_site == "Fallopian tube") &
      histology == "Mullerian mixed tumor"                               ~ "Other epithelial",
    primary_diagnosis_site == "Ovary" &
      histology == "Mesodermal mixed tumor"                              ~ "Other epithelial",
    (primary_diagnosis_site == "Ovary" |
       primary_diagnosis_site == "Fallopian tube") &
      histology == "Carcinosarcoma, NOS"                                 ~ "Other epithelial",
    primary_diagnosis_site == "Ovary" &
      histology == "Carcinoma, NOS"                                      ~ "Other epithelial",
    primary_diagnosis_site == "Ovary" &
      histology == "Mucinous adenocarcinoma"                             ~ "Mucinous",
    primary_diagnosis_site == "Ovary" &
      histology == "Mucinous cystadenocarcinoma, NOS"                    ~ "Mucinous",
    primary_diagnosis_site == "Ovary" &
      histology == "Squamous cell carcinoma, NOS"                        ~ "Other epithelial",
    primary_diagnosis_site == "Ovary" &
      histology == "Transitional cell carcinoma, NOS"                    ~ "Serous",
    primary_diagnosis_site == "Peritoneum, NOS" &
      histology == "Adenocarcinoma, NOS"                                 ~ "Other epithelial",
    (primary_diagnosis_site == "Ovary" |
       primary_diagnosis_site == "Fallopian tube") &
      histology == "Papillary adenocarcinoma, NOS"                       ~ "Other epithelial",
    primary_diagnosis_site == "Ovary" &
      histology == "Carcinoma, anaplastic, NOS"                          ~ "Other epithelial",
    primary_diagnosis_site == "Ovary" &
      histology == "Large cell carcinoma, NOS"                           ~ "Other epithelial",
    primary_diagnosis_site == "Ovary" &
      histology == "Small cell carcinoma, NOS"                           ~ "Other epithelial",
        
    #     
    # histology == "Serous cystadenocarcinoma, NOS"                        ~ "Serous",
    # histology == "Papillary serous cystadenocarcinoma"                   ~ "Serous",
    # histology == "Serous surface papillary carcinoma"                    ~ "Serous",
    # histology == "Carcinoma, metastatic, NOS"                            ~ "Other epithelial",
    # histology == "Neoplasm, malignant"                                  ~ "Exclude",
    # histology == "Adenocarcinoma, intestinal type" &
    #   histology_behavior == "Pseudomyxoma peritonei"                     ~ "Mucinous",
    # histology == "Clear cell adenocarcinofibroma"                           ~ "Exclude",
    # histology == "Serous endometrial/tubal intraepithelial carcinoma"         ~ "Exclude",
    # histology == "Granulosa cell tumor, malignant"                        ~ "Exclude",
    # histology == "Serous cystadenoma, borderline malignancy"         ~ "Serous but not invasive so probably exclude?",
    # histology == "Sertoli-leydig cell tumor, poorly differentiated"         ~ "Exclude",
    # histology == "Sertoli cell carcinoma"                           ~ "Exclude",
    # histology == "Adenocarcinoma in situ, NOS"                           ~ "Other epithelial",
    # 
    # 
    # 
    # 
    # 
    # 
    # 
    # 
    # 
    #   
    #   histology == "Adenocarcinoma, NOS"                                 ~ "Endometrioid",
    #   histology == "Mixed cell adenocarcinoma"                           ~ "Endometrioid",
    # 
    #   histology == "Serous cystadenoma, borderline malignancy"           ~ "Exclude",
    # 
    #   histology == "Neoplasm, malignant"                                 ~ "Serous",
    #   histology == "Neoplasm, malignant"                                 ~ "Serous",
    #   histology == "Other"                                               ~ "Endometrioid",
    #   histology == "Serous endometrial/tubal intraepithelial carcinoma"  ~ "Serous",
    # 
    # 
    # histology == "Mucinous cystic tumor of borderline malignancy"         ~ "Exclude",
    # histology == "Other"                                             ~ "Exclude",
    # histology == "Papillary carcinoma, NOS"                           ~ "Serous",
    # histology == "Sebaceous adenocarcinoma"                           ~ "Exclude",
    # histology == "Thecoma, malignant"                                    ~ "Exclude",
    # histology == "Dysgerminoma"                                          ~ "Exclude",
    # histology == "Mixed germ cell tumor"                                 ~ "Exclude",
    # histology == "Steroid cell tumor, malignant"                         ~ "Exclude",
    TRUE          ~ "?"
  )) %>% 
  mutate(grade_pathological = str_to_sentence(grade_pathological)) %>% 
  mutate(grade = case_when(
    (primary_diagnosis_site == "Ovary" |
       primary_diagnosis_site == "Peritoneum, NOS" |
       primary_diagnosis_site == "Fallopian tube") &
      histology == "Low-grade serous carcinoma"                          ~ "Low", # LGSOC
    (primary_diagnosis_site == "Ovary" |
       primary_diagnosis_site == "Peritoneum, NOS" |
       primary_diagnosis_site == "Fallopian tube") &
      histology == "High-grade serous carcinoma"                         ~ "High", # HGSOC
    histotype == "Serous" &
      grade_pathological %in% c(
        "Low grade", 
        "Well Differentiated")                                           ~ "Low",
    histotype == "Serous" &
      grade_pathological %in% c(
        "High grade", 
        "Moderately Differentiated",
        "Poorly Differentiated",
        "Undifferentiated")                                              ~ "High",
    histotype == "Serous" &
      grade_pathological %in% c(
        "Grade cannot be assessed; Unknown",
        "Site-specific grade system category",
        "Neoadjuvant/First Line Therapy Received",
        "Unknown/Not Applicable"
      )                                                                  ~ "High",
    histotype == "Exclude"                                               ~ "Exclude"#,
    # TRUE                                                                 ~ paste0(histotype, grade_pathological)
  ), grade = factor(grade, levels = c("Low",
                                      "High"))
  ) %>% 
  mutate(hgsoc_vs_others = case_when(
    grade == "High"                                    ~ "High-grade serous/carcinosarcoma",
    !is.na(grade) | 
      !is.na(histotype)                                ~ "Others"
  )) %>% 
  # select(avatar_key, 
  #        histology, histology_behavior, new_grade_pathological = grade_pathological,
  #        histology, histology_behavior, old_grade_pathological,
  #        histotype_christelle = histotype,
  #        grade_christelle = grade, hgsoc_vs_others_christelle = hgsoc_vs_others)
  # Age at last contact
  mutate(is_age_lastcontact_lastdate = case_when(
    age_at_last_contact >= age_at_lab                  ~ "Yes",
    !is.na(age_at_last_contact) & 
      is.na(age_at_lab)                                ~ "Yes",
    age_at_last_contact < age_at_lab                   ~ "No"
  )) %>% 
  mutate(is_deathage_last_ = case_when(
    age_at_death == age_at_last_contact                ~ "Yes",
    age_at_death > age_at_last_contact                 ~ "Yes",
    !is.na(age_at_last_contact) & 
      is.na(age_at_death)                              ~ "Yes",
    age_at_death < age_at_last_contact                 ~ "No"
  )) |> 
  # mutate(age_at_last_contact = case_when( # Not needed
  #   is_deathage_last_ == "Yes"                        ~ age_at_last_contact,
  #   is_deathage_last_ == "Same"                       ~ age_at_last_contact,
  #   is_deathage_last_ == "No"                         ~ age_at_death
  # )) %>%
  # Age at first treatment # Already done
  # mutate(first_treatment = case_when(
  #   (age_at_radiation_start_1 < age_at_surgery |
  #      (is.na(age_at_surgery) &
  #         !is.na(age_at_radiation_start_1))) &
  #     
  #     (age_at_radiation_start_1 < age_at_med_start_1 |
  #        (is.na(age_at_med_start_1) &
  #           !is.na(age_at_radiation_start_1)))        ~ "Radiation",
  #   
  #   (age_at_surgery < age_at_radiation_start_1 | 
  #      (is.na(age_at_radiation_start_1) &
  #         !is.na(age_at_surgery))) &
  #     
  #     (age_at_surgery < age_at_med_start_1 |
  #        (is.na(age_at_med_start_1) &
  #           !is.na(age_at_surgery)))          ~ "Surgery",
  #   
  #   (age_at_med_start_1 < age_at_radiation_start_1 | 
  #      (is.na(age_at_radiation_start_1) &
  #         !is.na(age_at_med_start_1))) &
  #     
  #     (age_at_med_start_1 < age_at_surgery |
  #        (is.na(age_at_surgery) &
  #           !is.na(age_at_med_start_1)))             ~ "Drugs",
  #   
  #   age_at_med_start_1 == age_at_radiation_start_1     ~ "Drugs",
  #   age_at_med_start_1 == age_at_surgery       ~ "Surgery"
  #   
  #   # (age_at_first_sct < age_at_radiation_start_1 | 
  #   #    (is.na(age_at_radiation_start_1) &
  #   #       !is.na(age_at_first_sct))) &
  #   #   
  #   #   (age_at_first_sct < age_at_surgery |
  #   #      (is.na(age_at_surgery) &
  #   #         !is.na(age_at_first_sct))) &
  #   #   
  #   #   (age_at_first_sct < age_at_med_start_1 |
  #   #      (is.na(age_at_med_start_1) &
  #   #         !is.na(age_at_first_sct)))             ~ "SCT"
  #   
  # )) %>% 
  # mutate(age_at_first_treatment = case_when(
  #   first_treatment == "Radiation"                  ~ age_at_radiation_start_1,
  #   first_treatment == "Surgery"                    ~ age_at_surgery,
  #   first_treatment == "Drugs"                      ~ age_at_med_start_1#,
  #   # first_treatment == "SCT"                        ~ age_at_first_sct
  # )) %>% 
  mutate(had_prior_radiation = case_when(
    age_at_radiation_start_1 > ever_first_rad_age       ~ "Yes"
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
    str_detect(regimen_drugname_1, "carboplatin")            ~ "Yes",
    has_medication_data == "Yes"                    ~ "No"
  )) %>%
  mutate(received_any_platin_regimen1 = case_when(
    str_detect(regimen_drugname_1, "platin")                 ~ "Yes",
    has_medication_data == "Yes"                    ~ "No"
  )) %>%
  mutate(received_any_taxel_regimen1 = case_when(
    str_detect(regimen_drugname_1, "taxel")                  ~ "Yes",
    has_medication_data == "Yes"                    ~ "No"
  )) %>%
  mutate(received_bev_regimen1 = case_when(
    str_detect(regimen_drugname_1, "bevacizumab")            ~ "Yes",
    has_medication_data == "Yes"                    ~ "No"
  )) |> 
  # OS----
  mutate(os_event = case_when(
    vital_status == "Alive"                          ~ 0,
    vital_status == "Lost to follow-up"              ~ 0,
    vital_status == "Dead"                           ~ 1,
    has_os_data == "Yes"                             ~ 98,
    !is.na(vital_status)                             ~ 99
  )) %>% 
  mutate(os_age = coalesce(age_at_death, age_at_last_contact)) %>% 
  mutate(os_time_from_dx_years = os_age - age_at_diagnosis) %>% 
  mutate(os_time_from_treatment_years = os_age - age_at_first_treatment) %>% 
  # RFS
  mutate(rfs_event = case_when(
    !is.na(prog_recur_ind)                          ~ 1,
    has_outcomes_data == "Yes" &
      is.na(prog_recur_ind)                         ~ 0
  )) %>% 
  mutate(rfs_age = age_at_prog_recur) %>% 
  mutate(rfs_time_from_dx_years = rfs_age - age_at_diagnosis) %>% 
  mutate(rfs_time_from_treatment_years = rfs_age - age_at_first_treatment) |> 
  # PFS - need to add death----
  # If no PFS data (event/age), I don't include their death and will be NA
  # as we don't know if there was a progression before
  mutate(pfs_event = case_when(
    rfs_event == 1                                  ~ 1,
    rfs_event == 0 &
      os_event == 1                                 ~ 1,
    rfs_event == 0                                  ~ 0
    # !is.na(pfs_event)                               ~ pfs_event,
  )) %>%
  mutate(pfs_age = case_when(
    rfs_event == 1                                  ~ rfs_age,
    rfs_event == 0 &
      os_event == 1                                 ~ os_age,
    rfs_event == 0                                  ~ rfs_age
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
  mutate(time_dx_to_first_treatment_days = (age_at_first_treatment - age_at_diagnosis) * 365) %>% 
  mutate(time_blood_to_first_treatment_days = (age_at_first_treatment - germline_collection_age) * 365)
  
# All patients who have outcome data but have missing pfs event have unknown data in Outcome raw file

# Save
write_csv(ovarian_data,
          paste0(here::here(), 
                 "/data/processed_data",
                 "/CHinOvary_FinalClinicalData_",
                 str_remove_all(today(), "-"), ".csv"))
write_rds(ovarian_data, 
          paste0(here::here(), 
                 "/data/processed_data",
                 "/CHinOvary_FinalClinicalData_",
                 str_remove_all(today(), "-"), ".rds"))

path_gillis <- fs::path("", "Volumes", "Gillis_Research",
                        "Lab_Data", "CHinOvary")
write_csv(ovarian_data, 
          paste0(path_gillis, 
                 "/ProcessedData",
                 "/CHinOvary_FinalClinicalData_",
                 str_remove_all(today(), "-"), ".csv"))

write_csv(ovarian_data, 
          paste0(path_gillis, 
                 "/ProcessedData",
                 "/CHinOvary_FinalClinicalData_",
                 str_remove_all(today(), "-"), ".rds"))

path_peres <- fs::path("", "Volumes", "Peres_Research",
                       "ORIEN analysis", "Ovary CH")

write_csv(ovarian_data, 
          paste0(path_peres, 
                 "/data/processed data",
                 "/CHinOvary_FinalClinicalData_",
                 str_remove_all(today(), "-"), ".csv"))

# End create variables

