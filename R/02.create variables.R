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
    "/ORIEN_data_with_ovarian_sample_2025-08-22.rds")) %>% 
  janitor::clean_names()


################################################################################# II ### Create new variables
gyn_data <- clinical_data %>% 
  # Create new var
  mutate(histotype = case_when(
    histology == "Serous cystadenocarcinoma, NOS"                     ~ "Serous",
    histology == "Endometrioid adenocarcinoma, NOS"                  ~ "Endometrioid",
    histology == "Papillary serous cystadenocarcinoma"                  ~ "Serous",
    histology == "Clear cell adenocarcinoma, NOS"                           ~ "Clear cell",
    histology == "Serous surface papillary carcinoma"                  ~ "Serous",
    histology == "Mixed cell adenocarcinoma"                         ~ "Other epithelial",
    histology == "Adenocarcinoma, NOS"                                  ~ "Other epithelial",
    histology == "Carcinosarcoma, NOS"                                  ~ "Other epithelial",
    histology == "Mucinous adenocarcinoma"                                  ~ "Mucinous",
    histology == "Carcinoma, NOS"                                         ~ "Other epithelial",
    histology == "Mullerian mixed tumor"                                  ~ "Other epithelial",
    histology == "Carcinoma, metastatic, NOS"                           ~ "Other epithelial",
    histology == "Mucinous cystadenocarcinoma, NOS"                  ~ "Mucinous",
    histology == "Neoplasm, malignant"                                  ~ "Exclude",
    histology == "Papillary adenocarcinoma, NOS"                           ~ "Serous",
    histology == "Squamous cell carcinoma, NOS"                           ~ "Other epithelial",
    histology == "Transitional cell carcinoma, NOS"                  ~ "Serous",
    histology == "Adenocarcinoma, intestinal type"                           ~ "Since this is intestinal type, it is probably a metastasis and could either be ovarian or colorectal primary. We could include for now as mucinous but may want to consider excluding.",
    histology == "Carcinoma, anaplastic, NOS"                           ~ "Other epithelial",
    histology == "Clear cell adenocarcinofibroma"                           ~ "Exclude",
    histology == "Mesodermal mixed tumor"                                      ~ "Other epithelial",
    histology == "Serous endometrial/tubal intraepithelial carcinoma"         ~ "Exclude",
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
        # "Grade cannot be assessed; Unknown",
        "Unknown/Not Applicable"
      )                                                ~ "High - Unknown Serous",
    TRUE                                               ~ grade_pathological
  )) %>% 
  mutate(hgsoc_vs_others = case_when(
    grade == "High"                                    ~ "High-grade serous/carcinosarcoma",
    !is.na(grade)                                      ~ "Others"
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
           ~ str_detect(., "carboplatin"))          ~ "Yes"
  )) %>%
  mutate(received_taxane_any_regimen = case_when(
    if_any(any_of(contains("regimen_")), 
           ~ str_detect(., "taxel"))                ~ "Yes"
  )) %>%
  mutate(received_bev_any_regimen = case_when(
    if_any(any_of(contains("regimen_")), 
           ~ str_detect(., "bevacizumab"))          ~ "Yes"
  )) %>% 
  mutate(received_carboplatin_regimen1 = case_when(
    str_detect(regimen_1, "carboplatin")            ~ "Yes"
  )) %>%
  mutate(received_taxane_regimen1 = case_when(
    str_detect(regimen_1, "taxel")                  ~ "Yes"
  )) %>%
  mutate(received_bev_regimen1 = case_when(
    str_detect(regimen_1, "bevacizumab")            ~ "Yes"
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

