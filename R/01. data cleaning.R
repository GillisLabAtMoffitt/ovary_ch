# Import Library
library(tidyverse)
library(lubridate)


################################################################################# I ### Load data
path_raw <- fs::path("", "Volumes", "Lab_Gillis",
                     "ORIEN_Avatar", "resources")
load("~/Documents/GitHub/Gillis/ch_orien/Rmd/R01 Gillis-Teng 2023-2024/initial_files.RData")

# The cleaning and merging steps are in the ch_orien project
# The clean data is located here :
# data <- read_rds(paste0("~/Documents/GitHub/Gillis/ch_orien",
#                         "/Rmd/R01 Gillis-Teng 2023-2024/cleaned_all_tumor_type_data_06242024.rds"))

ch_calls <- read_csv(paste0(
  here::here(),
  "/data/raw data",
  "/Ovarian_sample_list_with_CH_status_07.14.2025.csv"))

path_drug_class <- fs::path("", "Volumes", "Gillis_Research", "Lab_Data", "BreastCHEvolution")

drug_class <- 
  readxl::read_xlsx(paste0(
    path_drug_class, 
    "/ProcessedData/BreastCH_BoltonChemoDosingNatGenet2020.xlsx"),
    sheet = "Sup. Table 1", skip = 1) %>% 
  janitor::clean_names()

path <- fs::path("", "Volumes", "Gillis_Research", "Lab_Data", "CHinOvary")


################################################################################# II ### Data cleaning
# CH calle----
ch_calls <- ch_calls %>%
  rename(AvatarKey = ORIENAvatarKey) %>%
  select(-`...1`)
  # mutate(ch_status = c(rep(c("No CH", "CH"), 285), "No CH"),
  #        ch_status = factor(ch_status, levels = c("No CH", 
  #                                                 "CH")))

# sequencing----
TumorSequencing <- TumorSequencing %>%
  filter(TumorSequencingInd == "Yes") %>%
  select(AvatarKey, AgeAtTumorSequencing) %>%
  mutate(not_real_sequencing_age = case_when(
    AgeAtTumorSequencing == "Age 90 or older"   
                                          ~ "Age 90 or older"
  )) %>%
  mutate(across(c("AgeAtTumorSequencing"), ~ case_when(
    . == "Age 90 or older"                ~ 90,
    . == "Unknown/Not Applicable"         ~ NA_real_,
    TRUE                                  ~ as.numeric(.)
  )))

Tumor_wes <-
  ClinicalMolLinkage %>%
  filter(Tumor.Germline == "Tumor") %>%
  # left_join(., MSI_marker,
  #           by = c("ORIENAvatarKey", "WES" = "WES SLID")) %>% 
  select(ORIENAvatarKey, Tumor_WES = WES,
         tumor_collection_age = Age.At.Specimen.Collection,
         Disease.Type,
         tumor_SpecimenSiteOfCollection = SpecimenSiteOfCollection#,
         # MSI_high_score
  )
Germline_wes <-
  ClinicalMolLinkage %>%
  filter(Tumor.Germline == "Germline") %>%
  select(ORIENAvatarKey, Germline_WES = WES,
         germline_collection_age = Age.At.Specimen.Collection,
         germline_SpecimenSiteOfCollection = SpecimenSiteOfCollection)

# this create a germline-tumor file with all combination germline-tumor
# We selected the good combination to do the sequencing
# So the correct combination will be selected in the next merge
ClinicalMolLinkage <- full_join(Germline_wes, Tumor_wes,
                                by = "ORIENAvatarKey")

ch_calls <- ch_calls %>%
  left_join(., ClinicalMolLinkage,
            by = c("AvatarKey" = "ORIENAvatarKey", 
                   "Germline_WES", "Tumor_WES", "Disease.Type"
            )) %>%
  mutate(across(c("germline_collection_age",
                  "tumor_collection_age"), ~ case_when(
                    . == "Age 90 or older"                ~ 90,
                    . == "Unknown/Not Applicable"         ~ NA_real_,
                    TRUE                                  ~ as.numeric(.)
                  )))
rm(Tumor_wes, Germline_wes, ClinicalMolLinkage)


# cytogenetics cleaning
CytogeneticAbnormalities <- CytogeneticAbnormalities %>%
  filter(CytogenAbnormResult == "Positive")

CytogeneticAbnormalities <-
  CytogeneticAbnormalities %>%
  select(AvatarKey, CytogenAbnormName, CytogenAbnormInd) %>%
  distinct() %>%
  # group_by(AvatarKey, hormone_therapy_start_date) %>%
  pivot_wider(id_cols = c(AvatarKey),
              names_from = CytogenAbnormName,
              values_from = CytogenAbnormInd)


# demo ----
table(PatientMaster$Race)
table(PatientMaster$Ethnicity)

demographics <- PatientMaster %>%
  filter(str_detect(AvatarKey, paste0(ch_calls$AvatarKey, collapse = "|"))) %>% 
  mutate(Race = case_when(
    str_detect(Race, "American Indian")        ~ "American Indian or Alaska Native",
    Race == "White"                            ~ "White",
    Race == "Black"                            ~ "Black",
    str_detect(Race, "Asian") |
      str_detect(Race, "Cambodian") |
      Race == "Chinese" |
      Race == "Filipino" |
      Race == "Japanese" |
      Race == "Korean" |
      Race == "Pakistani" |
      Race == "Thai" |
      Race == "Vietnamese" |
      Race == "Laotian"                        ~ "Asian",
    Race == "Hawaiian" |
      Race == "Micronesian, NOS" |
      Race == "Pacific Islander, NOS" |
      Race == "Polynesian, NOS" |
      Race == "Samoan" |
      str_detect(Race, "Tongan")               ~ "Native Hawaiian or Other Pacific Islander",
    str_detect(Race, "Unknown")                ~ "Unknown",
    Race == "Other"                            ~ "Unknown",
    TRUE                                       ~ Race
  )) %>%
  mutate(Ethnicity = case_when(
    Ethnicity == "Spanish surname only"        ~ "Non-Hispanic",
    str_detect(Ethnicity, "Non-Hispanic")      ~ "Non-Hispanic",
    str_detect(Ethnicity, "Unknown")           ~ "Unknown",
    !is.na(Ethnicity)                          ~ "Hispanic",
    TRUE                                       ~ Ethnicity
  )) %>%
  mutate(race_eth = case_when(
    Race == "White" &
      Ethnicity == "Non-Hispanic"              ~ "White Non-Hispanic",
    Race == "Black" &
      Ethnicity == "Non-Hispanic"              ~ "Black Non-Hispanic",
    Race == "Others" &
      Ethnicity == "Non-Hispanic"              ~ "Others Non-Hispanic",
    Ethnicity == "Hispanic"                    ~ "Hispanic any race"
  )) #%>% 
  # mutate(Race = factor(Race, levels = c("White", "Black", "Asian",
  #                                       "American Indian or Alaska Native", 
  #                                       "Native Hawaiian or Other Pacific Islander",
  #                                       "Unknown")))

rm(PatientMaster)


# patient history - smoking
PatientHistory <- PatientHistory %>%
  filter(str_detect(AvatarKey, paste0(ch_calls$AvatarKey, collapse = "|"))) %>% 
  mutate(SmokingStatus = case_when(
    SmokingStatus == "Never"                   ~ "Never",
    SmokingStatus == "Current" |
      str_detect(SmokingStatus, "Ever")        ~ "Ever"
  ), SmokingStatus = factor(SmokingStatus, levels = c("Never", "Ever")))


# Diagnosis
Diagnosis %>%
  select(PrimaryDiagnosisSite, PrimaryDiagnosisSiteCode) %>%
  distinct() %>%
  arrange(PrimaryDiagnosisSiteCode, PrimaryDiagnosisSite)

Diagnosis1 <- Diagnosis %>%
  filter(str_detect(AvatarKey, paste0(ch_calls$AvatarKey, collapse = "|"))) %>% 
  group_by(AvatarKey) %>% 
  mutate(number_of_dx = n(), .before = AgeAtDiagnosis) %>% 
  ungroup() %>% 
  select(AvatarKey, number_of_dx,
         AgeAtFirstContact,
         AgeAtDiagnosis, YearOfDiagnosis,
         PrimaryDiagnosisSiteCode : Histology,
         ClinGroupStage, PathGroupStage,
         GradeClinical, GradePathological,
         CurrentlySeenForPrimaryOrRecurr,
         PerformStatusAtDiagnosis, 
         OtherStagingSystem, OtherStagingValue) %>%
  # mutate(ECOG = str_match(PerformStatusAtDiagnosis, "ECOG ([:digit:])")[,2]) %>%
  # mutate(Karnofsky = str_match(PerformStatusAtDiagnosis, "Karnofsky ([:digit:].*)%")[,2]) %>%
  # mutate(ClinGroupStage_1 = case_when(
  #   str_detect(ClinGroupStage, "IV")      ~ "IV",
  #   str_detect(ClinGroupStage, "III")     ~ "III",
  #   str_detect(ClinGroupStage, "II")      ~ "II",
  #   str_detect(ClinGroupStage, "^I")      ~ "I",
  #   TRUE                                  ~ ClinGroupStage
  # )) %>% 
  
  # Few missing age at dx
  # It is ok to use the age at first contact to fill missing age at dx
  mutate(not_real_age = case_when(
    AgeAtDiagnosis == "Age 90 or older"   ~ "Age 90 or older"
  )) %>%
  mutate(across(c("AgeAtDiagnosis", 
                  "AgeAtFirstContact"), ~ case_when(
                    . == "Age 90 or older"                ~ 90,
                    . == "Unknown/Not Applicable"         ~ NA_real_,
                    TRUE                                  ~ as.numeric(.)
                  ))) %>%
  mutate(dx_age_is_first_contact_age = case_when(
    is.na(AgeAtDiagnosis) &
      !is.na(AgeAtFirstContact)           ~ "Yes"
  ), .after = AgeAtDiagnosis) %>% 
  mutate(AgeAtDiagnosis = coalesce(AgeAtDiagnosis, AgeAtFirstContact)) %>% 
  arrange(AvatarKey, AgeAtDiagnosis) %>% 
  # Create var specific to ovarian dx
  group_by(AvatarKey) %>% 
  mutate(diagnosis_sequence = row_number(AvatarKey), .after = AvatarKey) %>%
  ungroup() %>% 
  mutate(had_ovarian_cancer = case_when(
    PrimaryDiagnosisSite %in% c(
      "Ovary", "Fallopian tube", 
      "Peritoneum, NOS", 
      "Specified parts of peritoneum"
    )                                     ~ "Yes",
    PrimaryDiagnosisSite ==
      "Overlapping lesion of female genital organs" &
      Histology ==
      "Serous cystadenocarcinoma, NOS"    ~ "Yes"
  )) %>% 
  mutate(age_at_ovarian_cancer = case_when(
    had_ovarian_cancer == "Yes"           ~ AgeAtDiagnosis
  )) %>% 
  group_by(AvatarKey, had_ovarian_cancer) %>% 
  mutate(age_at_first_ovarian_cancer = first(age_at_ovarian_cancer)) %>% 
  group_by(AvatarKey) %>% 
  fill(had_ovarian_cancer,
       age_at_first_ovarian_cancer,
       .direction = "updown") %>%
  ungroup() %>% 
  filter(had_ovarian_cancer == "Yes")
# write_csv(Diagnosis1 %>% filter(is.na(age_at_ovarian_cancer)) %>% arrange(AvatarKey, AgeAtDiagnosis),
#           "Patients with no ovarian cancer.csv")


Diagnosis2 <- Diagnosis1 %>% 
  # Select dx for the first ovarian cancer 
  mutate(is_first_ovarian_dx = case_when(
    age_at_first_ovarian_cancer == 
      age_at_ovarian_cancer             ~ "Yes",
    is.na(had_ovarian_cancer)           ~ "Never ovarian dx",
    !is.na(age_at_first_ovarian_cancer) ~ "No"
  )) %>% 
  # mutate(stage = case_when(
  #   is_first_ovarian_dx == "Yes"       ~ ClinGroupStage
  # )) %>% 
  # Summarize pre and post cancer info
  group_by(AvatarKey, diagnosis_sequence) %>% 
  mutate(pre_cancer_info = case_when(
    is_first_ovarian_dx == "No" &
      age_at_first_ovarian_cancer > 
      AgeAtDiagnosis                      ~ paste0(
        "age = ", AgeAtDiagnosis, ", site = ", PrimaryDiagnosisSite, ", histology = ", Histology)
  )) %>% 
  mutate(post_cancer_info = case_when(
    is_first_ovarian_dx == "No" &
      age_at_first_ovarian_cancer <= 
      AgeAtDiagnosis                      ~ paste0(
        "age = ", AgeAtDiagnosis, ", site = ", PrimaryDiagnosisSite, ", histology = ", Histology)
  )) %>% 
  group_by(AvatarKey) %>%
  mutate_at(c("pre_cancer_info", 
              "post_cancer_info"), ~ paste0(., collapse = "; ")) %>% 
  mutate_at(c("pre_cancer_info", 
              "post_cancer_info"), ~ str_remove_all(., "; NA|NA; ")) %>% 
  mutate_at(c("pre_cancer_info", 
              "post_cancer_info"), ~ na_if(., "NA")) %>% 
  ungroup() %>% 
  # Keep rows for the first ovarian dx
  filter(is_first_ovarian_dx == "Yes") %>% 
  # 1 patient has 2 rows for the same dx but different first contact date 
  distinct(AvatarKey, age_at_ovarian_cancer, .keep_all = TRUE)


# vitals os----
VitalStatus <- VitalStatus %>%
  filter(str_detect(AvatarKey, paste0(ch_calls$AvatarKey, collapse = "|"))) %>% 
  select(AvatarKey, VitalStatus, AgeAtLastContact, AgeAtDeath) %>%
  mutate(across(c("AgeAtLastContact", "AgeAtDeath"), ~ case_when(
    . == "Age 90 or older"                ~ 90,
    . == "Unknown/Not Applicable"         ~ NA_real_,
    TRUE                                  ~ as.numeric(.)
  ))) %>%
  mutate(has_os_data = "Yes")


# outcomes pfs----
# Outcomes__ <- Outcomes
# Outcomes <- Outcomes__
Outcomes <- Outcomes %>%
  filter(str_detect(AvatarKey, paste0(ch_calls$AvatarKey, collapse = "|"))) %>% 
  # mutate(is_ovary_outcome = case_when(
  #   OutcomesPrimaryDiagnosisSite %in% c(
  #     "Ovary", "Fallopian tube", 
  #     "Peritoneum, NOS", 
  #     "Specified parts of peritoneum",
  #     "Overlapping lesion of female genital organs")      ~ "Yes"
  # )) %>% 
  mutate(across(c("AgeAtProgRecur",
                  "AgeAtCurrentDiseaseStatus"), ~ case_when(
    . == "Age 90 or older"                ~ 90,
    . == "Unknown/Not Applicable"         ~ NA_real_,
    TRUE                                  ~ as.numeric(.)
  ))) %>%
  select(AvatarKey, ProgRecurInd, #is_ovary_outcome, pfs_event, 
         AgeAtProgRecur, 
         YearOfProgRecur, AgeAtCurrentDiseaseStatus, AgeAtPerformStatusMostRecent,
         OutcomesPrimaryDiagnosisSite, OutcomesPrimaryDiagnosisSiteCode, 
         everything()) %>% 
  mutate(age_at_disease_check = coalesce(AgeAtProgRecur, AgeAtCurrentDiseaseStatus), 
         .after = AgeAtCurrentDiseaseStatus) %>% 
  arrange(AvatarKey, age_at_disease_check)
  
# First select the outcome associated to the correct dx
Outcomes1 <- Outcomes %>%
  inner_join(., Diagnosis2 %>%
               select(AvatarKey, PrimaryDiagnosisSiteCode, PrimaryDiagnosisSite, number_of_dx),
             by = c("AvatarKey", "OutcomesPrimaryDiagnosisSiteCode" = "PrimaryDiagnosisSiteCode",
                    "OutcomesPrimaryDiagnosisSite" = "PrimaryDiagnosisSite"))
  # group_by(AvatarKey, pfs_event) %>% 
  # mutate(age_at_first_ovary_progression = case_when(
  #   pfs_event == 1                   ~ first(AgeAtProgRecur)
  # )) %>% 
  # arrange(AvatarKey, OutcomesPrimaryDiagnosisSite, desc(pfs_event), desc(AgeAtCurrentDiseaseStatus)) %>% 
  # mutate(age_at_last_disease_check = case_when(
  #   pfs_event == 0 &
  #     !is.na(AgeAtCurrentDiseaseStatus)      ~ first(AgeAtCurrentDiseaseStatus)
  # )) %>% 
  # group_by(AvatarKey) %>% 
  # fill(age_at_first_ovary_progression, age_at_last_disease_check, .direction = "updown") %>% 
  # ungroup() %>% 
  # arrange(AvatarKey, OutcomesPrimaryDiagnosisSite, desc(pfs_event)) %>% 
  # distinct(AvatarKey, .keep_all = TRUE) %>% 
  # mutate(pfs_age = coalesce(age_at_first_ovary_progression, age_at_last_disease_check)) %>% 
  # select(AvatarKey, contains("_"), everything())
  
patients_with_no_pfs_age <- Outcomes1 %>% filter(is.na(age_at_disease_check)) %>% 
  select(AvatarKey) %>% distinct()
patients_with_no_pfs_age <- paste0(patients_with_no_pfs_age$AvatarKey, collapse = "|")

Outcomes2 <- Outcomes %>%
  # For the patients with no diagnosis site
  filter((!str_detect(AvatarKey, paste0(Outcomes1$AvatarKey, collapse = "|")) ) |
           (str_detect(AvatarKey, patients_with_no_pfs_age) )
         ) %>% 
  filter(OutcomesPrimaryDiagnosisSite == "Unknown/Not Applicable") %>%
  left_join(., Diagnosis2 %>% select(AvatarKey, number_of_dx), by = "AvatarKey") %>% 
  # Make sure we keep unknown outcome only for patient that had only 1 dx
  # So the unknown outcome is associated to the ovary dx
  filter(number_of_dx == 1)

Outcomes_ <- Outcomes1 %>%
  # bind the data with a primary site known as first diagnosis
  # then the one without which show some No Metastasis
  bind_rows(., Outcomes2) %>%
  # PFS
  mutate(pfs_event = case_when(
    ProgRecurInd == "Progression" |
      ProgRecurInd == "Recurrence"        ~ 1,
    ProgRecurInd == "No"                  ~ 0,
    ProgRecurInd ==
      "Unknown/Not Applicable"            ~ NA_real_
  ), .after = ProgRecurInd) %>%
  # We want to take the first event == 1 but last event == 0
  arrange(AvatarKey, desc(pfs_event), age_at_disease_check) %>% 
  group_by(AvatarKey) %>% 
  mutate(pfs_age = case_when(
    pfs_event == 1                        ~ first(age_at_disease_check),
    pfs_event == 0                        ~ last(age_at_disease_check, na_rm = TRUE),
  ), .after = age_at_disease_check) %>% 
  ungroup() %>% 
  distinct(AvatarKey, .keep_all = TRUE) %>% 
  mutate(has_outcomes_data = "Yes") %>% 
  select(AvatarKey, pfs_event, pfs_age,
         OutcomesPrimaryDiagnosisSiteCode, OutcomesPrimaryDiagnosisSite,
         has_outcomes_data)

# Why so much missing values beside unknown 
# a <- Outcomes__ %>% filter(is.na(pfs_age) & ProgRecurInd != "Unknown/Not Applicable")
# Either is not associated with ovary or doesn't have age
rm(Outcomes1, Outcomes2)


# metastasis}
MetastaticDisease <- MetastaticDisease %>%
  filter(str_detect(AvatarKey, paste0(ch_calls$AvatarKey, collapse = "|"))) %>% 
  mutate(across(c("AgeAtMetastaticSite"), ~ case_when(
    . == "Age 90 or older"                ~ 90,
    . == "Unknown/Not Applicable"         ~ NA_real_,
    TRUE                                  ~ as.numeric(.)
  ))) %>%
  # Look at if patients would have yes and no in MetastaticDiseaseInd
  # There is none 
  # group_by(AvatarKey) %>% 
  # mutate(yes_and_no = dense_rank(MetastaticDiseaseInd)) %>% 
  # ungroup() %>% 
  arrange(AvatarKey, MetsDzPrimaryDiagnosisSite, AgeAtMetastaticSite) #%>%
  # distinct(AvatarKey, MetsDzPrimaryDiagnosisSite, .keep_all = TRUE)

MetastaticDisease1 <- MetastaticDisease %>%
  # For the patient with a diagnosis site
  inner_join(., Diagnosis2 %>%
               select(AvatarKey, PrimaryDiagnosisSiteCode, PrimaryDiagnosisSite),
             by = c("AvatarKey", "MetsDzPrimaryDiagnosisSiteCode" = "PrimaryDiagnosisSiteCode",
                    "MetsDzPrimaryDiagnosisSite" = "PrimaryDiagnosisSite")) %>%
  #keep the first record for each disease site
  arrange(AvatarKey, AgeAtMetastaticSite) %>%
  distinct(AvatarKey, .keep_all = TRUE)

MetastaticDisease2 <- MetastaticDisease %>%
  # For the patients with no diagnosis site
  filter((!str_detect(AvatarKey, 
                      paste0(MetastaticDisease1$AvatarKey, collapse = "|")))) %>% 
  filter(MetastaticDiseaseSite == "Unknown/Not Applicable") %>%
  left_join(., Diagnosis2 %>% select(AvatarKey, number_of_dx), by = "AvatarKey") %>% 
  # Make sure we keep unknown only for patient that had only 1 dx
  # So the unknown it is associated to the ovary dx
  filter(number_of_dx == 1) %>% 
  arrange(AvatarKey, AgeAtMetastaticSite) %>%
  distinct(AvatarKey, .keep_all = TRUE) %>% 
  select(-number_of_dx)

MetastaticDisease <- MetastaticDisease1 %>%
  # bind the data with a primary site known as first diagnosis
  # then the one without which show some No Metastasis
  bind_rows(., MetastaticDisease2) %>%
  select(AvatarKey, MetastaticDiseaseSiteCode, MetastaticDiseaseSite,
         had_metastasis = MetastaticDiseaseInd, AgeAtMetastaticSite, MetsDzPrimaryDiagnosisSite) %>%
  distinct(AvatarKey, .keep_all = TRUE) %>%
  mutate(has_metastasis_data = "Yes")

rm(MetastaticDisease1, MetastaticDisease2)


# medication----
Medications_ <- Medications %>%
  filter(str_detect(AvatarKey, paste0(ch_calls$AvatarKey, collapse = "|"))) %>% 
  filter(MedicationInd != "(Migrated) Cannot determine from available documentation") %>%
  # Create unique rows for identifying each drug line
  group_by(AvatarKey) %>% 
  mutate(drug_row_id = "drug_row_id_00", .after = AvatarKey) %>% 
  mutate(drug_row_id2 = (1000 + row_number()), .after = AvatarKey) %>% 
  ungroup() %>% 
  unite(drug_row_id, c(drug_row_id, drug_row_id2), sep = "", remove = TRUE) %>% 
  

# Medications1 <- Medications %>%
#   # For the patient with a diagnosis site
#   inner_join(., Diagnosis2 %>%
#                select(AvatarKey, PrimaryDiagnosisSiteCode, PrimaryDiagnosisSite),
#              by = c("AvatarKey", "MedPrimaryDiagnosisSiteCode" = "PrimaryDiagnosisSiteCode",
#                     "MedPrimaryDiagnosisSite" = "PrimaryDiagnosisSite"))
# 
# Medications2 <- Medications %>%
#   # For the patients with no diagnosis site
#   filter(MedPrimaryDiagnosisSite == "Unknown/Not Applicable") %>%
#   arrange(AvatarKey, AgeAtMedStart)
# 
# Medications <- Medications1 %>%
#   # bind the data with a primary site known as first diagnosis
#   # then the one without which show some No Metastasis
#   bind_rows(., Medications2) %>%
#   group_by(AvatarKey) %>%
#   mutate(n = dense_rank(MedPrimaryDiagnosisSiteCode), .after = MedPrimaryDiagnosisSiteCode) %>%
#   ungroup() %>%
#   filter(n == 1) %>%
#   select(AvatarKey, MedPrimaryDiagnosisSiteCode, MedPrimaryDiagnosisSite,
#          MedicationInd, AgeAtMedStart, Medication, AgeAtMedStop) %>%
  mutate(across(c("AgeAtMedStart",
                  "AgeAtMedStop"), ~ case_when(
    . == "Age 90 or older"                ~ 90,
    . == "Unknown/Not Applicable"         ~ NA_real_,
    TRUE                                  ~ as.numeric(.)
  ))) %>% 
  group_by(AvatarKey) %>% 
  mutate(ever_first_med_age = min(AgeAtMedStart, AgeAtMedStop)) %>% 
  ungroup()

# We reviewed the drugs name that ovarian patients for their ovarian cancer after the first data wrangling
# Here is a script to add drugs names into the Bolton categories
drug_class <- drug_class %>% 
  add_row(drug_name = "Denosumab", narrow_drug_class_cytotoxic_only = "targeted_therapy", general_drug_class = "targeted_therapy") %>% 
  add_row(drug_name = "Gemtuzumab Ozogamicin", narrow_drug_class_cytotoxic_only = "targeted_therapy", general_drug_class = "targeted_therapy") %>% 
  add_row(drug_name = "Altretamine", narrow_drug_class_cytotoxic_only = "alkylating_agent", general_drug_class = "cytotoxic_therapy") %>% 
  add_row(drug_name = "Interferon", narrow_drug_class_cytotoxic_only = "immune_therapy", general_drug_class = "immune_therapy") %>% 
  add_row(drug_name = "Trifluridine", narrow_drug_class_cytotoxic_only = "antimetabolite", general_drug_class = "cytotoxic_therapy") %>% 
  add_row(drug_name = "Raloxifene", narrow_drug_class_cytotoxic_only = "hormonal_therapy", general_drug_class = "hormonal_therapy") %>% 
  add_row(drug_name = "Lanreotide", narrow_drug_class_cytotoxic_only = "targeted_therapy", general_drug_class = "targeted_therapy")
# write_csv(drug_class,
#           paste0(here::here(), 
#                  "/data/processed data",
#                  "/CHinOvary_Updated_BoltonChemoDosing_20251014.csv"))
# write_csv(drug_class,
#           paste0(path, 
#                  "/ProcessedData",
#                  "/CHinOvary_Updated_BoltonChemoDosing_20251014.csv"))


# Here the script to update names to fit Bolton categories
Medications_ <- Medications_ %>% 
  mutate(Medication = str_remove(
    Medication, 
    " Hydrochloride|Liposomal | Sulfate| Phosphate| Citrate| Acetate| Camsylate| Disodium| Mesylate| Ditosylate| Tosylate| Tartrate"), 
    Medication = case_when(
      str_detect(Medication, "Paclitaxel")        ~ "Paclitaxel",
      Medication == "Interferon, NOS"             ~ "Interferon",
      TRUE                                        ~ Medication
    ))
  

















Medications_raw <- Medications_


# Separate patient who didn't receive drugs
Medications_never <- Medications_ %>% 
  filter(MedicationInd == "No") %>% 
  mutate(has_first_line = case_when(
    MedicationInd == "No"                           ~ "no drug received"
  ))

# Separate patient who did receive drugs
Medications_ <- Medications_ %>%
  filter(!str_detect(AvatarKey, paste0(Medications_never$AvatarKey, collapse = "|"))) %>% 
  # Some regimen is line unknown but with same age as a known line, use to fill it up
  mutate(MedLineRegimen_temp = case_when(
    MedLineRegimen %in% c("Unknown/Not Applicable", 
                          "Unknown/Not Reported")    ~ NA_character_,
    TRUE                                             ~ MedLineRegimen
  )) %>% 
  group_by(AvatarKey, AgeAtMedStart) %>% 
  fill(MedLineRegimen_temp, .direction = "updown") %>% 
  ungroup() %>% 
  mutate(MedLineRegimen = coalesce(MedLineRegimen_temp, MedLineRegimen)) %>% 
  select(- MedLineRegimen_temp) %>% 
  # Re-code Lines by taking into account the line and the age
  # this will help to rearrange the data properly instead of by age only
  mutate(regimen_line = case_when(
    MedLineRegimen %in% c("Unknown/Not Applicable", 
                          "Unknown/Not Reported")    ~ 999,
    MedLineRegimen == "Neoadjuvant Regimen"          ~ 0,
    MedLineRegimen == "First Line/Regimen"           ~ 1,
    MedLineRegimen == "Second Line/Regimen"          ~ 2,
    MedLineRegimen == "Third Line/Regimen"           ~ 3,
    MedLineRegimen == "Fourth Line/Regimen"          ~ 4,
    MedLineRegimen == "Fifth Line/Regimen"           ~ 5,
    MedLineRegimen == "Sixth Line/Regimen"           ~ 6,
    MedLineRegimen == "Seventh Line/Regimen"         ~ 7,
    MedLineRegimen == "Eighth Line/Regimen"          ~ 8,
    MedLineRegimen == "Ninth Line/Regimen"           ~ 9,
    MedLineRegimen == "Tenth Line/Regimen"           ~ 10,
    MedLineRegimen == "Eleventh Line/Regimen"        ~ 11,
    MedLineRegimen == "Twelfth Line/Regimen"         ~ 12,
    MedLineRegimen == "Maintenance"                  ~ 90,
    MedLineRegimen == "Palliative"                   ~ 91
  ), .after = MedLineRegimen) %>% 
  arrange(AvatarKey, regimen_line, AgeAtMedStart) %>% 
  group_by(AvatarKey, MedPrimaryDiagnosisSite) %>% 
  mutate(linenumber_n = dense_rank(interaction(regimen_line, AgeAtMedStart)), .after = MedLineRegimen) %>% 
  ungroup() %>% 
  # For the case when age is NA (instead of using drop in dense_rank)
  mutate(regimen_line = coalesce(linenumber_n, regimen_line)) %>% 
  select(-linenumber_n) %>% 
  arrange(AvatarKey, regimen_line, AgeAtMedStart)
  
# For the patient with a ovary diagnosis site
Medications1 <- Medications_ %>%
  inner_join(., Diagnosis2 %>%
               select(AvatarKey, PrimaryDiagnosisSiteCode, PrimaryDiagnosisSite),
             by = c("AvatarKey", "MedPrimaryDiagnosisSiteCode" = "PrimaryDiagnosisSiteCode",
                    "MedPrimaryDiagnosisSite" = "PrimaryDiagnosisSite")) %>%
  mutate(has_first_line = case_when(
    (MedLineRegimen == "Neoadjuvant Regimen" |
       MedLineRegimen == "First Line/Regimen")      ~ "has first lines"
  )) %>% 
  group_by(AvatarKey) %>% 
  fill(has_first_line, .direction = "updown") %>% 
  ungroup() %>% 
  mutate(has_first_line = case_when(
    is.na(has_first_line)                       ~ "no first lines",
    TRUE                                            ~ has_first_line
  )) 

# Medications1_no_missing <- Medications1 %>%
#   filter(has_first_line == "has first lines")
# there is 30 patients difference
Medications1_missing <- Medications1 %>% 
  filter(has_first_line == "no first lines")


# Rescue first medications dates if no first line when patient only has 1 dx and 
# a first medication info with unknown dx site
Medications2 <- Medications_ %>%
  # For the patients with no diagnosis site
  filter((str_detect(AvatarKey, 
                      paste0(Medications1_missing$AvatarKey, collapse = "|")))) %>% 
  filter(MedPrimaryDiagnosisSite == "Unknown/Not Applicable") %>%
  left_join(., Diagnosis2 %>% select(AvatarKey, number_of_dx), by = "AvatarKey") %>% 
  # Make sure we keep unknown only for patient that had only 1 dx
  # So the unknown it is associated to the ovary dx
  filter(number_of_dx == 1) %>% 
  # arrange(AvatarKey, AgeAtMedStart) %>%
  select(-number_of_dx)
Medications1 <- Medications1 %>% 
  bind_rows(., Medications2) %>% 
  distinct()

# Add Info for patients with no medication dx site but only has 1 dx
Medications3 <- Medications_ %>%
  filter((!str_detect(AvatarKey, 
                     paste0(Medications1$AvatarKey, collapse = "|")))) %>% 
  filter(MedPrimaryDiagnosisSite == "Unknown/Not Applicable") %>%
  left_join(., Diagnosis2 %>% select(AvatarKey, number_of_dx), by = "AvatarKey") %>% 
  # Make sure we keep unknown only for patient that had only 1 dx
  # So the unknown it is associated to the ovary dx
  filter(number_of_dx == 1) %>% 
  # arrange(AvatarKey, AgeAtMedStart) %>%
  select(-number_of_dx)

Medications_clean <- bind_rows(Medications_never, Medications1, Medications3) %>% 
  distinct() %>% 
  # Recode again all the regimen number after adding other unknown site
  # Some regimen is line unknown but with same age as a known line, use to fill it up
  mutate(MedLineRegimen_temp = case_when(
    MedLineRegimen %in% c("Unknown/Not Applicable", 
                          "Unknown/Not Reported")    ~ NA_character_,
    TRUE                                             ~ MedLineRegimen
  )) %>% 
  group_by(AvatarKey, AgeAtMedStart) %>% 
  fill(MedLineRegimen_temp, .direction = "updown") %>% 
  ungroup() %>% 
  mutate(MedLineRegimen = coalesce(MedLineRegimen_temp, MedLineRegimen)) %>% 
  select(- MedLineRegimen_temp) %>% 
  # Re-code Lines by taking into account the line and the age
  # this will help to rearrange the data properly instead of by age only
  mutate(regimen_line = case_when(
    MedLineRegimen %in% c("Unknown/Not Applicable", 
                          "Unknown/Not Reported")    ~ 999,
    MedLineRegimen == "Neoadjuvant Regimen"          ~ 0,
    MedLineRegimen == "First Line/Regimen"           ~ 1,
    MedLineRegimen == "Second Line/Regimen"          ~ 2,
    MedLineRegimen == "Third Line/Regimen"           ~ 3,
    MedLineRegimen == "Fourth Line/Regimen"          ~ 4,
    MedLineRegimen == "Fifth Line/Regimen"           ~ 5,
    MedLineRegimen == "Sixth Line/Regimen"           ~ 6,
    MedLineRegimen == "Seventh Line/Regimen"         ~ 7,
    MedLineRegimen == "Eighth Line/Regimen"          ~ 8,
    MedLineRegimen == "Ninth Line/Regimen"           ~ 9,
    MedLineRegimen == "Tenth Line/Regimen"           ~ 10,
    MedLineRegimen == "Eleventh Line/Regimen"        ~ 11,
    MedLineRegimen == "Twelfth Line/Regimen"         ~ 12,
    MedLineRegimen == "Maintenance"                  ~ 90,
    MedLineRegimen == "Palliative"                   ~ 91
  ), .after = MedLineRegimen) %>% 
  arrange(AvatarKey, regimen_line, AgeAtMedStart) %>% 
  group_by(AvatarKey, MedPrimaryDiagnosisSite) %>% 
  mutate(linenumber_n = dense_rank(interaction(regimen_line, AgeAtMedStart)), .after = MedLineRegimen) %>% 
  ungroup() %>% 
  # For the case when age is NA (instead of using drop in dense_rank)
  mutate(regimen_line = coalesce(linenumber_n, regimen_line)) %>% 
  select(-linenumber_n) %>% 
  arrange(AvatarKey, regimen_line, AgeAtMedStart)

# Add Bolton drug class
Medications_clean <- Medications_clean %>% 
  mutate(Medication = str_to_lower(Medication)) %>% 
  left_join(., drug_class, by = c("Medication" = "drug_name"))  %>% # will need to fix some name like gemcitabine doxorubicin...
  # mutate(received_carboplatin = case_when(
  #   str_detect(Medication, "carboplatin")        ~ "Yes"
  # )) %>% 
  # mutate(received_taxane = case_when(
  #   str_detect(Medication, "taxel")              ~ "Yes"
  # )) %>% 
  # mutate(received_bev = case_when(
  #   str_detect(Medication, "bevacizumab")        ~ "Yes"
  # ))
  mutate(has_first_line = case_when(
    (MedLineRegimen == "First Line/Regimen" |
       MedLineRegimen == "Neoadjuvant Regimen")     ~ "has first line"
  )) %>% 
  group_by(AvatarKey) %>% 
  fill(has_first_line, .direction = "updown") %>% 
  ungroup() %>% 
  mutate(has_first_line = case_when(
    is.na(has_first_line)                           ~ "no first line",
    TRUE                                            ~ has_first_line
  )) 

# Create regimen
Medications_long <- Medications_clean %>%
  group_by(AvatarKey, regimen_line, AgeAtMedStart, MedicationInd, AgeAtMedStop, 
           ever_first_med_age, has_first_line#,
           # received_carboplatin, received_taxane, received_bev
           ) %>%
  summarise_at(vars(Medication, 
                    narrow_drug_class_cytotoxic_only, 
                    general_drug_class), str_c, collapse = "; ") %>%
  group_by(AvatarKey, regimen_line, AgeAtMedStart, MedicationInd, 
           ever_first_med_age, has_first_line#,
           # received_carboplatin, received_taxane, received_bev
  ) %>%
  summarise_at(vars(Medication, 
                    narrow_drug_class_cytotoxic_only, 
                    general_drug_class, 
                    AgeAtMedStop), str_c, collapse = "; ") %>%
  ungroup() %>% 
  rename(regimen = Medication) %>% 
  arrange(AvatarKey, AgeAtMedStart)


Medications_long <- Medications_long %>% 
  full_join(., ch_calls %>% 
              select(AvatarKey, germline_collection_age), 
            by = "AvatarKey") %>% 
  mutate(regimen_specimen_sequence = case_when(
    AgeAtMedStart <= germline_collection_age        ~ "regimen before germline",
    AgeAtMedStart > germline_collection_age         ~ "regimen after germline"
  )) %>% 
  select(-c(germline_collection_age))

write_rds(Medications_long,
          paste0(here::here(),
                 "/data/processed data",
                 "/Medication long format_",
                 today(), ".rds"))
write_rds(Medications_long,
          paste0(here::here(),
                 "/data/processed data",
                 "/Medication long format_",
                 today(), ".csv"))

library(data.table)
Medications_wide <- dcast(setDT(Medications_long),
                     AvatarKey + MedicationInd + ever_first_med_age + has_first_line #+ received_carboplatin + received_taxane + received_bev 
                     ~ rowid(AvatarKey),
                     value.var = c("AgeAtMedStart", "regimen", "AgeAtMedStop")) %>% 
  select(AvatarKey, drugs_ever = MedicationInd, ever_first_med_age, has_first_line, 
         # received_carboplatin, received_taxane, received_bev,
         starts_with("AgeAtMedStart_"),
         starts_with("regimen_"), starts_with("AgeAtMedStop_"), everything()) %>%
  # select(AvatarKey, drugs_ever = MedicationInd, AgeAtMedStart_1 : AgeAtMedStart_20,
  #        Medication_1 : Medication_20, AgeAtMedStop_1 : AgeAtMedStop_20) %>%
  mutate(has_medication_data = "Yes") %>% 
  mutate(had_prior_medication = case_when(
    AgeAtMedStart_1 > ever_first_med_age       ~ "Yes"
  ))

write_rds(Medications_wide,
          paste0(here::here(),
                 "/data/processed data",
                 "/Medication wide format_",
                 today(), ".rds"))

Medications_clean <- Medications_clean %>% 
  unite(drug_row_id_2, c(AvatarKey, drug_row_id), remove = FALSE) %>% 
  mutate(drug_sequence_vs_ovarian_drugs = "2- drugs for ovarian cancer")

Medications_raw_1 <- Medications_raw %>% 
  # Filter all drugs for previous cancer and posterior cancer
  filter((str_detect(AvatarKey, 
                     paste0(Medications_clean$AvatarKey, collapse = "|")))) %>% 
  unite(drug_row_id_2, c(AvatarKey, drug_row_id), remove = FALSE) %>% 
  filter((!str_detect(drug_row_id_2, 
                      paste0(Medications_clean_temp$drug_row_id_2, collapse = "|")))) %>% 
  # Add bolton class
  mutate(Medication = str_to_lower(Medication)) %>% 
  left_join(., drug_class, by = c("Medication" = "drug_name"))  %>%
  # Create var for pre post cancer drugs
  left_join(., Medications_wide %>% select(AvatarKey, AgeAtMedStart_1)) %>% 
  mutate(drug_sequence_vs_ovarian_drugs = case_when(
    AgeAtMedStart < AgeAtMedStart_1                 ~ "1- drugs for pre cancer",
    AgeAtMedStart > AgeAtMedStart_1                 ~ "3- drugs for post cancer"
  )) %>% 
  mutate(regimen_line = case_when(
    drug_sequence_vs_ovarian_drugs == 
      "1- drugs for pre cancer"                     ~ -1,
    drug_sequence_vs_ovarian_drugs == 
      "3- drugs for post cancer"                    ~ 88,
  )) %>% 
  select(-AgeAtMedStart_1) %>% 
  # Add all drugs for ovarian cancer
  bind_rows(., Medications_clean) %>% 
  # Create filter var for drugs before/after blood
  left_join(., ch_calls %>% 
              select(AvatarKey, germline_collection_age), 
            by = "AvatarKey") %>% 
  mutate(regimen_specimen_sequence = case_when(
    AgeAtMedStart <= germline_collection_age        ~ "regimen before germline",
    AgeAtMedStart > germline_collection_age         ~ "regimen after germline"
  )) %>% 
  select(-c(germline_collection_age)) %>% 
  arrange(AvatarKey, regimen_line) %>% 
  select(AvatarKey, regimen_line, regimen_specimen_sequence, 
         Medication, AgeAtMedStart, AgeAtMedStop, ever_first_med_age,
         narrow_drug_class_cytotoxic_only, general_drug_class,
         everything(), -drug_row_id_2)
  
write_rds(Medications_raw_1,
          paste0(here::here(),
                 "/data/processed data",
                 "/All medication all cancer for ovarian patients_",
                 today(), ".rds"))

rm(Medications_, Medications_clean, 
   Medications_never,
   Medications1, Medications1_missing, Medications1_no_missing,
   Medications2, Medications3)

# sct----
# The first exploratory analysis showed that none of 
# the patients with ovarian cancer had SCT so I remove the clean up

rm(StemCellTransplant)


# radiation----
Radiation <- Radiation %>%
  filter(str_detect(AvatarKey, paste0(ch_calls$AvatarKey, collapse = "|"))) %>% 
  filter(RadiationTherapyInd != "(Migrated) Cannot determine from available documentation") %>% 
  mutate(across(c("AgeAtRadiationStart",
                  "AgeAtRadiationStop"), ~ case_when(
                    . == "Age 90 or older"                ~ 90,
                    . == "Unknown/Not Applicable"         ~ NA_real_,
                    TRUE                                  ~ as.numeric(.)
                  ))) %>% 
  mutate(RadDose = as.numeric(RadDose)) %>% 
  group_by(AvatarKey) %>% 
  mutate(ever_first_rad_age = min(AgeAtRadiationStart)) %>% 
  ungroup()

Radiation1 <- Radiation %>%
  # For the patient with a diagnosis site
  inner_join(., Diagnosis2 %>%
               select(AvatarKey, PrimaryDiagnosisSiteCode, PrimaryDiagnosisSite),
             by = c("AvatarKey", "RadPrimaryDiagnosisSiteCode" = "PrimaryDiagnosisSiteCode",
                    "RadPrimaryDiagnosisSite" = "PrimaryDiagnosisSite")) # %>%
  # We shouldn't order by age as it would exclude unknown age of first line
  # group_by(AvatarKey) %>% 
  # mutate(first_age_is_na = first(AgeAtRadiationStart), .before = AgeAtRadiationStart) %>% 
  # ungroup()

Radiation2 <- Radiation %>%
  # For the patients with no diagnosis site
  filter((!str_detect(AvatarKey, 
                      paste0(Radiation1$AvatarKey, collapse = "|")))) %>% 
  # There is no first line to rescue from Unknown dx for patient with only 1 dx so moving forward
  filter(RadPrimaryDiagnosisSite == "Unknown/Not Applicable") %>%
  left_join(., Diagnosis2 %>% select(AvatarKey, number_of_dx), by = "AvatarKey") %>% 
  # Make sure we keep unknown only for patient that had only 1 dx
  # So the unknown it is associated to the ovary dx
  filter(number_of_dx == 1) %>% 
  select(-number_of_dx)

Radiation_long <- Radiation1 %>%
  # bind the data with a known primary site as first diagnosis
  # then the one without 
  bind_rows(., Radiation2) %>% 
  full_join(., ch_calls %>% 
              select(AvatarKey, germline_collection_age), 
            by = "AvatarKey") %>% 
  mutate(radiation_specimen_sequence = case_when(
    AgeAtRadiationStart <= germline_collection_age        ~ "radiation before germline",
    AgeAtRadiationStart > germline_collection_age         ~ "radiation after germline"
  )) %>% 
  select(-c(germline_collection_age))

write_rds(Radiation_long,
          paste0(here::here(),
                 "/data/processed data",
                 "/Radiation long format_",
                 today(), ".rds"))
write_rds(Radiation_long,
          paste0(here::here(),
                 "/data/processed data",
                 "/Radiation long format_",
                 today(), ".csv"))

library(data.table)
Radiation_wide <- dcast(setDT(Radiation_long),
                          AvatarKey + RadiationTherapyInd + ever_first_rad_age ~ rowid(AvatarKey),
                          value.var = c("AgeAtRadiationStart", "RadDose", "AgeAtRadiationStop")) %>% 
  select(AvatarKey, radiation_ever = RadiationTherapyInd, ever_first_rad_age, 
         age_at_first_radiation = AgeAtRadiationStart_1, starts_with("AgeAtRadiationStart_"),
         starts_with("RadDose_"), starts_with("AgeAtMedStop_"), everything()) %>%
  mutate(has_radiation_data = "Yes") %>% 
  mutate(had_prior_radiation = case_when(
    age_at_first_radiation > ever_first_rad_age       ~ "Yes"
  ))

write_rds(Radiation_wide,
          paste0(here::here(),
                 "/data/processed data",
                 "/Radiation wide format_",
                 today(), ".rds"))


rm(Radiation1, Radiation2)


# surgery ---- 
SurgeryBiopsy_ <- SurgeryBiopsy %>%
  filter(str_detect(AvatarKey, paste0(ch_calls$AvatarKey, collapse = "|"))) %>% 
  filter(SurgeryBiopsyInd != "(Migrated) Cannot determine from available documentation") %>% 
  mutate(across(c("AgeAtSurgeryBiopsy"), ~ case_when(
                    . == "Age 90 or older"                ~ 90,
                    . == "Unknown/Not Applicable"         ~ NA_real_,
                    TRUE                                  ~ as.numeric(.)
                  ))) %>% 
  group_by(AvatarKey, SiteTherapeutic) %>% 
  mutate(ever_first_surgery_age = case_when(
    SiteTherapeutic == "Yes"                              ~ min(AgeAtSurgeryBiopsy)
  )) %>% 
  select(AvatarKey, SiteTherapeutic, everything(), -SurgeryBiopsyInd) %>% 
  # keep only 1 record if multiple non therapeutic surg
  mutate(keep_filter = row_number()) %>% 
  group_by(AvatarKey) %>% 
  mutate(has_surgery_yes = case_when(
    SiteTherapeutic == "Yes"    ~ "Yes"
  )) %>% 
  fill(has_surgery_yes, .direction = "updown") %>% 
  ungroup() %>% 
  filter(SiteTherapeutic == "Yes" |
           (SiteTherapeutic == "No" & keep_filter == 1 & is.na(has_surgery_yes))) %>% 
  select(-keep_filter, -has_surgery_yes)

SurgeryBiopsy1 <- SurgeryBiopsy_ %>%
  # For the patient with a diagnosis site
  inner_join(., Diagnosis2 %>%
               select(AvatarKey, PrimaryDiagnosisSiteCode, PrimaryDiagnosisSite),
             by = c("AvatarKey", "PrimaryDiagnosisSiteCode",
                    "PrimaryDiagnosisSite"))

SurgeryBiopsy2 <- SurgeryBiopsy_ %>%
  # For the patients with no diagnosis site
  filter((!str_detect(AvatarKey, 
                      paste0(SurgeryBiopsy1$AvatarKey, collapse = "|")))) %>% 
  left_join(., Diagnosis2 %>% select(AvatarKey, number_of_dx), by = "AvatarKey") %>% 
  # Make sure we keep unknown only for patient that had only 1 dx
  # So the unknown it is associated to the ovary dx
  filter(number_of_dx == 1) %>% 
  select(-number_of_dx)

SurgeryBiopsy_long <- SurgeryBiopsy1 %>%
  # bind the data with a known primary site as first diagnosis
  # then the one without 
  bind_rows(., SurgeryBiopsy2) %>% 
  full_join(., ch_calls %>% 
              select(AvatarKey, germline_collection_age), 
            by = "AvatarKey") %>% 
  mutate(radiation_specimen_sequence = case_when(
    AgeAtSurgeryBiopsy <= germline_collection_age        ~ "surgery before germline",
    AgeAtSurgeryBiopsy > germline_collection_age         ~ "surgery after germline"
  )) %>% 
  select(-c(germline_collection_age))

write_rds(SurgeryBiopsy_long,
          paste0(here::here(),
                 "/data/processed data",
                 "/Surgery long format_",
                 today(), ".rds"))
write_rds(SurgeryBiopsy_long,
          paste0(here::here(),
                 "/data/processed data",
                 "/Surgery long format_",
                 today(), ".csv"))

library(data.table)
SurgeryBiopsy_wide <- dcast(setDT(SurgeryBiopsy_long),
                        AvatarKey + SiteTherapeutic + ever_first_surgery_age ~ rowid(AvatarKey),
                        value.var = c("AgeAtSurgeryBiopsy")) %>% 
  select(AvatarKey, surgery_ever = SiteTherapeutic, ever_first_surgery_age, 
         age_at_first_surgery = `1`, AgeAtSurgeryBiopsy_2 = `2`,
         AgeAtSurgeryBiopsy_3 = `3`, AgeAtSurgeryBiopsy_4 = `4`,
         AgeAtSurgeryBiopsy_5 = `5`, AgeAtSurgeryBiopsy_6 = `6`,
         AgeAtSurgeryBiopsy_7 = `7`, AgeAtSurgeryBiopsy_8 = `8`,
         AgeAtSurgeryBiopsy_9 = `9`) %>%
  mutate(has_surgery_data = "Yes") %>% 
  mutate(had_prior_surgery = case_when(
    age_at_first_surgery > ever_first_surgery_age       ~ "Yes"
  ))

rm(SurgeryBiopsy1, SurgeryBiopsy2)


# physical----
PhysicalAssessment <- PhysicalAssessment %>%
  filter(str_detect(AvatarKey, paste0(ch_calls$AvatarKey, collapse = "|"))) %>% 
  mutate(across(c("AgeAtPhysicalExam"), ~ case_when(
    . == "Age 90 or older"                ~ 90,
    . == "Unknown/Not Applicable"         ~ NA_real_,
    TRUE                                  ~ as.numeric(.)
  ))) %>%
  mutate(across(c("BodyWeight", "BodyHeight", "BMI"), ~ as.numeric(.)
  )) %>%
  # Keep closest record to dx
  # Only 3 record for NA age that are the only one recorded for the 3 patients
  # Can order by age and fill missing
  arrange(AgeAtPhysicalExam) %>% 
  group_by(AvatarKey) %>% 
  # get the closest/earliest data if NA
  fill(BodyWeight, BodyHeight, BMI, .direction = "up") %>% 
  left_join(., Diagnosis2 %>%
              select(AvatarKey, AgeAtDiagnosis),
            by = c("AvatarKey")) %>%
  mutate(int = abs(AgeAtPhysicalExam - AgeAtDiagnosis)) %>%
  arrange(AvatarKey, int) %>%
  distinct(AvatarKey, .keep_all = TRUE) %>%
  select(AvatarKey, AgeAtPhysicalExam, BodyWeight, BodyHeight, BMI) %>%
  mutate(bmi_cat = case_when(
    BMI < 18.5                  ~ "Underweight",
    BMI >= 18.5 &
      BMI < 25                  ~ "Healthy",
    BMI >= 25 &
      BMI < 30                  ~ "Overweight",
    BMI >= 30                   ~ "Obese"
  )) %>%
  mutate(bmi_cat = factor(bmi_cat, levels = c("Underweight", "Healthy", "Overweight", "Obese")))


# tumor markers----
TumorMarker_ <- TumorMarker %>%
  filter(str_detect(AvatarKey, paste0(ch_calls$AvatarKey, collapse = "|"))) %>% 
  # filter(TumorMarkerInd == "Yes") %>% 
  mutate(across(c("AgeAtTumorMarkerTest"), ~ case_when(
    . == "Age 90 or older"                ~ 90,
    . == "Unknown/Not Applicable"         ~ NA_real_,
    TRUE                                  ~ as.numeric(.)
  ))) %>%
  mutate(across(where(is.character), ~ na_if(., "Unknown/Not Applicable"))) %>%
  mutate(across(c("TMarkerValueUOM"), ~ na_if(., "Unknown"))) %>% 
  mutate(across(c("TMarkerRangeIndicator"), ~ na_if(., "Not Available")))



TumorMarker_2 <- TumorMarker_ %>%
  mutate(marker_value = case_when(
    TMarkerResult == "Value"              ~ TMarkerResultValue
  )) %>%
  mutate(marker_unit = case_when(
    TMarkerResult == "Value"              ~ TMarkerValueUOM
  )) %>% 
  mutate(marker_interpretation = case_when(
    !is.na(TMarkerRangeIndicator)          ~ TMarkerRangeIndicator,
    # rescue a few interpretation for CA-125 + some NAs
    marker_value >= TMarkerLowRange &
      marker_value <= TMarkerHighRange     ~ "WNL (Within Normal Limits)",
    marker_value < TMarkerLowRange         ~ "Low",
    marker_value > TMarkerHighRange        ~ "High"
  )) %>% 
  mutate(marker_value = case_when(
    !is.na(marker_value) &
      !is.na(marker_interpretation)       ~ str_c(marker_value, marker_interpretation, sep = " : "),
    !is.na(marker_value)                  ~ marker_value,
    !is.na(marker_interpretation)         ~ marker_interpretation
  )) %>% 
  mutate(across(c("TMarkerResult"), ~ na_if(., "Value"))) %>% 
  mutate(marker_result = coalesce(TMarkerResult, marker_value)) %>% 
  rename(marker_test = TMarkerTest)
  
  
  
  # group_by(marker_test) %>%
  # fill(TMarkerLowRange, TMarkerHighRange, .direction = "updown") %>%
  # ungroup() %>%
  # # mutate(across(c("marker_value", "TMarkerLowRange",
  # #                 "TMarkerHighRange"), ~ as.numeric(.)
  # # )) %>%
  # mutate(marker_interpretation = case_when(
  #   TMarkerResult == "Value"              ~ marker_interpretation,
  #   TRUE                                  ~ TMarkerResult
  # )) %>%
TumorMarker_3 <- TumorMarker_2 %>%
  mutate(marker_test_type = case_when(
    str_detect(marker_test, "HER2")       ~ str_match(marker_test, "(HER2Neu) (.*)")[,3],
    TRUE                                  ~ NA_character_
  )) %>%
  mutate(marker_test = case_when(
    str_detect(marker_test, "HER2")       ~ str_match(marker_test, "(HER2Neu) (.*)")[,2],
    TRUE                                  ~ marker_test
  )) %>% 
  mutate(marker_result = case_when(
    is.na(marker_test_type)               ~ marker_result,
    !is.na(marker_test_type)              ~ str_c(marker_test_type, marker_result, sep = " : ")
  )) %>% 
  select(AvatarKey, TumorMarkerInd, TumorMarkerKey, AgeAtTumorMarkerTest, 
         marker_test, marker_unit, marker_result) %>% 
  # group_by(AvatarKey, marker_test) %>%
  # mutate(er_result = case_when(
  #   marker_test == "Estrogen receptor (ER)"        ~ first(marker_result)
  # )) %>%
  # mutate(pr_result = case_when(
  #   marker_test == "Progesterone receptor (PR)"    ~ first(marker_result)
  # )) %>%
  mutate(marker_result = case_when(
    # str_detect(marker_test, "HER2") & 
    #   (str_detect(marker_result, "Positive") | 
    #      str_detect(marker_result, "Negative"))    ~ marker_result,
    str_detect(marker_test, "HER2") & 
      str_detect(marker_result, "Amplified")       ~ "Positive",
    TRUE                                           ~ marker_result
  ))# %>% 
  # mutate(her_result = case_when(
  #   str_detect(marker_test, "HER2")                ~ first(clean_her2)
  # ), .after = marker_test) %>% 
  # ungroup()


TumorMarker_3 <- TumorMarker_3 %>% 
  inner_join(., ch_calls %>% select(AvatarKey, tumor_collection_age), 
             by = "AvatarKey") %>% 
  # Take age at marker >= tumor collection ############################################# ASk -------------
  mutate(int = AgeAtTumorMarkerTest - tumor_collection_age,
         int = case_when(
           int < 0                 ~ NA_real_,
           int >= 0                ~ int
         )) %>%
  arrange(AvatarKey, int) %>%
  distinct(AvatarKey, marker_test, .keep_all = TRUE)

# Add Flow results
TumorMarker_long <- TumorMarker_3 %>%
  mutate(marker_test = str_replace_all(str_to_lower(marker_test), " ", "_")) %>% 
  left_join(., TumorMarkerFlowPanel, 
             by = c("AvatarKey", "TumorMarkerKey")) %>% 
  mutate(FlowPanelName = str_replace_all(str_to_lower(FlowPanelName), " ", "_")) %>% 
  unite(marker_test, c(marker_test, FlowPanelName), sep = "_", na.rm = TRUE) %>% 
  mutate(marker_result = coalesce(marker_result, FlowPanelResult)) %>% 
  select(-c(FlowPanelResult, tumor_collection_age, int, 
            TumorMarkerKey, AgeAtTumorMarkerTest))

write_rds(TumorMarker_long,
          paste0(here::here(),
                 "/data/processed data",
                 "/Tumor marker long format_",
                 today(), ".rds"))
write_rds(TumorMarker_long,
          paste0(here::here(),
                 "/data/processed data",
                 "/Tumor marker long format_",
                 today(), ".csv"))

TumorMarker_wide <- TumorMarker_long %>%
  # mutate_at(c("AgeAtTumorMarkerTest",
  #             "marker_result",
  #             "marker_unit"), ~ as.character(.)) %>%
  # pivot_longer(cols = c(AgeAtTumorMarkerTest, marker_result, marker_unit),
  #              names_to = "dat", 
  #              values_to = "val") %>% 
  # group_by(AvatarKey, TumorMarkerInd, marker_test) %>%
  pivot_wider(id_cols = c(AvatarKey, TumorMarkerInd),
              names_from = marker_test,
              values_from = -c(AvatarKey, TumorMarkerInd, marker_test),
              names_glue = "{marker_test}_{.value}",
              names_vary = "slowest"
              )



# # MSI
# MSI_marker <- MSI_marker %>% 
#   mutate(MSI_high_score = case_when(
#     `MSIsensor2 Score` >= 20       ~ "Yes",
#     `MSIsensor2 Score` < 20        ~ "No"
#   )) %>% 
#   select(ORIENAvatarKey, `WES SLID`, MSI_high_score)


# family history
FamilyHistory <- FamilyHistory %>%
  inner_join(., Diagnosis2 %>%
               select(AvatarKey, PrimaryDiagnosisSiteCode, PrimaryDiagnosisSite),
             by = c("AvatarKey", "FamilyRelationPrimarySiteCode" = "PrimaryDiagnosisSiteCode",
                    "FamilyRelationPrimarySite" = "PrimaryDiagnosisSite")) %>% 
  distinct(AvatarKey, FamilyRelation, .keep_all = TRUE) %>%
  select(AvatarKey, CancerInFamilyInd, FamilyRelation, FamilyRelationAgeGroup) %>% 
  group_by(AvatarKey, CancerInFamilyInd) %>% 
  summarise_at(vars(FamilyRelation, 
                    FamilyRelationAgeGroup), str_c, collapse = " + ") %>%
  ungroup()


# last date available----
last_date <- bind_rows(
  Imaging %>% select(AvatarKey, age_at_lab = AgeAtImageScan),
  Labs %>% select(AvatarKey, age_at_lab = AgeAtLabResults)
) %>%
  mutate(across(c("age_at_lab"), ~ case_when(
    . == "Age 90 or older"                ~ 90,
    . == "Unknown/Not Applicable"         ~ NA_real_,
    TRUE                                  ~ as.numeric(.)
  ))) %>%
  arrange(AvatarKey, desc(age_at_lab)) %>%
  distinct(AvatarKey, .keep_all = TRUE)

rm(Imaging, Labs)





  
# Merge data----
clinical_data <- ch_calls %>%
  full_join(., demographics, by = "AvatarKey") %>%
  full_join(., PatientHistory %>%
              select(AvatarKey, AgeAtLastUpdate, SmokingStatus, AlcoholUse, PreExistHPV),
            by = "AvatarKey") %>%
  full_join(., VitalStatus, by = "AvatarKey") %>%
  inner_join(., Diagnosis2, by = "AvatarKey") %>% 
  left_join(., CytogeneticAbnormalities, by = "AvatarKey") %>%
  left_join(., PhysicalAssessment, by = "AvatarKey") %>%
  left_join(., Medications_wide, by = "AvatarKey") %>%
  # left_join(., StemCellTransplant, by = "AvatarKey") %>%
  left_join(., Radiation_wide, by = "AvatarKey") %>%
  left_join(., SurgeryBiopsy_wide, by = "AvatarKey") %>%
  left_join(., Outcomes_, by = "AvatarKey") %>%
  left_join(., MetastaticDisease, by = "AvatarKey") %>%
  left_join(., last_date, by = "AvatarKey") %>%
  left_join(., FamilyHistory, by = "AvatarKey") %>% 
  left_join(., TumorMarker_wide, by = "AvatarKey")

rm(demographics, VitalStatus, Diagnosis,
   CytogeneticAbnormalities, TumorMarker,
   PhysicalAssessment, Medications,
   StemCellTransplant, Radiation,
   SurgeryBiopsy, MetastaticDisease,
   Outcomes, last_date, FamilyHistory,
   ClinicalMolLinkage, ch_calls)

# Save
write_csv(clinical_data,
          paste0(here::here(), 
                 "/data/processed data",
                 "/ORIEN_data_with_ovarian_sample_",
                 today(), ".csv"))
write_rds(clinical_data, 
          paste0(here::here(), 
                 "/data/processed data",
                 "/ORIEN_data_with_ovarian_sample_",
                 today(), ".rds"))


# End data cleaning

