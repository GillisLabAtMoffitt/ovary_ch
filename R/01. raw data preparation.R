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
  "/Ovarian_sample_list_with_CH_status_06.24.2025.csv"))


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
            by = c("AvatarKey" = "ORIENAvatarKey", "Germline_WES", "Tumor_WES"
            )) %>%
  mutate(across(c("germline_collection_age",
                  "tumor_collection_age"), ~ case_when(
                    . == "Age 90 or older"                ~ 90,
                    . == "Unknown/Not Applicable"         ~ NA_real_,
                    TRUE                                  ~ as.numeric(.)
                  )))
rm(Tumor_wes, Germline_wes)


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

Diagnosis <- Diagnosis %>%
  select(AvatarKey, AgeAtDiagnosis, YearOfDiagnosis,
         PrimaryDiagnosisSiteCode : Histology,
         ClinGroupStage,
         CurrentlySeenForPrimaryOrRecurr,
         PerformStatusAtDiagnosis, 
         OtherStagingSystem, OtherStagingValue) %>%
  mutate(not_real_age = case_when(
    AgeAtDiagnosis == "Age 90 or older"   ~ "Age 90 or older"
  )) %>%
  mutate(across(c("AgeAtDiagnosis"), ~ case_when(
    . == "Age 90 or older"                ~ 90,
    . == "Unknown/Not Applicable"         ~ NA_real_,
    TRUE                                  ~ as.numeric(.)
  ))) %>%
  mutate(ECOG = str_match(PerformStatusAtDiagnosis, "ECOG ([:digit:])")[,2]) %>%
  mutate(Karnofsky = str_match(PerformStatusAtDiagnosis, "Karnofsky ([:digit:].*)%")[,2]) %>%
  arrange(AvatarKey, AgeAtDiagnosis) %>%
  group_by(AvatarKey) %>%
  mutate(diagnosis_sequence = row_number(AvatarKey), .after = AvatarKey) %>%
  mutate(has_subsequent_cancer = case_when(
    diagnosis_sequence > 1                ~ "Yes"
  )) %>%
  fill(has_subsequent_cancer, .direction = "updown") %>%
  ungroup()

subsequent_diagnosis <- Diagnosis %>%
  group_by(AvatarKey) %>%
  summarise_at(vars(PrimaryDiagnosisSite, PrimaryDiagnosisSiteCode,
                    AgeAtDiagnosis), str_c, collapse = "; ") %>%
  ungroup() %>%
  separate_wider_delim(cols = PrimaryDiagnosisSite, delim = "; ",
                       names = c("PrimaryDiagnosisSite", "subsequent_cancer_PrimaryDiagnosisSite"),
                       too_few = "align_start", too_many = "merge",
                       cols_remove = TRUE) %>%
  separate_wider_delim(cols = PrimaryDiagnosisSiteCode, delim = "; ",
                       names = c("PrimaryDiagnosisSiteCode", "subsequent_cancer_PrimaryDiagnosisSiteCode"),
                       too_few = "align_start", too_many = "merge",
                       cols_remove = TRUE) %>%
  separate_wider_delim(cols = AgeAtDiagnosis, delim = "; ",
                       names = c("AgeAtDiagnosis", "subsequent_cancer_AgeAtDiagnosis"),
                       too_few = "align_start", too_many = "merge",
                       cols_remove = TRUE) %>%
  select(-c(PrimaryDiagnosisSite, PrimaryDiagnosisSiteCode,
            AgeAtDiagnosis))

Diagnosis <- Diagnosis %>%
  distinct(AvatarKey, .keep_all = TRUE) %>%
  full_join(., subsequent_diagnosis,
            by = c("AvatarKey"))

rm(subsequent_diagnosis)


# vitals
VitalStatus <- VitalStatus %>%
  select(AvatarKey, VitalStatus, AgeAtLastContact, AgeAtDeath) %>%
  mutate(across(c("AgeAtLastContact", "AgeAtDeath"), ~ case_when(
    . == "Age 90 or older"                ~ 90,
    . == "Unknown/Not Applicable"         ~ NA_real_,
    TRUE                                  ~ as.numeric(.)
  ))) %>%
  mutate(has_os_data = "Yes")


# outcomes
Outcomes1 <- Outcomes %>%
  inner_join(., Diagnosis %>%
               select(AvatarKey, PrimaryDiagnosisSiteCode, PrimaryDiagnosisSite),
             by = c("AvatarKey", "OutcomesPrimaryDiagnosisSiteCode" = "PrimaryDiagnosisSiteCode",
                    "OutcomesPrimaryDiagnosisSite" = "PrimaryDiagnosisSite")) %>%
  arrange(AvatarKey, AgeAtProgRecur) %>%
  distinct(AvatarKey, .keep_all = TRUE)

Outcomes2 <- Outcomes %>%
  # For the patients with no diagnosis site
  filter(OutcomesPrimaryDiagnosisSite == "Unknown/Not Applicable") %>%
  arrange(AvatarKey, AgeAtProgRecur) %>%
  distinct(AvatarKey, .keep_all = TRUE)

Outcomes <- Outcomes1 %>%
  # bind the data with a primary site known as first diagnosis
  # then the one without which show some No Metastasis
  bind_rows(., Outcomes2) %>%
  select(AvatarKey, OutcomesPrimaryDiagnosisSiteCode, OutcomesPrimaryDiagnosisSite,
         ProgRecurInd, AgeAtProgRecur, RelapseStatus) %>%
  mutate(across(c("AgeAtProgRecur"), ~ case_when(
    . == "Age 90 or older"                ~ 90,
    . == "Unknown/Not Applicable"         ~ NA_real_,
    TRUE                                  ~ as.numeric(.)
  ))) %>%
  distinct(AvatarKey, .keep_all = TRUE) %>%
  mutate(has_outcomes_data = "Yes")

rm(Outcomes1, Outcomes2)


# metastasis}
MetastaticDisease1 <- MetastaticDisease %>%
  # For the patient with a diagnosis site
  inner_join(., Diagnosis %>%
               select(AvatarKey, PrimaryDiagnosisSiteCode, PrimaryDiagnosisSite),
             by = c("AvatarKey", "MetastaticDiseaseSiteCode" = "PrimaryDiagnosisSiteCode",
                    "MetastaticDiseaseSite" = "PrimaryDiagnosisSite")) %>%
  arrange(AvatarKey, AgeAtMetastaticSite) %>%
  distinct(AvatarKey, .keep_all = TRUE)

MetastaticDisease2 <- MetastaticDisease %>%
  # For the patients with no diagnosis site
  filter(MetastaticDiseaseSite == "Unknown/Not Applicable") %>%
  arrange(AvatarKey, AgeAtMetastaticSite) %>%
  distinct(AvatarKey, .keep_all = TRUE)

MetastaticDisease <- MetastaticDisease1 %>%
  # bind the data with a primary site known as first diagnosis
  # then the one without which show some No Metastasis
  bind_rows(., MetastaticDisease2) %>%
  select(AvatarKey, MetastaticDiseaseSiteCode, MetastaticDiseaseSite,
         had_metastasis = MetastaticDiseaseInd, AgeAtMetastaticSite, MetsDzPrimaryDiagnosisSite) %>%
  mutate(across(c("AgeAtMetastaticSite"), ~ case_when(
    . == "Age 90 or older"                ~ 90,
    . == "Unknown/Not Applicable"         ~ NA_real_,
    TRUE                                  ~ as.numeric(.)
  ))) %>%
  distinct(AvatarKey, .keep_all = TRUE) %>%
  mutate(has_metastasis_data = "Yes")

rm(MetastaticDisease1, MetastaticDisease2)


# medication----
Medications <- Medications %>%
  filter(MedicationInd != "(Migrated) Cannot determine from available documentation") %>%
    arrange(AvatarKey, AgeAtMedStart) %>% 

# Medications1 <- Medications %>%
#   # For the patient with a diagnosis site
#   inner_join(., Diagnosis %>%
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
  mutate(across(c("AgeAtMedStart"), ~ case_when(
    . == "Age 90 or older"                ~ 90,
    . == "Unknown/Not Applicable"         ~ NA_real_,
    TRUE                                  ~ as.numeric(.)
  )))

Medications1 <- Medications %>%
  # Create regimen
  arrange(AvatarKey, AgeAtMedStart, Medication, AgeAtMedStop) %>%
  group_by(AvatarKey, AgeAtMedStart, MedicationInd, AgeAtMedStop) %>%
  summarise_at(vars(Medication), str_c, collapse = "; ") %>%
  group_by(AvatarKey, AgeAtMedStart, MedicationInd) %>%
  summarise_at(vars(Medication, AgeAtMedStop), str_c, collapse = "; ") %>%
  ungroup() %>% 
  rename(regimen = Medication)

# write_rds(Medications1,
#           paste0(here::here(), 
#                  "/data/processed data",
#                  "/Medication long format_",
#                  today(), ".rds"))

Medications1 <- Medications1 %>% 
  full_join(., ch_calls %>% 
              select(AvatarKey, Disease.Type, germline_collection_age), 
            by = "AvatarKey") %>% 
  mutate(regimen_specimen_sequence = case_when(
    AgeAtMedStart <= germline_collection_age        ~ "regimen before germline",
    AgeAtMedStart > germline_collection_age         ~ "regimen after germline",
    is.na(AgeAtMedStart) &
      !is.na(Disease.Type)                          ~ "no age"
  )) %>% 
  #restrict to drugs before germline
  filter(regimen_specimen_sequence != "regimen after germline" |
           is.na(regimen_specimen_sequence)) %>% 
  select(-c(Disease.Type, germline_collection_age, regimen_specimen_sequence))


library(data.table)
Medications <- dcast(setDT(Medications1),
                     AvatarKey + MedicationInd ~ rowid(AvatarKey),
                     value.var = c("AgeAtMedStart", "regimen", "AgeAtMedStop")) %>%
  # 1 patient still has 200 rows 
  select(AvatarKey, drugs_ever = MedicationInd, starts_with("AgeAtMedStart_"),
         starts_with("regimen_"), starts_with("AgeAtMedStop_")) %>%
  # select(AvatarKey, drugs_ever = MedicationInd, AgeAtMedStart_1 : AgeAtMedStart_20,
  #        Medication_1 : Medication_20, AgeAtMedStop_1 : AgeAtMedStop_20) %>%
  mutate(has_medication_data = "Yes")

# rm(Medications1, Medications2)
# write_rds(Medications,
#           paste0(here::here(), 
#                  "/data/processed data",
#                  "/Medication wide format_",
#                  today(), ".rds"))

# sct----
StemCellTransplant1 <- StemCellTransplant %>%
  # For the patient with a diagnosis site
  inner_join(., Diagnosis %>%
               select(AvatarKey, PrimaryDiagnosisSiteCode, PrimaryDiagnosisSite),
             by = c("AvatarKey", "SCTPrimaryDiagnosisSiteCode" = "PrimaryDiagnosisSiteCode",
                    "SCTPrimaryDiagnosisSite" = "PrimaryDiagnosisSite")) %>%
  arrange(AvatarKey, AgeAtTransplant) %>%
  distinct(AvatarKey, .keep_all = TRUE)

StemCellTransplant2 <- StemCellTransplant %>%
  # For the patients with no diagnosis site
  filter(SCTPrimaryDiagnosisSite == "Unknown/Not Applicable") %>%
  arrange(AvatarKey, AgeAtTransplant) %>%
  distinct(AvatarKey, .keep_all = TRUE)

StemCellTransplant <- StemCellTransplant1 %>%
  # bind the data with a primary site known as first diagnosis
  # then the one without which show some No Metastasis
  bind_rows(., StemCellTransplant2) %>%
  select(AvatarKey, sct_ever = SCTInd, 
         age_at_first_sct = AgeAtTransplant,
         TransplantType, TransplantCellSource) %>%
  mutate(across(c("age_at_first_sct"), ~ case_when(
    . == "Age 90 or older"                ~ 90,
    . == "Unknown/Not Applicable"         ~ NA_real_,
    TRUE                                  ~ as.numeric(.)
  ))) %>%
  distinct(AvatarKey, .keep_all = TRUE) %>%
  mutate(has_sct_data = "Yes")

rm(StemCellTransplant1, StemCellTransplant2)


# radiation}
Radiation <- Radiation %>%
  filter(RadiationTherapyInd != "(Migrated) Cannot determine from available documentation")

Radiation1 <- Radiation %>%
  # For the patient with a diagnosis site
  inner_join(., Diagnosis %>%
               select(AvatarKey, PrimaryDiagnosisSiteCode, PrimaryDiagnosisSite),
             by = c("AvatarKey", "RadPrimaryDiagnosisSiteCode" = "PrimaryDiagnosisSiteCode",
                    "RadPrimaryDiagnosisSite" = "PrimaryDiagnosisSite")) %>%
  arrange(AvatarKey, AgeAtRadiationStart) %>%
  distinct(AvatarKey, .keep_all = TRUE)

Radiation2 <- Radiation %>%
  # For the patients with no diagnosis site
  filter(RadPrimaryDiagnosisSite == "Unknown/Not Applicable") %>%
  arrange(AvatarKey, AgeAtRadiationStart) %>%
  distinct(AvatarKey, .keep_all = TRUE)

Radiation <- Radiation1 %>%
  # bind the data with a primary site known as first diagnosis
  # then the one without which show some No Metastasis
  bind_rows(., Radiation2) %>%
  select(AvatarKey, RadPrimaryDiagnosisSiteCode, RadPrimaryDiagnosisSite,
         radiation_ever = RadiationTherapyInd, 
         age_at_first_radiation = AgeAtRadiationStart) %>%
  mutate(across(c("age_at_first_radiation"), ~ case_when(
    . == "Age 90 or older"                ~ 90,
    . == "Unknown/Not Applicable"         ~ NA_real_,
    TRUE                                  ~ as.numeric(.)
  ))) %>%
  distinct(AvatarKey, .keep_all = TRUE) %>%
  mutate(has_radiation_data = "Yes")

rm(Radiation1, Radiation2)


# surgery ---- 
################################################################################################################# Done
SurgeryBiopsy <- SurgeryBiopsy %>%
  filter(SurgeryBiopsyInd != "(Migrated) Cannot determine from available documentation")

SurgeryBiopsy1 <- SurgeryBiopsy %>%
  # For the patient with a diagnosis site
  inner_join(., Diagnosis %>%
               select(AvatarKey, PrimaryDiagnosisSiteCode, PrimaryDiagnosisSite),
             by = c("AvatarKey", "PrimaryDiagnosisSiteCode",
                    "PrimaryDiagnosisSite")) %>%
  arrange(AvatarKey, AgeAtSurgeryBiopsy) %>%
  distinct(AvatarKey, .keep_all = TRUE)

SurgeryBiopsy2 <- SurgeryBiopsy %>%
  # For the patients with no diagnosis site
  filter(PrimaryDiagnosisSite == "Unknown/Not Applicable") %>%
  arrange(AvatarKey, AgeAtSurgeryBiopsy) %>%
  distinct(AvatarKey, .keep_all = TRUE)

SurgeryBiopsy <- SurgeryBiopsy1 %>%
  # bind the data with a primary site known as first diagnosis
  # then the one without which show some No Metastasis
  bind_rows(., SurgeryBiopsy2) %>%
  select(AvatarKey, SurgPrimaryDiagnosisSiteCode = PrimaryDiagnosisSiteCode,
         SurgPrimaryDiagnosisSite = PrimaryDiagnosisSite,
         surgerybiopsy_ever = SurgeryBiopsyInd,
         surgery_ever = SiteTherapeutic, 
         age_at_first_surgerybiopsy = AgeAtSurgeryBiopsy) %>%
  mutate(across(c("age_at_first_surgerybiopsy"), ~ case_when(
    . == "Age 90 or older"                ~ 90,
    . == "Unknown/Not Applicable"         ~ NA_real_,
    TRUE                                  ~ as.numeric(.)
  ))) %>%
  distinct(AvatarKey, .keep_all = TRUE) %>%
  mutate(surgery_ever = case_when(
    surgery_ever == "Unknown/Not Applicable" &
      !is.na(age_at_first_surgerybiopsy)  ~ "Yes",
    surgery_ever == "Unknown/Not Applicable" &
      is.na(age_at_first_surgerybiopsy)   ~ "No",
    TRUE                                  ~ surgery_ever
  )) %>%
  mutate(has_surgery_data = "Yes")

rm(SurgeryBiopsy1, SurgeryBiopsy2)


# physical}
PhysicalAssessment <- PhysicalAssessment %>%
  left_join(., Diagnosis %>%
              select(AvatarKey, AgeAtDiagnosis),
            by = c("AvatarKey")) %>%
  mutate(across(c("AgeAtPhysicalExam"), ~ case_when(
    . == "Age 90 or older"                ~ 90,
    . == "Unknown/Not Applicable"         ~ NA_real_,
    TRUE                                  ~ as.numeric(.)
  ))) %>%
  mutate(int = abs(AgeAtPhysicalExam - AgeAtDiagnosis)) %>%
  arrange(AvatarKey, int) %>%
  distinct(AvatarKey, .keep_all = TRUE) %>%
  select(AvatarKey, AgeAtPhysicalExam, BodyWeight, BodyHeight, BMI) %>%
  mutate(across(c("BodyWeight", "BodyHeight", "BMI"), ~ as.numeric(.)
  )) %>%
  mutate(bmi_cat = case_when(
    BMI < 18.5                  ~ "Underweight",
    BMI >= 18.5 &
      BMI < 25                  ~ "Healthy",
    BMI >= 25 &
      BMI < 30                  ~ "Overweight",
    BMI >= 30                   ~ "Obese"
  )) %>%
  mutate(bmi_cat = factor(bmi_cat, levels = c("Underweight", "Healthy", "Overweight", "Obese")))


# tumor markers}
TumorMarker <- TumorMarker %>%
  mutate(across(c("AgeAtTumorMarkerTest"), ~ case_when(
    . == "Age 90 or older"                ~ 90,
    . == "Unknown/Not Applicable"         ~ NA_real_,
    TRUE                                  ~ as.numeric(.)
  ))) %>%
  mutate(across(where(is.character), ~ na_if(., "Unknown/Not Applicable"))) %>%
  rename(marker_test = TMarkerTest, marker_unit = TMarkerValueUOM,
         marker_value = TMarkerResultValue) %>%
  group_by(marker_test) %>%
  fill(TMarkerLowRange, TMarkerHighRange, .direction = "updown") %>%
  ungroup() %>%
  mutate(across(c("marker_value", "TMarkerLowRange",
                  "TMarkerHighRange"), ~ as.numeric(.)
  )) %>%
  mutate(marker_interpretation = case_when(
    marker_value >= TMarkerLowRange &
      marker_value <= TMarkerHighRange     ~ "WNL (Within Normal Limits)",
    marker_value < TMarkerLowRange         ~ "Low",
    marker_value > TMarkerHighRange        ~ "High"
  ), marker_interpretation = coalesce(TMarkerRangeIndicator, marker_interpretation)) %>%
  mutate(marker_interpretation = case_when(
    TMarkerResult == "Value"              ~ marker_interpretation,
    TRUE                                  ~ TMarkerResult
  )) %>%
  mutate(marker_test = case_when(
    str_detect(marker_test, "HER2")       ~ "HER2",
    TRUE                                  ~ marker_test
  )) %>% 
  group_by(AvatarKey, marker_test) %>%
  mutate(er_result = case_when(
    marker_test == "Estrogen receptor (ER)"        ~ first(TMarkerResult)
  )) %>%
  mutate(pr_result = case_when(
    marker_test == "Progesterone receptor (PR)"    ~ first(TMarkerResult)
  )) %>%
  mutate(clean_her2 = case_when(
    str_detect(marker_test, "HER2") & 
      (TMarkerResult == "Positive" | 
         TMarkerResult == "Negative")                ~ TMarkerResult,
    str_detect(marker_test, "HER2") & 
      TMarkerResult == "Amplified"                 ~ "Positive"
  ), .after = marker_test) %>% 
  mutate(her_result = case_when(
    str_detect(marker_test, "HER2")                ~ first(clean_her2)
  ), .after = marker_test) %>% 
  mutate(ER_PR_HER2_test = case_when(
    marker_test == "Estrogen receptor (ER)" |
      marker_test == "Progesterone receptor (PR)" |
      marker_test == "HER2"                        ~ "ER_PR_HER2_status"
  )) %>% 
  mutate(er_receptor = case_when(
    er_result == "Positive"         ~ "ER+",
    er_result == "Negative"         ~ "ER-"
  )) %>% 
  mutate(pr_receptor = case_when(
    pr_result == "Positive"         ~ "PR+",
    pr_result == "Negative"         ~ "PR-"
  )) %>% 
  mutate(her_receptor = case_when(
    her_result == "Positive"         ~ "HER2+",
    her_result == "Negative"         ~ "HER2-"
  )) %>% 
  group_by(AvatarKey, ER_PR_HER2_test) %>% 
  fill(er_receptor, pr_receptor, her_receptor, .direction = "downup") %>% 
  ungroup() %>% 
  unite(ER_PR_HER2_status, c(er_receptor, pr_receptor, her_receptor), sep = "/", remove = FALSE, na.rm = TRUE) %>% 
  # TMarkerPercentStainResultValue is not helpful
  mutate(marker_test = coalesce(ER_PR_HER2_test, marker_test)) %>% 
  mutate(marker_interpretation = coalesce(ER_PR_HER2_status, marker_interpretation)) %>% 
  select(AvatarKey, AgeAtTumorMarkerTest, marker_test,
         marker_interpretation, marker_value, marker_unit)

TumorMarker <- TumorMarker %>%
  arrange(AvatarKey, marker_test, AgeAtTumorMarkerTest) %>%
  distinct(AvatarKey, marker_test, .keep_all = TRUE) %>%
  pivot_wider(id_cols = AvatarKey,
              names_from = marker_test,
              values_from = marker_interpretation)


# # MSI
# MSI_marker <- MSI_marker %>% 
#   mutate(MSI_high_score = case_when(
#     `MSIsensor2 Score` >= 20       ~ "Yes",
#     `MSIsensor2 Score` < 20        ~ "No"
#   )) %>% 
#   select(ORIENAvatarKey, `WES SLID`, MSI_high_score)


# family history
FamilyHistory <- FamilyHistory %>%
  arrange(AvatarKey, desc(CancerInFamilyInd)) %>%
  distinct(AvatarKey, .keep_all = TRUE) %>%
  select(AvatarKey, CancerInFamilyInd)


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
all_tumor_data <- ch_calls %>%
  full_join(., demographics, by = "AvatarKey") %>%
  full_join(., PatientHistory %>%
              select(AvatarKey, SmokingStatus),
            by = "AvatarKey") %>%
  full_join(., VitalStatus, by = "AvatarKey") %>%
  full_join(., Diagnosis, by = "AvatarKey") %>%
  full_join(., CytogeneticAbnormalities, by = "AvatarKey") %>%
  full_join(., TumorMarker, by = "AvatarKey") %>%
  full_join(., PhysicalAssessment, by = "AvatarKey") %>%
  full_join(., Medications, by = "AvatarKey") %>%
  full_join(., StemCellTransplant, by = "AvatarKey") %>%
  full_join(., Radiation, by = "AvatarKey") %>%
  full_join(., SurgeryBiopsy, by = "AvatarKey") %>%
  full_join(., Outcomes, by = "AvatarKey") %>%
  full_join(., MetastaticDisease, by = "AvatarKey") %>%
  full_join(., last_date, by = "AvatarKey") %>%
  full_join(., FamilyHistory, by = "AvatarKey")

rm(demographics, VitalStatus, Diagnosis,
   CytogeneticAbnormalities, TumorMarker,
   PhysicalAssessment, Medications,
   StemCellTransplant, Radiation,
   SurgeryBiopsy, MetastaticDisease,
   Outcomes, last_date, FamilyHistory,
   ClinicalMolLinkage, ch_calls)


all_tumor_data <- all_tumor_data %>% 
  # Age at last contact
  mutate(is_agelastcontact_last_date = case_when(
    AgeAtLastContact >= age_at_lab                  ~ "Yes",
    AgeAtLastContact < age_at_lab                   ~ "No"
  )) %>% 
  mutate(age_last_contact = case_when(
    is_agelastcontact_last_date == "Yes"            ~ AgeAtLastContact,
    is_agelastcontact_last_date == "No"             ~ age_at_lab
  )) %>% 
  # Age at first treatment
  mutate(first_treatment = case_when(
    (age_at_first_radiation < age_at_first_surgerybiopsy |
       (is.na(age_at_first_surgerybiopsy) &
          !is.na(age_at_first_radiation))) &
      
      (age_at_first_radiation < AgeAtMedStart_1 |
         (is.na(AgeAtMedStart_1) &
            !is.na(age_at_first_radiation))) &
      
      (age_at_first_radiation < age_at_first_sct | 
         (is.na(age_at_first_sct) &
            !is.na(age_at_first_radiation)))              ~ "Radiation",
    
    (age_at_first_surgerybiopsy < age_at_first_radiation | 
       (is.na(age_at_first_radiation) &
          !is.na(age_at_first_surgerybiopsy))) &
      
      (age_at_first_surgerybiopsy < AgeAtMedStart_1 |
         (is.na(AgeAtMedStart_1) &
            !is.na(age_at_first_surgerybiopsy))) &
      
      (age_at_first_surgerybiopsy < age_at_first_sct | 
         (is.na(age_at_first_sct) &
            !is.na(age_at_first_surgerybiopsy)))          ~ "Surgery",
    
    (AgeAtMedStart_1 < age_at_first_radiation | 
       (is.na(age_at_first_radiation) &
          !is.na(AgeAtMedStart_1))) &
      
      (AgeAtMedStart_1 < age_at_first_surgerybiopsy |
         (is.na(age_at_first_surgerybiopsy) &
            !is.na(AgeAtMedStart_1))) &
      
      (AgeAtMedStart_1 < age_at_first_sct | 
         (is.na(age_at_first_sct) &
            !is.na(AgeAtMedStart_1)))             ~ "Drugs",
    
    (age_at_first_sct < age_at_first_radiation | 
       (is.na(age_at_first_radiation) &
          !is.na(age_at_first_sct))) &
      
      (age_at_first_sct < age_at_first_surgerybiopsy |
         (is.na(age_at_first_surgerybiopsy) &
            !is.na(age_at_first_sct))) &
      
      (age_at_first_sct < AgeAtMedStart_1 |
         (is.na(AgeAtMedStart_1) &
            !is.na(age_at_first_sct)))             ~ "SCT"
    
  )) %>% 
  mutate(age_at_first_treatment = case_when(
    first_treatment == "Radiation"                  ~ age_at_first_radiation,
    first_treatment == "Surgery"                    ~ age_at_first_surgerybiopsy,
    first_treatment == "Drugs"                      ~ AgeAtMedStart_1,
    first_treatment == "SCT"                        ~ age_at_first_sct,
  )) %>% 
  # OS
  mutate(os_event = case_when(
    VitalStatus == "Alive"                          ~ 0,
    VitalStatus == "Lost to follow-up"              ~ 0,
    VitalStatus == "Dead"                           ~ 1
  )) %>% 
  mutate(os_age = coalesce(AgeAtDeath, age_last_contact)) %>% 
  mutate(os_time_from_dx_years = os_age - AgeAtDiagnosis) %>% 
  mutate(os_time_from_treatment_years = os_age - age_at_first_treatment) %>% 
  # PFS
  mutate(pfs_event = case_when(
    ProgRecurInd == "No"                            ~ 0,
    ProgRecurInd == "Progression"                   ~ 1,
    ProgRecurInd == "Recurrence"                    ~ 1
  )) %>% 
  mutate(pfs_age = coalesce(AgeAtProgRecur, AgeAtLastContact)) %>% 
  mutate(pfs_time_from_dx_years = pfs_age - AgeAtDiagnosis) %>% 
  mutate(pfs_time_from_treatment_years = pfs_age - age_at_first_treatment) %>% 
  # Metastasis
  mutate(metastase_event = case_when(
    had_metastasis == "Yes"                         ~ 1,
    !is.na(has_metastasis_data)                     ~ 0
  )) %>% 
  mutate(met_age = coalesce(AgeAtMetastaticSite, AgeAtLastContact)) %>% 
  mutate(met_time_from_dx_years = met_age - AgeAtDiagnosis) %>% 
  mutate(met_time_from_treatment_years = met_age - age_at_first_treatment) %>% 
  mutate(time_dx_to_first_treatment = age_at_first_treatment - AgeAtDiagnosis)


# Save
write_csv(all_tumor_data,
          paste0(here::here(), 
                 "/data/processed data",
                 "/all_tumor_data_",
                 today(), ".csv"))
write_rds(all_tumor_data, 
          paste0(here::here(), 
                 "/data/processed data",
                 "/all_tumor_data_",
                 today(), ".rds"))

path_save <- fs::path("", "Volumes", "Gillis_Research",
                      "Christelle Colin-Leitzinger", "Ovary CH",
                      "ovary_ch")
write_csv(all_tumor_data, 
          paste0(path_save, 
                 "/data/processed data",
                 "/all_tumor_data_",
                 today(), ".csv"))

path_save <- fs::path("", "Volumes", "Peres_Research",
                      "ORIEN analysis", "Ovary CH")
write_csv(all_tumor_data, 
          paste0(path_save, 
                 "/data/processed data",
                 "/all_tumor_data_",
                 today(), ".csv"))

gyn_data <- all_tumor_data %>% 
  filter(Disease.Type == "GYN - Ovarian Cancer")

write_csv(gyn_data, 
          paste0(path_save, 
                 "/data/processed data",
                 "/ORIEN_data_with_ovarian_sample",
                 today(), ".csv"))

write_csv(gyn_data,
          paste0(here::here(), 
                 "/data/processed data",
                 "/ORIEN_data_with_ovarian_sample",
                 today(), ".csv"))
write_rds(gyn_data,
          paste0(here::here(), 
                 "/data/processed data",
                 "/ORIEN_data_with_ovarian_sample",
                 today(), ".rds"))


# End data cleaning

