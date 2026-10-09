# Import Library
library(tidyverse)
library(lubridate)


################################################################################# I ### Load data
# # Load each file in its own dataframe
path_orien <- fs::path("", "Volumes", "Gillis_Research",
                     "Lab_Data", "CHinORIEN")

load(paste0(dirname(here::here()),
            "/ch_orien/data/raw_data",
            # path_orien,
            # "/RawData/23PRJ127MCC_NormalizedFiles_20260730",
            "/CHinORIEN_AllTumorTypes_Files_20260730.RData"))

# path <- fs::path("", "Volumes", "Gillis_Research",
#                      "Lab_Data")

path_gillis <- fs::path("", "Volumes", "Gillis_Research", "Lab_Data", "CHinOvary")
# ch_calls <- read_csv(paste0(
#   path,
#   "/ProcessedData",
#   "/CHinOvary_Samplelist_CHstatus_07.14.2025.csv"))
ch_calls <- read_csv(paste0(here::here(),
                            "/data/raw_data",
                            # path_gillis,
                            # "/ProcessedData",
                            "/CHinOvary_Samplelist_CHstatus_07.14.2025.csv"))

# ClinicalMolLinkage file in Normalized folder doesn't have WESid
# Use the file from first folder
old_ClinicalMolLinkage <- read_csv(paste0(here::here(),
                                      "/data/raw_data",
                                      # dirname(path_gillis), "/CHinORIEN",
                                      # "/RawData",
                                      # "/23PRJ127MCC_NormalizedFiles_20260730",
                                      "/OldVersion/CHinOvary_ClinicalMolLinkage_V4.csv"))

# I updated the drug_class file several times # Not run
# path_drug_class <- fs::path("", "Volumes", "Gillis_Research", "Lab_Data", "BreastCHEvolution")
# oldversion_drug_class <- 
#   readxl::read_xlsx(paste0(
#     path_drug_class, 
#     "/ProcessedData/BreastCH_BoltonChemoDosingNatGenet2020.xlsx"),
#     sheet = "Sup. Table 1", skip = 1) %>% 
#   janitor::clean_names()
# Load updated file
path_drug <- fs::path("", "Volumes", "Gillis_Research",
                 "Lab_Data", "SharedResources")
drug_class <- 
  read.csv(paste0(#here::here(), 
                  # "/data/processed_data",
                  path_drug, "/BoltonDrugCategories",
                  "/CHevolution_Updated_BoltonChemoDosing_20260713.csv"))
removed_drug <- 
  read.csv(paste0(#here::here(), 
                  # "/data/processed_data",
                  path_drug, "/NotCancerDrugs",
                  "/CHinOvary_NotCancerDrug_ToRemovedFromData_20260617.csv"))


################################################################################# II ### Data cleaning
# CH calle----
ch_calls <- ch_calls %>%
  rename(AvatarKey = ORIENAvatarKey) %>%
  select(-`...1`)

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
  old_ClinicalMolLinkage %>%
  filter(Tumor.Germline == "Tumor") %>%
  # left_join(., MSI_marker,
  #           by = c("ORIENAvatarKey", "WES" = "WES SLID")) %>% 
  select(ORIENAvatarKey, Tumor_WES = WES,
         tumor_collection_age = Age.At.Specimen.Collection,
         tumor_collection_disease_type= Disease.Type,
         tumor_SpecimenSiteOfCollection = SpecimenSiteOfCollection#,
         # MSI_high_score
  )
Germline_wes <-
  old_ClinicalMolLinkage %>%
  filter(Tumor.Germline == "Germline") %>%
  select(ORIENAvatarKey, Germline_WES = WES,
         germline_collection_age = Age.At.Specimen.Collection,
         germline_SpecimenSiteOfCollection = SpecimenSiteOfCollection)

# this create a germline-tumor file with all combination germline-tumor
# We selected the good combination to do the sequencing
# So the correct combination will be selected in the next merge
ClinicalMolLinkage1 <- full_join(Germline_wes, Tumor_wes,
                                by = "ORIENAvatarKey")

ch_calls <- ch_calls %>%
  left_join(., ClinicalMolLinkage1,
            by = c("AvatarKey" = "ORIENAvatarKey", 
                   "Germline_WES", "Tumor_WES", "Disease.Type" = "tumor_collection_disease_type"
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
  mutate(race_all = case_when(
    str_detect(Race, "American Indian")        ~ "American Indian or Alaska Native",
    Race == "White"                            ~ "White",
    Race == "Black or African American"        ~ "Black",
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
      Race == "Native Hawaiian" |
      Race == "Micronesian, NOS" |
      Race == "Pacific Islander, NOS" |
      Race == "Polynesian, NOS" |
      Race == "Samoan" |
      str_detect(Race, "Tongan")               ~ "Native Hawaiian or Other Pacific Islander",
    str_detect(Race, "Unknown")                ~ NA_character_,
    Race == "Some other race"                  ~ "Other",
    TRUE                                       ~ Race
  )) %>%
  mutate(race = case_when(
    race_all == "White"                                              ~ "White",
    race_all == "Black"                                              ~ "Black",
    race_all == "Asian"                                              ~ "Asian",
    !is.na(race_all)                                                 ~ "Other"
  )) |> 
  mutate(Ethnicity = case_when(
    Ethnicity == "Spanish surname only"        ~ "Non-Hispanic",
    str_detect(Ethnicity, "Non-Hispanic")      ~ "Non-Hispanic",
    str_detect(Ethnicity, "Unknown")           ~ NA_character_,
    !is.na(Ethnicity)                          ~ "Hispanic",
    TRUE                                       ~ Ethnicity
  )) %>%
  mutate(raceeth = case_when(
    race == "White" &
      Ethnicity == "Non-Hispanic"              ~ "White Non-Hispanic",
    race == "Black" &
      Ethnicity == "Non-Hispanic"              ~ "Black Non-Hispanic",
    race == "Others" &
      Ethnicity == "Non-Hispanic"              ~ "Others Non-Hispanic",
    Ethnicity == "Hispanic"                    ~ "Hispanic, any race",
    TRUE                                       ~ paste(Race, Ethnicity)
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
  ), SmokingStatus = factor(SmokingStatus, levels = c("Never", "Ever"))) |> 
  mutate(AlcoholUse = case_when(
    AlcoholUse == "Never"                                            ~ "Never",
    AlcoholUse == "Current" |
      str_detect(AlcoholUse, "Ever")                                 ~ "Ever"
  ), AlcoholUse = factor(AlcoholUse, levels = c("Never", "Ever")))


# Diagnosis----
Diagnosis_save <- Diagnosis
Diagnosis <- Diagnosis_save


Diagnosis %>%
  filter(str_detect(AvatarKey, paste0(ch_calls$AvatarKey, collapse = "|"))) %>% 
  select(PrimaryDiagnosisSite, PrimaryDiagnosisSiteCode) %>%
  distinct() %>%
  arrange(PrimaryDiagnosisSiteCode, PrimaryDiagnosisSite)

Diagnosis1 <- Diagnosis %>%
  filter(str_detect(AvatarKey, paste0(ch_calls$AvatarKey, collapse = "|"))) %>% 
  group_by(AvatarKey) %>% 
  mutate(number_of_dx = n(), .before = AgeAtDiagnosis) %>% 
  ungroup() |> 
  
  # Few missing age at dx
  # It is ok to use the age at first contact to fill missing age at dx
  mutate(not_real_dxage = case_when(
    AgeAtDiagnosis == "Age 90 or older"   ~ "Age 90 or older"
  )) %>%
  mutate(across(c("AgeAtDiagnosis", 
                  "AgeAtFirstContact"), ~ case_when(
                    . == "Age 90 or older"                ~ 90,
                    . == "Unknown/Not Applicable"         ~ NA_real_,
                    TRUE                                  ~ as.numeric(.)
                  ))) %>%
  mutate(dx_age_is_first_contact_age = case_when(
    !is.na(AgeAtDiagnosis)                ~ "Age is from AgeAtDiagnosis",
    is.na(AgeAtDiagnosis) &
      !is.na(AgeAtFirstContact)           ~ "Age is from AgeAtFirstContact"
  ), .after = AgeAtDiagnosis) |> 
  mutate(AgeAtDiagnosis = coalesce(AgeAtDiagnosis, AgeAtFirstContact)) %>% 
  arrange(AvatarKey, AgeAtDiagnosis) %>% 
  # Create var specific to ovarian dx
  group_by(AvatarKey) %>% 
  mutate(diagnosis_sequence = row_number(AvatarKey), .after = AvatarKey) %>%
  ungroup() %>%
  # Find ovarian cancer
  left_join(., ch_calls |> 
              select(AvatarKey, Histology.Behavior), 
            by = c("AvatarKey")) |> 
  mutate(is_ovarian_cancer = case_when(
    PrimaryDiagnosisSite == "Ovary" & 
      (Histology == "Serous carcinoma, NOS" |
         Histology == "Endometrioid adenocarcinoma, NOS" |
         Histology == "Low-grade serous carcinoma" |
         Histology == "High-grade serous carcinoma" |
         Histology == "Clear cell adenocarcinoma, NOS" |
         Histology == "Mixed cell adenocarcinoma" |
         Histology == "Carcinosarcoma, NOS" |
         Histology == "Carcinoma, NOS" |
         Histology == "Mucinous adenocarcinoma" |
         Histology == "Mucinous cystadenocarcinoma, NOS" |
         Histology == "Mullerian mixed tumor" |
         Histology == "Squamous cell carcinoma, NOS" |
         Histology == "Transitional cell carcinoma, NOS" |
         Histology == "Carcinoma, anaplastic, NOS" |
         Histology == "Clear cell adenocarcinofibroma" |
         Histology == "Large cell carcinoma, NOS" |
         Histology == "Mesodermal mixed tumor" |
         Histology == "Papillary adenocarcinoma, NOS" |
         Histology == "Small cell carcinoma, NOS")                ~ "Yes",
    PrimaryDiagnosisSite == "Fallopian tube" & 
      (Histology == "Serous carcinoma, NOS" |
         Histology == "High-grade serous carcinoma" |
         Histology == "Low-grade serous carcinoma" |
         Histology == "Carcinosarcoma, NOS" |
         Histology == "Mullerian mixed tumor" |
         Histology == "Papillary adenocarcinoma, NOS")            ~ "Yes",
    PrimaryDiagnosisSite == "Ovary" & 
      Histology == "Adenocarcinoma, NOS" &
      (Histology.Behavior == "84413 Serous cystadenocarcinoma, NOS" |
         Histology.Behavior == "81403 Adenocarcinoma, NOS")       ~ "Yes",
    PrimaryDiagnosisSite == "Peritoneum, NOS" & 
      (Histology == "Serous carcinoma, NOS" |
         Histology == "High-grade serous carcinoma" |
         Histology == "Adenocarcinoma, NOS" |
         Histology == "Endometrioid adenocarcinoma, NOS" |
         Histology == "Low-grade serous carcinoma")               ~ "Yes",
    PrimaryDiagnosisSite == "Specified parts of peritoneum" & 
      Histology == "Serous carcinoma, NOS"                        ~ "Yes"
  ), .after = Histology) |> 
  select(-Histology.Behavior) |> 
  group_by(AvatarKey, is_ovarian_cancer) |> 
  mutate(number_of_ovarian_dx = case_when(
    is_ovarian_cancer == "Yes"                                    ~ n()
  ), .before = AgeAtDiagnosis) %>% 
  ungroup()

Diagnosis_ovarian <- Diagnosis1 %>%
  # Keep the first ovarian cancer
  # Will create a cancer ID so we can identify the pre and post cancer
  filter(is_ovarian_cancer == "Yes") |> 
  distinct(AvatarKey, .keep_all = TRUE) |> 
  select(AvatarKey, number_of_dx, diagnosis_sequence,
         AgeAtFirstContact,
         AgeAtDiagnosis, not_real_dxage, AgeAtDiagnosisFlag, dx_age_is_first_contact_age,
         YearOfDiagnosis,
         PrimaryDiagnosisSiteCode : Histology,
         ClinGroupStage, PathGroupStage,
         GradeClinical, GradePathological,
         CurrentlySeenForPrimaryOrRecurr,
         PerformStatusAtDiagnosis, 
         OtherStagingSystem, OtherStagingValue, everything(),
         -c(is_ovarian_cancer, AgeAtFirstContactFlag)) |> 
  unite(ovarian_dx_id, c("AvatarKey", "diagnosis_sequence"), sep = "_", remove = FALSE)

pre_post_diagnosis <- Diagnosis1 %>%
  unite(ovarian_dx_id, c("AvatarKey", "diagnosis_sequence"), sep = "_", remove = FALSE) %>%
  # Merge by the ovarian dx id I created
  full_join(., Diagnosis_ovarian |> 
              select(ovarian_dx_id, diagnosis_sequence2 = diagnosis_sequence) |> 
              mutate(is_actual_ovarian_dx = "Yes"),
            by = "ovarian_dx_id") |> 
  # create var for pre and post dx
  group_by(AvatarKey) |> 
  fill(diagnosis_sequence2, .direction = "updown") |> 
  ungroup() |> 
  mutate(pre_or_post = case_when(
    diagnosis_sequence < diagnosis_sequence2           ~ "pre",
    diagnosis_sequence > diagnosis_sequence2           ~ "post",
  )) |> 
  # Remove first ovarian cancer from the pre and post
  filter(pre_or_post == "pre" | pre_or_post == "post") |> 
  select(AvatarKey, 
         AgeAtDiagnosis, PrimaryDiagnosisSite, Histology, pre_or_post) |> 
  unite(cancer_info_seperated_byslash, c("AgeAtDiagnosis" : "Histology"), sep = "; ", remove = TRUE) |> 
  group_by(AvatarKey, pre_or_post) |> 
  summarise_at(vars(cancer_info_seperated_byslash), str_c, collapse = " / ") |> 
  ungroup() |> 
  pivot_wider(id_cols = AvatarKey, 
              names_from = "pre_or_post", 
              values_from = cancer_info_seperated_byslash, 
              names_glue = "{pre_or_post}_{.value}")

Diagnosis <- Diagnosis_ovarian %>%
  full_join(., pre_post_diagnosis,
            by = "AvatarKey") |> 
  select(-ovarian_dx_id) |> 
  mutate(ECOG = str_match(PerformStatusAtDiagnosis, "ECOG ([:digit:])")[,2], 
         .after = PerformStatusAtDiagnosis) %>%
  mutate(Karnofsky = str_match(PerformStatusAtDiagnosis, "Karnofsky ([:digit:].*)%")[,2], 
         .after = PerformStatusAtDiagnosis) %>%
  # Stage
  mutate(ClinGroupStage = case_when(
    str_detect(ClinGroupStage, "IV")                        ~ "4",
    str_detect(ClinGroupStage, "III")                       ~ "3",
    str_detect(ClinGroupStage, "II")                        ~ "2",
    str_detect(ClinGroupStage, "I")                         ~ "1",
    str_detect(ClinGroupStage, "0")                         ~ "0"
  )) %>% 
  mutate(PathGroupStage = case_when(
    str_detect(PathGroupStage, "IV")                        ~ "4",
    str_detect(PathGroupStage, "III")                       ~ "3",
    str_detect(PathGroupStage, "II")                        ~ "2",
    str_detect(PathGroupStage, "I")                         ~ "1",
    str_detect(PathGroupStage, "0")                         ~ "0"
  )) %>% 
  group_by(AvatarKey) %>% 
  mutate(stage = max(ClinGroupStage, PathGroupStage, na.rm = TRUE), .after = PathGroupStage) %>%
  ungroup()


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


# Medication----
Medications_ <- Medications %>%
  filter(str_detect(AvatarKey, paste0(Diagnosis$AvatarKey, collapse = "|"))) %>% 
  mutate(across(c("AgeAtMedStart",
                  "AgeAtMedStop"), ~ case_when(
                    . == "Age 90 or older"                ~ 90,
                    . == "Unknown/Not Applicable"         ~ NA_real_,
                    TRUE                                  ~ as.numeric(.)
                  ))) %>% 
  group_by(AvatarKey) %>% 
  mutate(ever_first_med_age = min(AgeAtMedStart, AgeAtMedStop)) %>% 
  ungroup() |> 
  # Here the script to update names to fit Bolton categories
  mutate(Medication = str_remove(
    Medication, 
    " Hydrochloride|Liposomal | Sulfate| Phosphate| Citrate| Acetate| Camsylate| Disodium| Mesylate| Ditosylate| Tosylate| Tartrate| \\(Leucovorin\\)"), 
    Medication = case_when(
      str_detect(Medication, "Paclitaxel")        ~ "Paclitaxel",
      str_detect(Medication, "^Interferon")       ~ "Interferon",
      Medication == "Bevacizumab-maly"            ~ "Bevacizumab",
      TRUE                                        ~ Medication
    ))

# Now that the drug names are updated, I can join the data with the Bolton categories
Medications1 <- Medications_ %>% 
  mutate(Medication = str_to_lower(Medication)) %>% 
  # Remove non cancer drugs from data so their date start will not count
  left_join(., removed_drug, by = c("Medication" = "drug_name")) %>% 
  filter(is.na(drug_type)) %>% select(-drug_type) %>% 
  # Add Bolton chemo cat
  left_join(., drug_class, by = c("Medication" = "drug_name")) %>% 
  # and code an overall chemotherapy variable
  mutate(is_chemotherapy = case_when(
    !is.na(narrow_drug_class_cytotoxic_only)        ~ "Yes",
    Medication == "chemo, nos"                      ~ "Yes"
  ), .after = MedicationInd) 

# Separate patient who didn't receive drugs
Medications_never <- Medications1 %>% 
  filter(MedicationInd == "No") %>% 
  mutate(has_first_line = "Never received any drug", .after = MedicationInd) %>%
  mutate(is_chemotherapy =  "Never received any drug", .after = MedicationInd)

Medications_yes <- Medications1 %>% 
  # Filter patients who received drugs
  filter(MedicationInd == "Yes") |> 
  # Fix age for a couple of patients
  mutate(AgeAtMedStart = coalesce(AgeAtMedStart, AgeAtMedStop)) %>% 
  group_by(AvatarKey) %>% 
  mutate(ever_first_med_age = min(AgeAtMedStart, na.rm = TRUE),
         ever_first_med_age = na_if(ever_first_med_age, Inf)) %>% 
  ungroup() %>% 
  # recode line
  # There are regimen line are unknown 
  # but some have the same age as other row for which the line is known, use to fill it up
  mutate(MedLineRegimen = case_when( # m
    MedLineRegimen %in% c("Unknown/Not Applicable",
                          "Unknown/Not Reported")    ~ NA_character_,
    TRUE                                             ~ MedLineRegimen
  )) %>%
  group_by(AvatarKey, AgeAtMedStart) %>% 
  fill(MedLineRegimen, .direction = "updown") %>% 
  ungroup() %>% 
  # Do the reverse
  group_by(AvatarKey, MedLineRegimen) %>% 
  fill(AgeAtMedStart, .direction = "updown") %>% 
  ungroup() %>% 
  # Organize vars
  select(-c(MedReasonNoneGiven : MedPrimaryDiagnosisSite), 
         everything(), MedReasonNoneGiven : MedPrimaryDiagnosisSite) %>% 
  distinct()
  
# Start a clean Chemotherapy data
Medications_yes <- Medications_yes %>% 
  filter(MedicationInd == "Yes") %>% 
  mutate(regimen_line = case_when(
    is.na(MedLineRegimen)                            ~ 999,
    str_detect(MedLineRegimen, "First") &
      str_detect(MedLineRegimen, "Neoadjuvant")      ~ -1,
    str_detect(MedLineRegimen, "First") &
      str_detect(MedLineRegimen, "Adjuvant")         ~ 1,
    str_detect(MedLineRegimen, "Second")             ~ 2,
    str_detect(MedLineRegimen, "Third")              ~ 3,
    str_detect(MedLineRegimen, "Fourth")             ~ 4,
    str_detect(MedLineRegimen, "Fifth")              ~ 5,
    str_detect(MedLineRegimen, "Sixth")              ~ 6,
    str_detect(MedLineRegimen, "Seventh")            ~ 7,
    str_detect(MedLineRegimen, "Eighth")             ~ 8,
    str_detect(MedLineRegimen, "Ninth")              ~ 9,
    str_detect(MedLineRegimen, "Tenth")              ~ 10,
    str_detect(MedLineRegimen, "Eleventh")           ~ 11,
    str_detect(MedLineRegimen, "Twelfth")            ~ 12,
    MedLineRegimen == "Consolidation"                ~ 90,
    MedLineRegimen == "Maintenance"                  ~ 91,
    MedLineRegimen == "Palliative"                   ~ 92,
    TRUE                                             ~ 1000 # none?
  ), .after = MedLineRegimen) %>% 
  
  # Manually code line if thee are no line info at all
  arrange(AvatarKey, AgeAtMedStart, Medication) %>% 
  # mutate(no_line_info = case_when(
  #   is.na(MedLineRegimen)                            ~ 999
  # )) %>% 
  # group_by(AvatarKey) %>% 
  # fill(no_line_info, .direction = "updown") %>% 
  # ungroup() %>% 
  # mutate(has_some_sort_line_info = case_when(
  #   no_line_info == 999 &
  #     !is.na(MedLineRegimen)                         ~ "has_some_sort_line_info"
  # )) %>% 
  # group_by(AvatarKey) %>% 
  # fill(has_some_sort_line_info, .direction = "updown") %>% 
  # ungroup() %>% 
  group_by(AvatarKey) %>% 
  mutate(regimen_line = dense_rank(interaction(regimen_line, AgeAtMedStart)), .after = regimen_line) %>% 
  ungroup() %>% 
  # mutate(regimen_line = case_when(
  #   no_line_info == 999 &
  #     is.na(has_some_sort_line_info)                 ~ regimen_line2,
  #   TRUE                                             ~ regimen_line
  # )) %>% 
  # select(-c(no_line_info, has_some_sort_line_info, regimen_line2))
  arrange(AvatarKey, regimen_line, AgeAtMedStart, Medication)

write_rds(Medications_yes,
          paste0(here::here(),
                 "/data/processed_data",
                 "/All_PatientsMedication_LongFormat_",
                 str_remove_all(today(), "-"), ".rds"))

# For the patient with a ovary diagnosis site
Medications_yes <- Medications_yes %>% 
  left_join(., Diagnosis %>%
              select(AvatarKey, AgeAtDiagnosis, PrimaryDiagnosisSiteCode, PrimaryDiagnosisSite, number_of_dx, 
                     diagnosis_sequence, pre_cancer_info_seperated_byslash,
                     post_cancer_info_seperated_byslash),
            by = c("AvatarKey"#, "MedPrimaryDiagnosisSiteCode" = "PrimaryDiagnosisSiteCode",
                   # "MedPrimaryDiagnosisSite" = "PrimaryDiagnosisSite"
            )) |> 
  # select(AvatarKey, AgeAtMedStart, MedPrimaryDiagnosisSiteCode, MedPrimaryDiagnosisSite, 
  #        AgeAtDiagnosis, number_of_dx, Medication,
  #        diagnosis_sequence, pre_cancer_info_seperated_byslash,
  #        post_cancer_info_seperated_byslash) |> 
  mutate(filter_number_of_dx = number_of_dx, .after = number_of_dx) |> 
  group_by(AvatarKey) |> 
  fill(filter_number_of_dx, .direction = "updown") |> 
  ungroup() |> 
  mutate(isfor_ovarian_dx = case_when(
    MedPrimaryDiagnosisSite %in% c("Ovary", "Fallopian tube", "Peritoneum, NOS", "Specified parts of peritoneum")    ~ "Yes"
  )) |> 
  group_by(AvatarKey) %>% 
  mutate(drug_sequence = row_number(AvatarKey), .after = AvatarKey) %>%
  ungroup() %>% 
  unite(drug_id, c("AvatarKey", "drug_sequence"), sep = "_", remove = FALSE)



ovarian_medications <- Medications_yes |> 
  filter(!is.na(isfor_ovarian_dx) |
           (MedPrimaryDiagnosisSite == "Unknown/Not Applicable" & number_of_dx == 1) |
           (MedPrimaryDiagnosisSite == "Unknown/Not Applicable" & filter_number_of_dx == 1)) |> 
  group_by(AvatarKey) %>% 
  mutate(first_ovarian_drug_age = min(AgeAtMedStart, na.rm = TRUE),
         first_ovarian_drug_age = na_if(first_ovarian_drug_age, Inf)) %>% 
  ungroup()

not_ovarian <- Medications_yes |> 
  filter(!str_detect(drug_id, paste0(ovarian_medications$drug_id, collapse = "$|^"))) %>% 
  left_join(., ovarian_medications |> 
              select(AvatarKey, first_ovarian_drug_age, drug_sequence2 = drug_sequence) |>
              distinct(AvatarKey, .keep_all = TRUE),
            by = "AvatarKey") |> 
  # create var for pre and post dx
  mutate(pre_age = str_match(pre_cancer_info_seperated_byslash, "(\\d+\\.\\d+); ")[,2],
         pre_age = as.numeric(pre_age)) |> 
  mutate(post_age = str_match(post_cancer_info_seperated_byslash, "(\\d+\\.\\d+); |(\\d+); ")[,2],
         post_age = as.numeric(post_age)) |> 
  mutate(post_age2 = str_match(post_cancer_info_seperated_byslash, "(\\d+); ")[,2],
         post_age2 = as.numeric(post_age2)) |> 
  mutate(post_age = coalesce(post_age, post_age2)) |> select(-post_age2) |> 
  mutate(pre_or_post_note = case_when(
    diagnosis_sequence == 1 &
      AgeAtMedStart < post_age               ~ "ovarian",
    diagnosis_sequence == 1 &
      AgeAtMedStart >= post_age              ~ "post",
    
    diagnosis_sequence != 1 &
      AgeAtMedStart < AgeAtDiagnosis &
      AgeAtMedStart >= pre_age               ~ "most likely pre",
    diagnosis_sequence != 1 &
      AgeAtMedStart >= AgeAtDiagnosis &
      (AgeAtMedStart < post_age |
         is.na(post_age))                    ~ "most likely ovarian",
    diagnosis_sequence != 1 &
      AgeAtMedStart > AgeAtDiagnosis &
      AgeAtMedStart >= post_age              ~ "most likely post",
    diagnosis_sequence != 1 &
      AgeAtMedStart < AgeAtDiagnosis         ~ "most likely pre"
    # drug_sequence < drug_sequence2           ~ "pre",
    # drug_sequence > drug_sequence2           ~ "post"
  )) |> 
  mutate(pre_or_post = case_when(
    diagnosis_sequence == 1 &
      AgeAtMedStart < post_age               ~ "ovarian",
    diagnosis_sequence == 1 &
      AgeAtMedStart >= post_age              ~ "post",
    
    diagnosis_sequence != 1 &
      AgeAtMedStart < AgeAtDiagnosis &
      AgeAtMedStart >= pre_age               ~ "pre",
    diagnosis_sequence != 1 &
      AgeAtMedStart >= AgeAtDiagnosis &
      (AgeAtMedStart < post_age |
         is.na(post_age))                    ~ "ovarian",
    diagnosis_sequence != 1 &
      AgeAtMedStart > AgeAtDiagnosis &
      AgeAtMedStart >= post_age              ~ "post",
    diagnosis_sequence != 1 &
      AgeAtMedStart < AgeAtDiagnosis         ~ "pre"
  )) |> 
  # For 1 missing AgeAtMedStart
  group_by(AvatarKey, MedPrimaryDiagnosisSite) |> 
  fill(pre_or_post, pre_or_post_note, .direction = "updown") |> 
  ungroup()

more_ovarian_drug <- not_ovarian |> 
  filter(pre_or_post == "ovarian")

not_ovarian <- not_ovarian |> 
  filter(pre_or_post == "pre" | pre_or_post == "post") |> 
  select(AvatarKey, 
         AgeAtMedStart, Medication, MedPrimaryDiagnosisSite, pre_or_post) |> 
  unite(other_cancer_drug_seperated_byslash, c("AgeAtMedStart" : "MedPrimaryDiagnosisSite"), sep = "; ", remove = TRUE) |> 
  group_by(AvatarKey, pre_or_post) |> 
  summarise_at(vars(other_cancer_drug_seperated_byslash), str_c, collapse = " / ") |> 
  ungroup() |> 
  pivot_wider(id_cols = AvatarKey, 
              names_from = "pre_or_post", 
              values_from = other_cancer_drug_seperated_byslash, 
              names_glue = "{pre_or_post}_{.value}")

ovarian_medications <- ovarian_medications %>%
  bind_rows(., more_ovarian_drug) %>%
  full_join(., not_ovarian,
            by = "AvatarKey") |> 
  select(-c(drug_id, AgeAtDiagnosis, number_of_dx, filter_number_of_dx,
            diagnosis_sequence, isfor_ovarian_dx,
            pre_cancer_info_seperated_byslash,
            post_cancer_info_seperated_byslash, first_ovarian_drug_age,
            drug_sequence2, pre_age, post_age, pre_or_post_note, pre_or_post
  ))


Medications_final <- bind_rows(ovarian_medications, Medications_never) %>% 
  distinct(AvatarKey, AgeAtMedStart, Medication, AgeAtMedStop, .keep_all = TRUE) %>% 
  mutate(ever_received_chemotherapy = case_when(
    is_chemotherapy == "Yes"               ~ "Ever"
  )) %>% 
  group_by(AvatarKey, ever_received_chemotherapy) %>% 
  mutate(year_first_chemo = case_when(
    is_chemotherapy == "Yes"               ~ first(YearOfMedStart)
  ), .after = YearOfMedStart) %>% 
  
  group_by(AvatarKey) %>% 
  fill(ever_received_chemotherapy, year_first_chemo, .direction = "updown") %>% 
  mutate(year_first_medication = first(YearOfMedStart), .after = YearOfMedStart) %>% 
  group_by(AvatarKey, is_chemotherapy) %>% 
  mutate(first_chemo_age = case_when(
    is_chemotherapy == "Yes"            ~ first(AgeAtMedStart)
  ), .after = is_chemotherapy) %>% 
  ungroup() %>% 
  group_by(AvatarKey) %>% 
  fill(first_chemo_age, .direction = "updown") %>% 
  ungroup()

write_rds(Medications_final,
          paste0(here::here(),
                 "/data/processed_data",
                 "/OVCA_PatientsMedication_LongFormat_",
                 str_remove_all(today(), "-"), ".rds"))


# Create regimen - keep as "long" data
Medications_regimen <- Medications_final %>%
  # Same age start and stop
  group_by(AvatarKey, regimen_line, 
           AgeAtMedStart, MedicationInd, AgeAtMedStop, 
           ever_received_chemotherapy, first_chemo_age, year_first_chemo,
           ever_first_med_age, year_first_medication, 
           pre_other_cancer_drug_seperated_byslash, post_other_cancer_drug_seperated_byslash
  ) %>%
  summarise_at(vars(Medication, 
                    narrow_drug_class_cytotoxic_only, 
                    general_drug_class), str_c, collapse = "; ") %>%
  # same age start
  group_by(AvatarKey, regimen_line, 
           AgeAtMedStart, MedicationInd, 
           ever_received_chemotherapy, first_chemo_age, year_first_chemo,
           ever_first_med_age, year_first_medication, 
           pre_other_cancer_drug_seperated_byslash, post_other_cancer_drug_seperated_byslash
  ) %>%
  summarise_at(vars(Medication, 
                    narrow_drug_class_cytotoxic_only, 
                    general_drug_class,
                    AgeAtMedStop), str_c, collapse = "; ") %>%
  mutate(general_drug_class = paste(unique(unlist(strsplit(general_drug_class, "; "))), 
                                    collapse = "; ")) %>% 
  mutate(narrow_drug_class_cytotoxic_only = paste(unique(unlist(strsplit(narrow_drug_class_cytotoxic_only, "; "))), 
                                                  collapse = "; ")) %>% 
  # # remove regimen line - not right if NAs
  # group_by(AvatarKey, AgeAtMedStart, MedicationInd, 
  #          ever_first_med_age, year_first_medication
  # ) %>%
  # summarise_at(vars(Medication, 
  #                   AgeAtMedStop), str_c, collapse = "; ") %>%
  ungroup() %>% 
  rename(regimen_drugname = Medication) %>% 
  arrange(AvatarKey, regimen_line, AgeAtMedStart)
  
write_rds(Medications_regimen,
          paste0(here::here(),
                 "/data/processed_data",
                 "/Medication long format_",
                 str_remove_all(today(), "-"), ".rds"))
write_csv(Medications_regimen,
          paste0(here::here(),
                 "/data/processed_data",
                 "/Medication long format_",
                 str_remove_all(today(), "-"), ".csv"))

library(data.table)
Medications_wide <- dcast(setDT(Medications_regimen),
                          AvatarKey + MedicationInd + 
                            ever_received_chemotherapy + first_chemo_age + year_first_chemo +
                            ever_first_med_age + year_first_medication + 
                            pre_other_cancer_drug_seperated_byslash + post_other_cancer_drug_seperated_byslash
                          ~ rowid(AvatarKey),
                          value.var = c("AgeAtMedStart", "regimen_drugname", "AgeAtMedStop",
                                        "narrow_drug_class_cytotoxic_only", "general_drug_class")) %>% 
  rename(drug_ever = MedicationInd) %>%
  mutate(has_medication_data = "Yes")

write_rds(Medications_wide,
          paste0(here::here(),
                 "/data/processed_data",
                 "/Regimen wide format_",
                 str_remove_all(today(), "-"), ".rds"))
write_rds(Medications_wide,
          paste0(here::here(),
                 "/data/processed_data",
                 "/Regimen wide format_",
                 str_remove_all(today(), "-"), ".csv"))

# sct----
# The first exploratory analysis showed that none of 
# the patients with ovarian cancer had SCT so I remove the clean up

rm(StemCellTransplant)


# radiation----
Radiation <- Radiation %>%
  filter(str_detect(AvatarKey, paste0(ch_calls$AvatarKey, collapse = "|"))) %>% 
  mutate(across(c("AgeAtRadiationStart",
                  "AgeAtRadiationStop"), ~ case_when(
                    . == "Age 90 or older"                ~ 90,
                    . == "Unknown/Not Applicable"         ~ NA_real_,
                    TRUE                                  ~ as.numeric(.)
                  ))) %>% 
  mutate(RadDose = as.numeric(RadDose)) %>% 
  group_by(AvatarKey) %>% 
  mutate(ever_first_rad_age = min(AgeAtRadiationStart, na.rm = TRUE),
         ever_first_rad_age = na_if(ever_first_rad_age, Inf)) %>% 
  mutate(year_first_radiation = first(YearOfRadiationStart), .after = YearOfRadiationStart) %>% 
  ungroup() %>% 
  distinct(AvatarKey, AgeAtRadiationStart, RadDose, AgeAtRadiationStop, .keep_all = TRUE)

Radiation1 <- Radiation |> 
  filter(RadPrimaryDiagnosisSite %in% 
           c("Ovary", "Fallopian tube", "Peritoneum, NOS", 
             "Specified parts of peritoneum", 
             "Unknown/Not Applicable")) %>%
  # For the patient with a diagnosis site
  inner_join(., Diagnosis %>%
               select(AvatarKey, AgeAtDiagnosis, PrimaryDiagnosisSiteCode, PrimaryDiagnosisSite),
             by = c("AvatarKey"#, "RadPrimaryDiagnosisSiteCode" = "PrimaryDiagnosisSiteCode",
                    # "RadPrimaryDiagnosisSite" = "PrimaryDiagnosisSite"
                    )) #%>%
  # mutate(check = case_when(
  #   AgeAtRadiationStart < AgeAtDiagnosis     ~ "exclude"
  # ))
  # We shouldn't order by age as it would exclude unknown age of first line
  # group_by(AvatarKey) %>% 
  # mutate(first_age_is_na = first(AgeAtRadiationStart), .before = AgeAtRadiationStart) %>% 
  # ungroup()

# Radiation2 <- Radiation %>%
#   # For the patients with no diagnosis site
#   filter((!str_detect(AvatarKey, 
#                       paste0(Radiation1$AvatarKey, collapse = "|")))) %>% 
#   # There is no first line to rescue from Unknown dx for patient with only 1 dx so moving forward
#   filter(RadPrimaryDiagnosisSite == "Unknown/Not Applicable") %>%
#   left_join(., Diagnosis2 %>% select(AvatarKey, number_of_dx), by = "AvatarKey") %>% 
#   # Make sure we keep unknown only for patient that had only 1 dx
#   # So the unknown it is associated to the ovary dx
#   filter(number_of_dx == 1) %>% 
#   select(-number_of_dx)
# 
# Radiation_long <- Radiation1 %>%
#   # bind the data with a known primary site as first diagnosis
#   # then the one without 
#   bind_rows(., Radiation2) %>% 
#   full_join(., ch_calls %>% 
#               select(AvatarKey, germline_collection_age), 
#             by = "AvatarKey") %>% 
#   mutate(radiation_specimen_sequence = case_when(
#     AgeAtRadiationStart <= germline_collection_age        ~ "radiation before germline",
#     AgeAtRadiationStart > germline_collection_age         ~ "radiation after germline"
#   )) %>% 
#   select(-c(germline_collection_age))

# write_rds(Radiation_long,
#           paste0(here::here(),
#                  "/data/processed data",
#                  "/Radiation long format_",
#                  str_remove_all(today(), "-"), ".rds"))
# write_rds(Radiation_long,
#           paste0(here::here(),
#                  "/data/processed data",
#                  "/Radiation long format_",
#                  str_remove_all(today(), "-"), ".csv"))

library(data.table)
Radiation_wide <- dcast(setDT(Radiation1),
                          AvatarKey + RadiationTherapyInd + ever_first_rad_age ~ rowid(AvatarKey),
                          value.var = c("AgeAtRadiationStart", "RadDose", "AgeAtRadiationStop")) %>% 
  select(AvatarKey, radiation_ever = RadiationTherapyInd, ever_first_rad_age, 
         starts_with("AgeAtRadiationStart_"),
         starts_with("RadDose_"), starts_with("AgeAtMedStop_"), everything()) %>%
  mutate(has_radiation_data = "Yes") %>% 
  mutate(had_prior_radiation = case_when(
    AgeAtRadiationStart_1 > ever_first_rad_age       ~ "Yes"
  ))

# write_rds(Radiation_wide,
#           paste0(here::here(),
#                  "/data/processed data",
#                  "/Radiation wide format_",
#                  str_remove_all(today(), "-"), ".rds"))


rm(Radiation1, Radiation2)


# surgery ---- 
SurgeryBiopsy <- SurgeryBiopsy %>%
  filter(str_detect(AvatarKey, paste0(Diagnosis$AvatarKey, collapse = "|"))) |> 
  mutate(across(c("AgeAtSurgeryBiopsy"), ~ case_when(
                    . == "Age 90 or older"                ~ 90,
                    . == "Unknown/Not Applicable"         ~ NA_real_,
                    TRUE                                  ~ as.numeric(.)
                  ))) %>% 
  group_by(AvatarKey, SiteTherapeutic) %>% 
  mutate(ever_first_surgery_age = case_when(
    SiteTherapeutic == "Yes"                              ~ min(AgeAtSurgeryBiopsy)
  )) %>%
  ungroup() |> 
  select(AvatarKey, SiteTherapeutic, everything()) 

surgery_ever <- SurgeryBiopsy |> 
  filter(SiteTherapeutic == "Yes") |> 
  filter(PrimaryDiagnosisSite %in% 
           c("Ovary", "Fallopian tube", "Peritoneum, NOS", 
             "Specified parts of peritoneum", 
             "Unknown/Not Applicable")) |> 
  select(-c(PrimaryDiagnosisSite, PrimaryDiagnosisSiteCode)) %>%
  inner_join(., Diagnosis %>%
               select(AvatarKey, AgeAtDiagnosis),
             by = c("AvatarKey"#, "RadPrimaryDiagnosisSiteCode" = "PrimaryDiagnosisSiteCode",
                    # "RadPrimaryDiagnosisSite" = "PrimaryDiagnosisSite"
             )) |> 
  # mutate(check = case_when(
  #   AgeAtSurgeryBiopsy < AgeAtDiagnosis     ~ "exclude"
  # )) |> 
  filter(AgeAtSurgeryBiopsy >= AgeAtDiagnosis) |> 
  arrange(AvatarKey, AgeAtSurgeryBiopsy) |> 
  distinct(AvatarKey, .keep_all = TRUE) |> 
  select(-c(SurgeryBiopsyInd, AgeAtDiagnosis)) |> 
  rename(AgeAtSurgery = AgeAtSurgeryBiopsy) |> 
  mutate(surgery_ever = "Yes", .after = 1)

#   # keep only 1 record if multiple non therapeutic surg
#   mutate(keep_filter = row_number()) %>% 
#   group_by(AvatarKey) %>% 
#   mutate(has_surgery_yes = case_when(
#     SiteTherapeutic == "Yes"    ~ "Yes"
#   )) %>% 
#   fill(has_surgery_yes, .direction = "updown") %>% 
#   ungroup() %>% 
#   filter(SiteTherapeutic == "Yes" |
#            (SiteTherapeutic == "No" & keep_filter == 1 & is.na(has_surgery_yes))) %>% 
#   select(-keep_filter, -has_surgery_yes)
# 
# SurgeryBiopsy1 <- SurgeryBiopsy_ %>%
#   # For the patient with a diagnosis site
#   inner_join(., Diagnosis2 %>%
#                select(AvatarKey, PrimaryDiagnosisSiteCode, PrimaryDiagnosisSite),
#              by = c("AvatarKey", "PrimaryDiagnosisSiteCode",
#                     "PrimaryDiagnosisSite"))
# 
# SurgeryBiopsy2 <- SurgeryBiopsy_ %>%
#   # For the patients with no diagnosis site
#   filter((!str_detect(AvatarKey, 
#                       paste0(SurgeryBiopsy1$AvatarKey, collapse = "|")))) %>% 
#   left_join(., Diagnosis2 %>% select(AvatarKey, number_of_dx), by = "AvatarKey") %>% 
#   # Make sure we keep unknown only for patient that had only 1 dx
#   # So the unknown it is associated to the ovary dx
#   filter(number_of_dx == 1) %>% 
#   select(-number_of_dx)
# 
# SurgeryBiopsy_long <- SurgeryBiopsy1 %>%
#   # bind the data with a known primary site as first diagnosis
#   # then the one without 
#   bind_rows(., SurgeryBiopsy2) %>% 
#   full_join(., ch_calls %>% 
#               select(AvatarKey, germline_collection_age), 
#             by = "AvatarKey") %>% 
#   mutate(surgery_specimen_sequence = case_when(
#     AgeAtSurgeryBiopsy <= germline_collection_age        ~ "surgery before germline",
#     AgeAtSurgeryBiopsy > germline_collection_age         ~ "surgery after germline"
#   )) %>% 
#   select(-c(germline_collection_age))
# 
# # write_rds(SurgeryBiopsy_long,
# #           paste0(here::here(),
# #                  "/data/processed data",
# #                  "/Surgery long format_",
# #                  str_remove_all(today(), "-"), ".rds"))
# # write_rds(SurgeryBiopsy_long,
# #           paste0(here::here(),
# #                  "/data/processed data",
# #                  "/Surgery long format_",
# #                  str_remove_all(today(), "-"), ".csv"))
# 
# library(data.table)
# SurgeryBiopsy_wide <- dcast(setDT(SurgeryBiopsy_long),
#                         AvatarKey + SiteTherapeutic + ever_first_surgery_age ~ rowid(AvatarKey),
#                         value.var = c("AgeAtSurgeryBiopsy")) %>% 
#   select(AvatarKey, surgery_ever = SiteTherapeutic, ever_first_surgery_age, 
#          age_at_first_surgery = `1`, AgeAtSurgeryBiopsy_2 = `2`,
#          AgeAtSurgeryBiopsy_3 = `3`, AgeAtSurgeryBiopsy_4 = `4`,
#          AgeAtSurgeryBiopsy_5 = `5`, AgeAtSurgeryBiopsy_6 = `6`,
#          AgeAtSurgeryBiopsy_7 = `7`, AgeAtSurgeryBiopsy_8 = `8`,
#          AgeAtSurgeryBiopsy_9 = `9`) %>%
#   mutate(has_surgery_data = "Yes") %>% 
#   mutate(had_prior_surgery = case_when(
#     age_at_first_surgery > ever_first_surgery_age       ~ "Yes"
#   ))
# 
# rm(SurgeryBiopsy1, SurgeryBiopsy2)
surgery_never <- SurgeryBiopsy %>%
  filter((SiteTherapeutic == "No" &
            SurgeryBiopsyInd == "Yes") | 
           SurgeryBiopsyInd == "No" #|
           # !str_detect(AvatarKey, paste0(surgery_ever$AvatarKey, collapse = "|"))
  ) |> 
  select(AvatarKey) |> 
  distinct() |> 
  mutate(surgery_ever = "No")

surgery_ever <- bind_rows(surgery_ever, surgery_never) |> 
  distinct(AvatarKey, .keep_all = TRUE) |> 
  mutate(has_surgery_data = "Yes")

rm(surgery_never, SurgeryBiopsy)

# metastasis
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
  inner_join(., Diagnosis %>%
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
  left_join(., Diagnosis %>% select(AvatarKey, number_of_dx), by = "AvatarKey") %>% 
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
  left_join(., Diagnosis %>%
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
  inner_join(., ch_calls %>% select(AvatarKey, germline_collection_age), 
             by = "AvatarKey") %>% 
  # Take age at marker closest to blood collection disregarding if before or after
  mutate(int = abs(AgeAtTumorMarkerTest - germline_collection_age)) |> 
  arrange(AvatarKey, int) |> 
  distinct(AvatarKey, marker_test, .keep_all = TRUE) |> 
  select(-germline_collection_age)

# Add Flow results
TumorMarker_long <- TumorMarker_3 %>%
  mutate(marker_test = str_replace_all(str_to_lower(marker_test), " ", "_")) %>% 
  left_join(., TumorMarkerFlowPanel, 
            by = c("AvatarKey", "TumorMarkerKey")) %>% 
  mutate(FlowPanelName = str_replace_all(str_to_lower(FlowPanelName), " ", "_")) %>% 
  unite(marker_test, c(marker_test, FlowPanelName), sep = "_", na.rm = TRUE) %>% 
  mutate(marker_result = coalesce(marker_result, FlowPanelResult)) %>% 
  select(-c(FlowPanelResult, int, 
            TumorMarkerKey, AgeAtTumorMarkerTest))

# write_rds(TumorMarker_long,
#           paste0(here::here(),
#                  "/data/processed data",
#                  "/Tumor marker long format_",
#                  str_remove_all(today(), "-"), ".rds"))
# write_rds(TumorMarker_long,
#           paste0(here::here(),
#                  "/data/processed data",
#                  "/Tumor marker long format_",
#                  str_remove_all(today(), "-"), ".csv"))

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
  inner_join(., Diagnosis %>%
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


################################################################################# II ### First merge to call first treatment
# This is needed to pick the right progression / recurrence event
data <- ch_calls %>%
  full_join(., demographics, by = "AvatarKey") %>%
  inner_join(., Diagnosis, by = "AvatarKey") %>% 
  left_join(., VitalStatus, by = "AvatarKey") %>%
  left_join(., Medications_wide, by = "AvatarKey") %>%
  # left_join(., StemCellTransplant, by = "AvatarKey") %>%
  left_join(., Radiation_wide, by = "AvatarKey") %>%
  left_join(., surgery_ever, by = "AvatarKey") %>%
  left_join(., PatientHistory %>%
              select(AvatarKey, AgeAtLastUpdate, SmokingStatus, AlcoholUse, PreExistHPV),
            by = "AvatarKey") %>%
  left_join(., PhysicalAssessment, by = "AvatarKey") %>%
  left_join(., MetastaticDisease, by = "AvatarKey") %>%
  left_join(., last_date, by = "AvatarKey") %>%
  left_join(., FamilyHistory, by = "AvatarKey") %>% 
  left_join(., CytogeneticAbnormalities, by = "AvatarKey") %>%
  left_join(., TumorMarker_wide, by = "AvatarKey")

# rm(demographics, VitalStatus, 
#    CytogeneticAbnormalities, TumorMarker,
#    PhysicalAssessment, Medications,
#    Radiation,
#    MetastaticDisease,
#    last_date, FamilyHistory,
#    ch_calls)

  
treatment <- data |> 
  mutate(treatment_sequence = case_when(
    is.na(has_medication_data)                             ~ NA_character_,
    is.na(has_surgery_data)                                ~ NA_character_,
    AgeAtMedStart_1 <= AgeAtSurgery                  ~ "drug/surgery",
    AgeAtMedStart_1 > AgeAtSurgery                   ~ "surgery/drug",
    drug_ever == "No" &
      surgery_ever == "No"                                 ~ "No drug or surgery",
    drug_ever == "Yes" &
      surgery_ever == "No"                                 ~ "drug only",
    surgery_ever == "Yes" &
      drug_ever == "No"                                    ~ "surgery only",
    is.na(AgeAtMedStart_1) &
      is.na(AgeAtSurgery)                            ~ "missing drug and surgery age",
    is.na(AgeAtMedStart_1)                                 ~ "missing drug age",
    is.na(AgeAtSurgery)                              ~ "missing surgery age",
    TRUE      ~ paste(AgeAtMedStart_1, AgeAtSurgery, sep = "; ")
  )) |> 
  mutate(first_treatment = case_when(
    str_detect(treatment_sequence, "^surg")                ~ "Surgery",
    str_detect(treatment_sequence, "^drug")                ~ "Drugs",
  )) |> 
  mutate(age_at_first_treatment = case_when(
    first_treatment == "Surgery"                           ~ AgeAtSurgery,
    first_treatment == "Drugs"                             ~ AgeAtMedStart_1
  )) |> 
  mutate(upfront_treatment = case_when(
    treatment_sequence == "surgery/drug"                   ~ "Upfront surgery",
    treatment_sequence == "drug/surgery"                   ~ "Upfront systemic"
  ))


################################################################################# II ### Data cleaning part 2
# outcomes PFS----
Outcomes_save <- Outcomes
Outcomes <- Outcomes_save

Outcomes <- Outcomes %>%
  filter(str_detect(AvatarKey, paste0(treatment$AvatarKey, collapse = "|"))) %>% 
  mutate(across(everything(), ~ 
                  str_replace_all(., 
                                  "Unknown/Not Applicable|Unknown/Not Performed|Age Unknown/Not Recorded",
                                  NA_character_))) |> 
  distinct() |> 
  purrr::keep(~!all(is.na(.))) |> 
  mutate(not_real_pfsage = case_when(
    AgeAtProgRecur == "Age 90 or older"   ~ "Age 90 or older"
  )) %>%
  mutate(across(c("AgeAtProgRecur",
                  "AgeAtCurrentDiseaseStatus"), ~ case_when(
                    . == "Age 90 or older"                ~ 90,
                    . == "Unknown/Not Applicable"         ~ NA_real_,
                    TRUE                                  ~ as.numeric(.)
                  ))) %>%
  select(AvatarKey, ProgRecurInd, 
         AgeAtProgRecur, 
         YearOfProgRecur, AgeAtCurrentDiseaseStatus, AgeAtPerformStatusMostRecent,
         OutcomesPrimaryDiagnosisSite, OutcomesPrimaryDiagnosisSiteCode, 
         everything())

Outcomes1 <- Outcomes |> 
  # Select the outcome associated to the correct dx
  mutate(is_ovarian_pfs_info = case_when(
    OutcomesPrimaryDiagnosisSite %in% c("Ovary", "Fallopian tube", "Peritoneum, NOS", "Specified parts of peritoneum")    ~ "Yes"
  )) %>%
  # Combine age at event + age at no event for arranging
  mutate(age_at_disease_check = coalesce(AgeAtProgRecur, AgeAtCurrentDiseaseStatus), 
         .after = AgeAtCurrentDiseaseStatus) %>% 
  arrange(AvatarKey, age_at_disease_check) %>% 
  # Select the correct outcome date data - after dx and before over dx?
  inner_join(., treatment %>%
               select(AvatarKey, AgeAtDiagnosis, number_of_dx, 
                      age_at_first_treatment#,
                      # post_cancer_info_seperated_byslash
               ),
             by = "AvatarKey") |> 
  # Need to be after dx
  filter(age_at_disease_check > AgeAtDiagnosis) |> 
  # Before post dx? Do not run
  # Some event after post cancer but all of them are noted for ovarian cancer - keep in
  # mutate(post_age = str_match(post_cancer_info_seperated_byslash, "(\\d+\\.\\d+); |(\\d+); ")[,2],
  #        post_age = as.numeric(post_age)) |> 
  # mutate(post_age2 = str_match(post_cancer_info_seperated_byslash, "(\\d+); ")[,2],
  #        post_age2 = as.numeric(post_age2)) |> 
  # mutate(post_age = coalesce(post_age, post_age2)) |> select(-post_age2) |> 
  
  # Check outcomes info for missing primary site 
  # group_by(AvatarKey) |>
  # fill(is_ovarian_pfs_info, .direction = "updown") |>
  # ungroup() |> 
  
  # Select the outcome associated to the correct dx
  filter(is_ovarian_pfs_info == "Yes" | 
           # 3 patients have outcomes info for missing site but they all have only 1 dx 
           # keep info into data and age meet dx age
           (is.na(is_ovarian_pfs_info) & number_of_dx == 1) |
           # Patients with outcomes data for another site but none for ovarian
           # Those outcomes are truly for the other site - exclude - none
           (is.na(is_ovarian_pfs_info) & is.na(OutcomesPrimaryDiagnosisSite))
         ) |> 
  select(-c())

Outcomes2 <- Outcomes1 |> 
  # Need to be after first treatment
  filter(age_at_disease_check > age_at_first_treatment)

Outcomes_yes <- Outcomes2 %>%
  # Select Yes event
  filter(!is.na(ProgRecurInd)) |> 
  mutate(AgeAtProgRecur = coalesce(AgeAtProgRecur, AgeAtCurrentDiseaseStatus)) |> 
  arrange(AvatarKey, AgeAtProgRecur) %>% 
  # select first event
  distinct(AvatarKey, .keep_all = TRUE)

pfs_free <- Outcomes2 %>%
  # Select records stating No event happened
  filter(is.na(ProgRecurInd)) |> 
  # Remove patients who had an event
  filter(!str_detect(AvatarKey, paste0(Outcomes_yes$AvatarKey, collapse = "|"))) |> 
  # For no event - we want to pick the latest date of disease check 
  arrange(AvatarKey, desc(age_at_disease_check)) |> 
  distinct(AvatarKey, .keep_all = TRUE) |> 
  mutate(ProgRecurInd = "No")

Outcomes <- bind_rows(Outcomes_yes, pfs_free) |> 
  mutate(has_outcomes_data = "Yes") |> 
  distinct(AvatarKey, .keep_all = TRUE) |> 
  select(-c(AgeAtDiagnosis, number_of_dx,
            age_at_first_treatment,
            age_at_first_treatment))


################################################################################# III ### Merging final data
ovarian_data <- treatment %>%
  left_join(., Outcomes, by = "AvatarKey") 


# Save
write_csv(ovarian_data,
          paste0(here::here(), 
                 "/data/processed_data",
                 "/CHinOvary_CleanedClinicalData_",
                 str_remove_all(today(), "-"), ".csv"))
write_rds(ovarian_data, 
          paste0(here::here(), 
                 "/data/processed_data",
                 "/CHinOvary_CleanedClinicalData_",
                 str_remove_all(today(), "-"), ".rds"))
write_csv(ovarian_data,
          paste0(path_gillis,
                 "/ProcessedData",
                 "/CHinOvary_CleanedClinicalData_",
                 str_remove_all(today(), "-"), ".csv"))

path_peres <- fs::path("", "Volumes", "Peres_Research",
                      "ORIEN analysis", "Ovary CH")

write_csv(ovarian_data, 
          paste0(path_peres, 
                 "/data/processed data",
                 "/CHinOvary_CleanedClinicalData_",
                 str_remove_all(today(), "-"), ".csv"))
# End data cleaning

