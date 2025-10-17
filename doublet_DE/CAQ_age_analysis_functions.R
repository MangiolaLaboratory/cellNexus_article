library(Seurat)
library(purrr)
library(magrittr)
library(cellNexus)
library(dplyr)
library(readr)
library(forcats)
library(glue)
# devtools::load_all("~/labHead/cellNexus/")

#--------------------#
# QUERY FOR DHARMESH
#--------------------#

age_bin <- function(age_days, sex) {
  # Convert age in days to age in years
  age_years <- age_days / 365.25

  # Initialise an empty vector to store the results
  age_bins <- vector("character", length(age_years))

  # Define average thresholds for "unknown" sex based on midpoint between male and female stages
  unknown_thresholds <- c(3, 13, 20, 38, 52)

  # Loop through each element to assign the appropriate bin based on sex and age
  for (i in seq_along(age_years)) {
    if (sex[i] == "male") {
      age_bins[i] <- dplyr::case_when(
        age_years[i] < 3 ~ "Infancy",
        age_years[i] < 13 ~ "Childhood",
        age_years[i] < 21 ~ "Adolescence",
        age_years[i] < 40 ~ "Young Adulthood",
        age_years[i] < 55 ~ "Middle Age",
        age_years[i] >= 55 ~ "Senior",
        TRUE ~ NA_character_
      )
    } else if (sex[i] == "female") {
      age_bins[i] <- dplyr::case_when(
        age_years[i] < 3 ~ "Infancy",
        age_years[i] < 13 ~ "Childhood",
        age_years[i] < 19 ~ "Adolescence",
        age_years[i] < 36 ~ "Young Adulthood",
        age_years[i] < 50 ~ "Middle Age",
        age_years[i] >= 50 ~ "Senior",
        TRUE ~ NA_character_
      )
    } else if (sex[i] == "unknown") {
      age_bins[i] <- dplyr::case_when(
        age_years[i] < unknown_thresholds[1] ~ "Infancy",
        age_years[i] < unknown_thresholds[2] ~ "Childhood",
        age_years[i] < unknown_thresholds[3] ~ "Adolescence",
        age_years[i] < unknown_thresholds[4] ~ "Young Adulthood",
        age_years[i] < unknown_thresholds[5] ~ "Middle Age",
        age_years[i] >= unknown_thresholds[5] ~ "Senior",
        TRUE ~ NA_character_
      )
    } else {
      stop("Each element of 'sex' must be either 'male', 'female', or 'unknown'.")
    }
  }

  return(age_bins)
}


coarse_age_bin <- function(age_days, sex) {
  # Convert age in days to age in years
  age_years <- age_days / 365.25
  
  # Initialise an empty vector to store the results
  age_bins <- vector("character", length(age_years))
  
  # Define average thresholds for "unknown" sex based on midpoint between male and female stages
  unknown_thresholds <- c(18, 60)
  
  # Loop through each element to assign the appropriate bin based on sex and age
  for (i in seq_along(age_years)) {
    if (sex[i] == "male") {
      age_bins[i] <- dplyr::case_when(
        age_years[i] < 18 ~ "0-17",
        age_years[i] < 60 ~ "18-59",
        age_years[i] >= 60 ~ "60+",
        TRUE ~ NA_character_
      )
    } else if (sex[i] == "female") {
      age_bins[i] <- dplyr::case_when(
        age_years[i] < 18 ~ "0-17",
        age_years[i] < 60 ~ "18-59",
        age_years[i] >= 60 ~ "60+",
        TRUE ~ NA_character_
      )
    } else if (sex[i] == "unknown") {
      age_bins[i] <- dplyr::case_when(
        age_years[i] < unknown_thresholds[1] ~ "0-17",
        age_years[i] < unknown_thresholds[2] ~ "18-59",
        age_years[i] >= unknown_thresholds[2] ~ "60+",
        TRUE ~ NA_character_
      )
    } else {
      stop("Each element of 'sex' must be either 'male', 'female', or 'unknown'.")
    }
  }
  
  return(age_bins)
}

edit_covariates_from_dharmesh = function(tbl, disease_tbl) {
  ethnicity_grouped <- tribble(
    ~self_reported_ethnicity, ~ethnicity_groups,
    "unknown", "Other/Unknown",
    "European", "European",
    "Korean", "East Asian",
    "Asian", "East Asian",
    "Japanese", "East Asian",
    "African American", "African",
    "Hispanic or Latin American", "Hispanic/Latin American",
    "Singaporean Chinese", "East Asian",
    "Han Chinese", "East Asian",
    "Singaporean Indian", "South Asian",
    "Singaporean Malay", "Other/Unknown",
    "British", "European",
    "African", "African",
    "South Asian", "South Asian",
    "European American", "European",
    "East Asian", "East Asian",
    "American", "Other/Unknown",
    "African American or Afro-Caribbean", "African",
    "Oceanian", "Native American & Pacific Islander",
    "Jewish Israeli", "Middle Eastern & North African",
    "Chinese", "East Asian",
    "South East Asian", "Other/Unknown",
    "Greater Middle Eastern  (Middle Eastern or North African or Persian)", "Middle Eastern & North African",
    "Native American", "Native American & Pacific Islander",
    "Pacific Islander", "Native American & Pacific Islander",
    "Finnish", "European",
    "Bangladeshi", "South Asian",
    "Native American,Hispanic or Latin American", "Hispanic/Latin American",
    "Irish", "European",
    "Iraqi", "Middle Eastern & North African",
    "European,Asian", "European"
  )

  assay_data_grouped <- tribble(
    ~assay, ~assay_groups,
    "10x 3' v2", "10x Genomics 3",
    "10x 3' v3", "10x Genomics 3",
    "10x 5' v2", "10x Genomics 5",
    "10x 5' v1", "10x Genomics 5",
    "MARS-seq", "Plate based Technologies",
    "10x 3' transcription profiling", "10x Genomics 3",
    "10x 5' transcription profiling", "10x Genomics 5",
    "Smart-seq2", "Smart seq",
    "microwell-seq", "Microwell Technologies",
    "TruDrop", "TruDrop",
    "Drop-seq", "Drop based Technologies",
    "Seq-Well S3", "Microwell Technologies",
    "GEXSCOPE technology", "Other Technologies",
    "Seq-Well", "Microwell Technologies",
    "sci-RNA-seq", "Other Technologies",
    "10x 3' v1", "10x Genomics 3",
    "BD Rhapsody Whole Transcriptome Analysis", "Other Technologies",
    "BD Rhapsody Targeted mRNA", "Other Technologies",
    "CEL-seq2", "Plate based Technologies",
    "SPLiT-seq", "Other Technologies",
    "STRT-seq", "Plate based Technologies",
    "inDrop", "Drop based Technologies",
    "Smart-seq v4", "Smart seq",
    "ScaleBio single cell RNA sequencing", "Other Technologies"
  )


  disease_data_grouped <- tribble(
    ~disease, ~disease_groups,

    # Normal control
    "normal", "Normal",

    # Isolated Diseases
    "COVID-19", "COVID-19 related",
    "post-COVID-19 disorder", "COVID-19 related",
    "long COVID-19", "COVID-19 related",
    "glioblastoma", "Glioblastoma",
    "lung adenocarcinoma", "Lung Adenocarcinoma",
    "systemic lupus erythematosus", "Systemic Lupus Erythematosus",

    # Infectious and Immune-related Diseases (other than COVID-19)
    "Crohn disease", "Infectious and Immune-related Diseases",
    "Crohn ileitis", "Infectious and Immune-related Diseases",
    "pneumonia", "Infectious and Immune-related Diseases",
    "common variable immunodeficiency", "Infectious and Immune-related Diseases",
    "toxoplasmosis", "Infectious and Immune-related Diseases",
    "Plasmodium malariae malaria", "Infectious and Immune-related Diseases",
    "type 1 diabetes mellitus", "Infectious and Immune-related Diseases",
    "influenza", "Infectious and Immune-related Diseases",
    "chronic rhinitis", "Infectious and Immune-related Diseases",
    "periodontitis", "Infectious and Immune-related Diseases",
    "localized scleroderma", "Infectious and Immune-related Diseases",
    "lymphangioleiomyomatosis", "Infectious and Immune-related Diseases",
    "listeriosis", "Infectious and Immune-related Diseases",

    # Cancer (other than isolated cancers)
    "squamous cell lung carcinoma", "Cancer",
    "small cell lung carcinoma", "Cancer",
    "non-small cell lung carcinoma", "Cancer",
    "breast carcinoma", "Cancer",
    "breast cancer", "Cancer",
    "luminal B breast carcinoma", "Cancer",
    "luminal A breast carcinoma", "Cancer",
    "triple-negative breast carcinoma", "Cancer",
    "gastric cancer", "Cancer",
    "colorectal cancer", "Cancer",
    "colon sessile serrated adenoma/polyp", "Cancer",
    "follicular lymphoma", "Cancer",
    "B-cell acute lymphoblastic leukemia", "Cancer",
    "B-cell non-Hodgkin lymphoma", "Cancer",
    "acute myeloid leukemia", "Cancer",
    "acute promyelocytic leukemia", "Cancer",
    "plasma cell myeloma", "Cancer",
    "clear cell renal carcinoma", "Cancer",
    "nonpapillary renal cell carcinoma", "Cancer",
    "basal cell carcinoma", "Cancer",
    "colorectal neoplasm", "Cancer",
    "adenocarcinoma", "Cancer",
    "chromophobe renal cell carcinoma", "Cancer",
    "neuroendocrine carcinoma", "Cancer",
    "lung large cell carcinoma", "Cancer",
    "tongue cancer", "Cancer",
    "Wilms tumor", "Cancer",
    "pleomorphic carcinoma", "Cancer",
    "blastoma", "Cancer",

    # Neurodegenerative and Neurological Disorders
    "dementia", "Neurodegenerative and Neurological Disorders",
    "Alzheimer disease", "Neurodegenerative and Neurological Disorders",
    "Parkinson disease", "Neurodegenerative and Neurological Disorders",
    "amyotrophic lateral sclerosis", "Neurodegenerative and Neurological Disorders",
    "multiple sclerosis", "Neurodegenerative and Neurological Disorders",
    "Down syndrome", "Neurodegenerative and Neurological Disorders",
    "trisomy 18", "Neurodegenerative and Neurological Disorders",
    "frontotemporal dementia", "Neurodegenerative and Neurological Disorders",
    "temporal lobe epilepsy", "Neurodegenerative and Neurological Disorders",
    "Lewy body dementia", "Neurodegenerative and Neurological Disorders",
    "amyotrophic lateral sclerosis 26 with or without frontotemporal dementia", "Neurodegenerative and Neurological Disorders",

    # Respiratory Conditions
    "pulmonary fibrosis", "Respiratory Conditions",
    "respiratory system disorder", "Respiratory Conditions",
    "chronic obstructive pulmonary disease", "Respiratory Conditions",
    "cystic fibrosis", "Respiratory Conditions",
    "interstitial lung disease", "Respiratory Conditions",
    "hypersensitivity pneumonitis", "Respiratory Conditions",
    "non-specific interstitial pneumonia", "Respiratory Conditions",
    "aspiration pneumonia", "Respiratory Conditions",
    "pulmonary emphysema", "Respiratory Conditions",
    "pulmonary sarcoidosis", "Respiratory Conditions",

    # Cardiovascular Diseases
    "myocardial infarction", "Cardiovascular Diseases",
    "acute myocardial infarction", "Cardiovascular Diseases",
    "dilated cardiomyopathy", "Cardiovascular Diseases",
    "heart failure", "Cardiovascular Diseases",
    "arrhythmogenic right ventricular cardiomyopathy", "Cardiovascular Diseases",
    "congenital heart disease", "Cardiovascular Diseases",
    "non-compaction cardiomyopathy", "Cardiovascular Diseases",
    "cardiomyopathy", "Cardiovascular Diseases",
    "heart disorder", "Cardiovascular Diseases",

    # Metabolic and Other Disorders
    "type 2 diabetes mellitus", "Metabolic and Other Disorders",
    "chronic kidney disease", "Metabolic and Other Disorders",
    "digestive system disorder", "Metabolic and Other Disorders",
    "primary sclerosing cholangitis", "Metabolic and Other Disorders",
    "gastritis", "Metabolic and Other Disorders",
    "acute kidney failure", "Metabolic and Other Disorders",
    "tubular adenoma", "Metabolic and Other Disorders",
    "benign prostatic hyperplasia", "Metabolic and Other Disorders",
    "opiate dependence", "Metabolic and Other Disorders",
    "gingivitis", "Metabolic and Other Disorders",
    "hyperplastic polyp", "Metabolic and Other Disorders",
    "clonal hematopoiesis", "Metabolic and Other Disorders",
    "epilepsy", "Metabolic and Other Disorders",
    "age related macular degeneration 7", "Metabolic and Other Disorders",
    "kidney benign neoplasm", "Metabolic and Other Disorders",
    "malignant pancreatic neoplasm", "Metabolic and Other Disorders",
    "cataract", "Metabolic and Other Disorders",
    "macular degeneration", "Metabolic and Other Disorders",
    "hydrosalpinx", "Metabolic and Other Disorders",
    "tubulovillous adenoma", "Metabolic and Other Disorders",
    "gastric intestinal metaplasia", "Metabolic and Other Disorders",
    "Barrett esophagus", "Metabolic and Other Disorders",

    # Other Diseases
    "injury", "Other Diseases",
    "anencephaly", "Other Diseases",
    "primary biliary cholangitis", "Other Diseases",
    "keloid", "Other Diseases",
    "kidney oncocytoma", "Other Diseases",
    "respiratory failure", "Other Diseases",
    "pilocytic astrocytoma", "Other Diseases"
  )

  disease_data_grouped =
    disease_data_grouped |>
    select(-disease_groups) |>
    left_join(disease_tbl) |>
    mutate(disease_groups = if_else(disease_groups |> is.na(), "other", disease_groups))

  age_bin_table =
    tbl |>
    distinct(age_days, sex) |>
    filter(!age_days |> is.na()) |>
    mutate(sex = if_else(sex |> is.na(), "unknown", sex)) |>
    as_tibble() |>
    mutate(age_bin = age_bin(age_days, sex))

  tbl |>
    # TECH
    left_join(assay_data_grouped, copy = TRUE) |>
    # DISEASE
    left_join(disease_data_grouped, copy = TRUE) |>
    # TEMPORARY. de-group pancreas and liver
    mutate(tissue_groups = case_when(
      tissue %in% c("gallbladder") ~ "gallbladder",
      tissue %in% c("pancreas", "exocrine pancreas") ~ "pancreas",
      tissue %in% c("liver", "caudate lobe of liver", "hepatic cecum") ~ "liver",
      TRUE ~ tissue_groups
    )) |>
    # SEX edit
    mutate(sex = if_else(sex |> is.na(), "unknown", sex)) |>
    # Age
    filter(age_days > 365) |>
    left_join(age_bin_table, copy = TRUE) |>
    # ETHNICITY
    left_join(ethnicity_grouped, copy = TRUE) |>
    dplyr::select(cell_id, sample_id, donor_id, dataset_id, file_id_cellNexus_single_cell, title, collection_id, age_days, age_bin, sex, ethnicity_groups, tissue_groups, tissue, assay_groups, cell_type_unified_ensemble, cell_type, disease_groups, atlas_id) |>
    as_tibble() |>
    # Set intercept
    mutate(
      ethnicity_groups = fct_relevel(ethnicity_groups, "European"),
      assay_groups = fct_relevel(assay_groups, "10x Genomics 3"),
      disease_groups = fct_relevel(disease_groups, "Normal"),
      age_bin = fct_relevel(age_bin, "Adolescence")
    ) |>
    # Center based on adolescence
    mutate(age_days_scaled = age_days |> scale(center = 15 * 365) |> as.numeric())
}

edit_covariates_from_stefano = function(tbl){
  
  tissue_grouped = list(
    
    # Respiratory System
    "respiratory system" = c(
      "lung", "lung parenchyma", "alveolus of lung",  "bronchus",
      "respiratory airway", "pleura", "pleural effusion", "middle lobe of right lung",
      "upper lobe of left lung", "lower lobe of left lung", "upper lobe of right lung",
      "lower lobe of right lung", "lingula of left lung", "right lung", "left lung"
    ),
    
    "trachea" = c( "epithelium of trachea", "trachea"),
    
    # Cardiovascular System
    "cardiovascular system" = c(
      "heart", "heart left ventricle", "heart right ventricle", "cardiac ventricle",
      "cardiac atrium", "right cardiac atrium", "left cardiac atrium", "apex of heart",
      "aorta", "coronary artery", 
      "venous blood", "anterior wall of left ventricle", "myocardium", "interventricular septum", "ventricular tissue", "basal zone of heart"
    ),
    
    "vasculature" = c("kidney blood vessel", "artery", "vein", "vasculature", "mesenteric artery"),
    
    # Umbilical Cord Blood
    "umbilical cord blood" = "umbilical cord blood",
    
    # Oesophagus
    "oesophagus" = c(
      "esophagus", "lower esophagus", "esophagus muscularis mucosa",
      "submucosal esophageal gland",
      
      # Epithelium
      "epithelium of esophagus"
    ),
    
    # Stomach
    "stomach" = c(
      "stomach", "body of stomach", "cardia of stomach"
    ),
    
    # Small Intestine
    "small intestine" = c(
      "small intestine", "duodenum", "jejunum", "ileum",
      
      # Epithelium
      "epithelium of small intestine", "jejunal epithelium", "ileal epithelium",
      "submucosa of ileum", "lamina propria of small intestine"
    ),
    
    # Large Intestine
    "large intestine" = c(
      "large intestine", "colon", "left colon", "right colon",
      "sigmoid colon", "descending colon", "transverse colon",
      "ascending colon", "hepatic flexure of colon", "caecum",
      "rectum", "appendix", "vermiform appendix",
      
      # epithelium
      "colonic epithelium", "submucosa of ascending colon", "lamina propria of large intestine",
      "mucosa of colon", "lamina propria of mucosa of colon", "caecum epithelium"
    ),
    
    # Digestive System (General)
    "digestive system (general)" = c(
      "intestine", "hindgut", "lamina propria", "mucosa"
    ),
    
    # Nasal, Oral, and Pharyngeal Regions
    "nasal, oral, and pharyngeal regions" = c(
      "nasal cavity", "nasopharynx", "oral mucosa", "tongue", "anterior part of tongue",
      "posterior part of tongue", "gingiva", "nose", "saliva"
    ),
    
    # Cerebral Lobes and Cortical Areas
    "cerebral lobes and cortical areas" = c(
      "frontal lobe", "left frontal lobe", "right frontal lobe", "primary motor cortex",
      "dorsolateral prefrontal cortex", "superior frontal gyrus", "orbitofrontal cortex",
      "medial orbital frontal cortex", "Broca's area", "prefrontal cortex",
      "temporal lobe", "left temporal lobe", "right temporal lobe", 
      "angular gyrus", "entorhinal cortex",
      "parietal lobe", "left parietal lobe", "right parietal lobe", "primary somatosensory cortex",
      "occipital lobe", "right occipital lobe", "primary visual cortex",
      "occipital cortex", "insular cortex", "parietal cortex", "temporal cortex",
      "frontal cortex", "Brodmann (1909) area 4", "temporoparietal junction",
      "middle temporal gyrus", "cingulate cortex", "brain", "brain white matter", "cerebral cortex", "cerebral nuclei"
    ),
    
    # Limbic and Basal Systems
    "limbic and basal systems" = c(
      "anterior cingulate cortex", "anterior cingulate gyrus", "hippocampal formation",
      "hypothalamus", "thalamic complex", "dentate nucleus", "basal ganglion",
      "caudate nucleus", "putamen", "substantia nigra pars compacta",
      "lateral ganglionic eminence", "medial ganglionic eminence",
      "caudal ganglionic eminence", "ganglionic eminence"
    ),
    
    # Brainstem and Cerebellar Structures
    "brainstem and cerebellar structures" = c(
      "pons", "midbrain", "myelencephalon", "telencephalon", "forebrain",
      "cerebellum", "cerebellum vermis lobule", "cerebellar cortex",
      "hemisphere part of cerebellar posterior lobe", "white matter of cerebellum"
    ),
    
    # General Brain and Major Structures
    "general brain and major structures" = c(
      "spinal cord", "neural tube", "cervical spinal cord white matter"
    ),
    
    # Muscular System (Skeletal Muscles)
    "muscular system (skeletal muscles)" = c(
      "rectus abdominis muscle", "gastrocnemius", "muscle of abdomen", "muscle organ",
      "muscle tissue", "pelvic diaphragm muscle", "skeletal muscle tissue", "muscle of pelvic diaphragm"
    ),
    
    # Connective Tissue
    "connective tissue" = c(
      "connective tissue", "tendon of semitendinosus", "vault of skull", "bone spine",
      "rib"
    ),
    
    # Adipose Tissue
    "adipose tissue" = c(
      "adipose tissue", "subcutaneous adipose tissue", "visceral abdominal adipose tissue",
      "perirenal fat", "omental fat pad", "subcutaneous abdominal adipose tissue",
      "abdominal adipose tissue"
    ),
    
    # Endocrine System
    "endocrine system" = c(
      "thyroid gland", "adrenal tissue", "adrenal gland", "islet of Langerhans",
      "endocrine pancreas", "pineal gland"
    ),
    
    # Lymphatic System
    "lymphatic system" = c(
      "lymph node", "mesenteric lymph node", "thoracic lymph node",
      "cervical lymph node", "bronchopulmonary lymph node", "tonsil", "inguinal lymph node"
    ),
    
    # Integumentary System (Skin)
    "integumentary system (skin)" = c(
      "skin of abdomen", "skin of forearm", "skin of scalp", "skin of face", "skin of leg",
      "skin of chest", "skin of back", "skin of hip", "skin of body", "skin of cheek",
      "skin of temple", "skin of shoulder", "skin of external ear", "skin of trunk",
      "skin of prepuce of penis", "skin epidermis", "arm skin", "lower leg skin",
      "hindlimb skin", "zone of skin", "dermis", "skin of nose", "skin of forehead",
      "skin of pes", "axilla"
    ),
    
    # Gastrointestinal Accessory Organs
    "gallbladder" =  "gallbladder",
    
    # Gastrointestinal Accessory Organs
    "pancreas" = c( "pancreas", "exocrine pancreas" ),
    
    # Gastrointestinal Accessory Organs
    "liver" = c( "liver", "caudate lobe of liver", "hepatic cecum" ),
    
    # Spleen
    "spleen" = "spleen",
    
    # Thymus
    "thymus" = "thymus",
    
    # Blood
    "blood" = "blood",
    
    # Bone Marrow
    "bone marrow" = "bone marrow",
    
    # Female Reproductive System
    "female reproductive system" = c(
      "uterus", "myometrium", "fallopian tube", "ampulla of uterine tube",
      "fimbria of uterine tube", "uterine cervix", "endometrium",
      "decidua", "decidua basalis", "placenta", "yolk sac", "isthmus of fallopian tube"
    ),
    "ovary" = "ovary", 
    
    # Male Reproductive System
    "male reproductive system (other)" = c(
      "testis", "gonad"
    ),
    
    # Prostate
    "prostate" = c(
      "prostate gland", "transition zone of prostate", "peripheral zone of prostate"
    ),
    
    # Renal System
    "renal system" = c(
      "kidney", "cortex of kidney", "renal medulla", "renal papilla",
      "renal pelvis", "ureter", "bladder organ"
    ),
    
    # Miscellaneous Glands
    "miscellaneous glands" = c(
      "parotid gland", "lacrimal gland", "sublingual gland", "mammary gland",
      "chorionic villus"
    ),
    
    # Eye and Visual-Related Structures
    "sensory-related structures" = c(
      "retina",
      "retinal neural layer",
      "macula lutea",
      "macula lutea proper",
      "sclera",
      "trabecular meshwork",
      "conjunctiva",
      "pigment epithelium of eye",
      "cornea",
      "iris",
      "ciliary body",
      "peripheral region of retina",
      "eye trabecular meshwork",
      "perifoveal part of retina",
      "choroid plexus",
      "lens of camera-type eye",
      "corneo-scleral junction",
      "fovea centralis",
      "eye",
      "inner ear",
      "vestibular system",
      "primary auditory cortex"
    ),
    
    # Digestive Tract Junctions and Connections
    "digestive tract junctions and connections" = c(
      "esophagogastric junction", "duodeno-jejunal junction", "hepatopancreatic ampulla",
      "hepatopancreatic duct", "pyloric antrum"
    ),
    
    # Peritoneal and Abdominal Cavity Structures
    "peritoneal and abdominal cavity structures" = c(
      "peritoneum", "omentum", "retroperitoneum", "mesentery"
    ),
    
    # Breast
    "breast" = c(
      "breast", "upper outer quadrant of breast"
    )
  ) |> 
    enframe(name ="tissue_groups") |> 
    distinct() |> 
    unnest(value) |>
    rename(tissue = value) |>
    mutate()
  
  ethnicity_grouped <- tribble(
    ~self_reported_ethnicity, ~ethnicity_groups,
    "unknown", "Other/Unknown",
    "European", "European",
    "Korean", "East Asian",
    "Asian", "East Asian",
    "Japanese", "Japanese",
    "African American", "African",
    "Hispanic or Latin American", "Hispanic/Latin American",
    "Singaporean Chinese", "East Asian",
    "Han Chinese", "East Asian",
    "Singaporean Indian", "South Asian",
    "Singaporean Malay", "Other/Unknown",
    "British", "European",
    "African", "African",
    "South Asian", "South Asian",
    "European American", "European",
    "East Asian", "East Asian",
    "American", "Other/Unknown",
    "African American or Afro-Caribbean", "African",
    "Oceanian", "Native American & Pacific Islander",
    "Jewish Israeli", "Middle Eastern & North African",
    "Chinese", "East Asian",
    "South East Asian", "Other/Unknown",
    "Greater Middle Eastern  (Middle Eastern or North African or Persian)", "Middle Eastern & North African",
    "Native American", "Native American & Pacific Islander",
    "Pacific Islander", "Native American & Pacific Islander",
    "Finnish", "European",
    "Bangladeshi", "South Asian",
    "Native American,Hispanic or Latin American", "Hispanic/Latin American",
    "Irish", "European",
    "Iraqi", "Middle Eastern & North African",
    "European,Asian", "European"
  )
  
  assay_data_grouped <- tribble(
    ~assay, ~assay_groups,
    "10x 3' v2", "10x Genomics 3",
    "10x 3' v3", "10x Genomics 3",
    "10x 5' v2", "10x Genomics 5",
    "10x 5' v1", "10x Genomics 5",
    "MARS-seq", "Plate based Technologies",
    "10x 3' transcription profiling", "10x Genomics 3",
    "10x 5' transcription profiling", "10x Genomics 5",
    "Smart-seq2", "Smart seq",
    "microwell-seq", "Microwell Technologies",
    "TruDrop", "TruDrop",
    "Drop-seq", "Drop based Technologies",
    "Seq-Well S3", "Microwell Technologies",
    "GEXSCOPE technology", "Other Technologies",
    "Seq-Well", "Microwell Technologies",
    "sci-RNA-seq", "Other Technologies",
    "10x 3' v1", "10x Genomics 3",
    "BD Rhapsody Whole Transcriptome Analysis", "Other Technologies",
    "BD Rhapsody Targeted mRNA", "Other Technologies",
    "CEL-seq2", "Plate based Technologies",
    "SPLiT-seq", "Other Technologies",
    "STRT-seq", "Plate based Technologies",
    "inDrop", "Drop based Technologies",
    "Smart-seq v4", "Smart seq",
    "ScaleBio single cell RNA sequencing", "Other Technologies"
  )
  
  
  disease_data_grouped <- tribble(
    ~disease, ~disease_groups,
    
    # Normal control
    "normal", "Normal",
    
    # Isolated Diseases
    "COVID-19", "COVID-19 related",
    "post-COVID-19 disorder", "COVID-19 related",
    "long COVID-19", "COVID-19 related",
    "glioblastoma", "Glioblastoma",
    "lung adenocarcinoma", "Lung Adenocarcinoma",
    "systemic lupus erythematosus", "Systemic Lupus Erythematosus",
    
    # Infectious and Immune-related Diseases (other than COVID-19)
    "Crohn disease", "Infectious and Immune-related Diseases",
    "Crohn ileitis", "Infectious and Immune-related Diseases",
    "pneumonia", "Infectious and Immune-related Diseases",
    "common variable immunodeficiency", "Infectious and Immune-related Diseases",
    "toxoplasmosis", "Infectious and Immune-related Diseases",
    "Plasmodium malariae malaria", "Infectious and Immune-related Diseases",
    "type 1 diabetes mellitus", "Infectious and Immune-related Diseases",
    "influenza", "Infectious and Immune-related Diseases",
    "chronic rhinitis", "Infectious and Immune-related Diseases",
    "periodontitis", "Infectious and Immune-related Diseases",
    "localized scleroderma", "Infectious and Immune-related Diseases",
    "lymphangioleiomyomatosis", "Infectious and Immune-related Diseases",
    "listeriosis", "Infectious and Immune-related Diseases",
    
    # Cancer (other than isolated cancers)
    "squamous cell lung carcinoma", "Cancer",
    "small cell lung carcinoma", "Cancer",
    "non-small cell lung carcinoma", "Cancer",
    "breast carcinoma", "Cancer",
    "breast cancer", "Cancer",
    "luminal B breast carcinoma", "Cancer",
    "luminal A breast carcinoma", "Cancer",
    "triple-negative breast carcinoma", "Cancer",
    "gastric cancer", "Cancer",
    "colorectal cancer", "Cancer",
    "colon sessile serrated adenoma/polyp", "Cancer",
    "follicular lymphoma", "Cancer",
    "B-cell acute lymphoblastic leukemia", "Cancer",
    "B-cell non-Hodgkin lymphoma", "Cancer",
    "acute myeloid leukemia", "Cancer",
    "acute promyelocytic leukemia", "Cancer",
    "plasma cell myeloma", "Cancer",
    "clear cell renal carcinoma", "Cancer",
    "nonpapillary renal cell carcinoma", "Cancer",
    "basal cell carcinoma", "Cancer",
    "colorectal neoplasm", "Cancer",
    "adenocarcinoma", "Cancer",
    "chromophobe renal cell carcinoma", "Cancer",
    "neuroendocrine carcinoma", "Cancer",
    "lung large cell carcinoma", "Cancer",
    "tongue cancer", "Cancer",
    "Wilms tumor", "Cancer",
    "pleomorphic carcinoma", "Cancer",
    "blastoma", "Cancer",
    
    # Neurodegenerative and Neurological Disorders
    "dementia", "Neurodegenerative and Neurological Disorders",
    "Alzheimer disease", "Neurodegenerative and Neurological Disorders",
    "Parkinson disease", "Neurodegenerative and Neurological Disorders",
    "amyotrophic lateral sclerosis", "Neurodegenerative and Neurological Disorders",
    "multiple sclerosis", "Neurodegenerative and Neurological Disorders",
    "Down syndrome", "Neurodegenerative and Neurological Disorders",
    "trisomy 18", "Neurodegenerative and Neurological Disorders",
    "frontotemporal dementia", "Neurodegenerative and Neurological Disorders",
    "temporal lobe epilepsy", "Neurodegenerative and Neurological Disorders",
    "Lewy body dementia", "Neurodegenerative and Neurological Disorders",
    "amyotrophic lateral sclerosis 26 with or without frontotemporal dementia", "Neurodegenerative and Neurological Disorders",
    
    # Respiratory Conditions
    "pulmonary fibrosis", "Respiratory Conditions",
    "respiratory system disorder", "Respiratory Conditions",
    "chronic obstructive pulmonary disease", "Respiratory Conditions",
    "cystic fibrosis", "Respiratory Conditions",
    "interstitial lung disease", "Respiratory Conditions",
    "hypersensitivity pneumonitis", "Respiratory Conditions",
    "non-specific interstitial pneumonia", "Respiratory Conditions",
    "aspiration pneumonia", "Respiratory Conditions",
    "pulmonary emphysema", "Respiratory Conditions",
    "pulmonary sarcoidosis", "Respiratory Conditions",
    
    # Cardiovascular Diseases
    "myocardial infarction", "Cardiovascular Diseases",
    "acute myocardial infarction", "Cardiovascular Diseases",
    "dilated cardiomyopathy", "Cardiovascular Diseases",
    "heart failure", "Cardiovascular Diseases",
    "arrhythmogenic right ventricular cardiomyopathy", "Cardiovascular Diseases",
    "congenital heart disease", "Cardiovascular Diseases",
    "non-compaction cardiomyopathy", "Cardiovascular Diseases",
    "cardiomyopathy", "Cardiovascular Diseases",
    "heart disorder", "Cardiovascular Diseases",
    
    # Metabolic and Other Disorders
    "type 2 diabetes mellitus", "Metabolic and Other Disorders",
    "chronic kidney disease", "Metabolic and Other Disorders",
    "digestive system disorder", "Metabolic and Other Disorders",
    "primary sclerosing cholangitis", "Metabolic and Other Disorders",
    "gastritis", "Metabolic and Other Disorders",
    "acute kidney failure", "Metabolic and Other Disorders",
    "tubular adenoma", "Metabolic and Other Disorders",
    "benign prostatic hyperplasia", "Metabolic and Other Disorders",
    "opiate dependence", "Metabolic and Other Disorders",
    "gingivitis", "Metabolic and Other Disorders",
    "hyperplastic polyp", "Metabolic and Other Disorders",
    "clonal hematopoiesis", "Metabolic and Other Disorders",
    "epilepsy", "Metabolic and Other Disorders",
    "age related macular degeneration 7", "Metabolic and Other Disorders",
    "kidney benign neoplasm", "Metabolic and Other Disorders",
    "malignant pancreatic neoplasm", "Metabolic and Other Disorders",
    "cataract", "Metabolic and Other Disorders",
    "macular degeneration", "Metabolic and Other Disorders",
    "hydrosalpinx", "Metabolic and Other Disorders",
    "tubulovillous adenoma", "Metabolic and Other Disorders",
    "gastric intestinal metaplasia", "Metabolic and Other Disorders",
    "Barrett esophagus", "Metabolic and Other Disorders",
    
    # Other Diseases
    "injury", "Other Diseases",
    "anencephaly", "Other Diseases",
    "primary biliary cholangitis", "Other Diseases",
    "keloid", "Other Diseases",
    "kidney oncocytoma", "Other Diseases",
    "respiratory failure", "Other Diseases",
    "pilocytic astrocytoma", "Other Diseases"
  )
  
  temp_path = tempdir()
  system(glue("~/bin/rclone copy box_adelaide:/minh_immune_map_disease/disease_data_grouped_further.csv {temp_path}/"))
  
  disease_data_grouped =
    disease_data_grouped |>
    left_join(
      read_csv(glue("{temp_path}/disease_data_grouped_further.csv")) |>
        rename(disease_groups_further = disease_groups)
    ) |>
    mutate(disease_groups = if_else(!disease_groups_further |> is.na(), disease_groups_further, disease_groups)) |>
    select(disease,  disease_groups)
  
  tbl |> 
    
    # TISSUE
    select(-any_of("tissue_groups")) |> 
    left_join(tissue_grouped, copy=TRUE) |> 
    
    # TECH
    left_join(assay_data_grouped, copy=TRUE) |> 
    

    
    
    # TEMPORARY. de-group pancreas and liver
    mutate(tissue_groups = case_when(
      
      tissue %in% c("gallbladder") ~ "gallbladder",
      tissue %in% c("pancreas", "exocrine pancreas") ~ "pancreas",
      tissue %in% c("liver", "caudate lobe of liver", "hepatic cecum" ) ~ "liver",
      TRUE ~ tissue_groups
    )) |> 
    
    # SEX edit
    mutate(sex = if_else(sex |> is.na(), "unknown", sex)) |> 
    
    # Age
    filter(age_days > 365) |> 
    mutate(age_years = age_days / 365) |> 
    mutate(age_bin = dplyr::case_when(
      age_years < 3 ~ "Infancy",
      age_years < 12 ~ "Childhood",
      age_years < 20 ~ "Adolescence",
      age_years < 40 ~ "Young Adulthood",
      age_years < 50 ~ "Middle Age",
      age_years < 60 ~ "Senior_50",
      age_years < 70 ~ "Senior_60",
      age_years >= 70 ~ "Senior_70",
      TRUE ~ NA_character_
    )) |> 
    mutate(age_decade = ceiling(age_years/10) |> as.character()) |> 
    
    # left_join(age_bin_table, copy=TRUE) |> 
    
    # ETHNICITY
    left_join(ethnicity_grouped, copy=TRUE) |> 
    
    dplyr::select(
      sample_id, donor_id, dataset_id, title, collection_id, age_days, age_bin, age_decade, sex, 
      ethnicity_groups, tissue_groups, tissue, assay_groups, cell_type_unified_ensemble,
      cell_type, disease_groups, is_immune
    ) |> 
    as_tibble() |> 
    
    # Center based on adolescence
    mutate(age_days_scaled = age_days  |> scale(center = 50*365) |> as.numeric()) 
  
}     

ethnicity_grouped <- tribble(
  ~self_reported_ethnicity, ~ethnicity_groups,
  "unknown", "Other/Unknown",
  "European", "European",
  "Korean", "East Asian",
  "Asian", "East Asian",
  "Japanese", "Japanese",
  "African American", "African",
  "Hispanic or Latin American", "Hispanic/Latin American",
  "Singaporean Chinese", "East Asian",
  "Han Chinese", "East Asian",
  "Singaporean Indian", "South Asian",
  "Singaporean Malay", "Other/Unknown",
  "British", "European",
  "African", "African",
  "South Asian", "South Asian",
  "European American", "European",
  "East Asian", "East Asian",
  "American", "Other/Unknown",
  "African American or Afro-Caribbean", "African",
  "Oceanian", "Native American & Pacific Islander",
  "Jewish Israeli", "Middle Eastern & North African",
  "Chinese", "East Asian",
  "South East Asian", "Other/Unknown",
  "Greater Middle Eastern  (Middle Eastern or North African or Persian)", "Middle Eastern & North African",
  "Native American", "Native American & Pacific Islander",
  "Pacific Islander", "Native American & Pacific Islander",
  "Finnish", "European",
  "Bangladeshi", "South Asian",
  "Native American,Hispanic or Latin American", "Hispanic/Latin American",
  "Irish", "European",
  "Iraqi", "Middle Eastern & North African",
  "European,Asian", "European"
)

disease_data_grouped <- tribble(
  ~disease, ~disease_groups,
  
  # Normal control
  "normal", "Normal",
  
  # Isolated Diseases
  "COVID-19", "COVID-19 related",
  "post-COVID-19 disorder", "COVID-19 related",
  "long COVID-19", "COVID-19 related",
  "glioblastoma", "Glioblastoma",
  "lung adenocarcinoma", "Lung Adenocarcinoma",
  "systemic lupus erythematosus", "Systemic Lupus Erythematosus",
  
  # Infectious and Immune-related Diseases (other than COVID-19)
  "Crohn disease", "Infectious and Immune-related Diseases",
  "Crohn ileitis", "Infectious and Immune-related Diseases",
  "pneumonia", "Infectious and Immune-related Diseases",
  "common variable immunodeficiency", "Infectious and Immune-related Diseases",
  "toxoplasmosis", "Infectious and Immune-related Diseases",
  "Plasmodium malariae malaria", "Infectious and Immune-related Diseases",
  "type 1 diabetes mellitus", "Infectious and Immune-related Diseases",
  "influenza", "Infectious and Immune-related Diseases",
  "chronic rhinitis", "Infectious and Immune-related Diseases",
  "periodontitis", "Infectious and Immune-related Diseases",
  "localized scleroderma", "Infectious and Immune-related Diseases",
  "lymphangioleiomyomatosis", "Infectious and Immune-related Diseases",
  "listeriosis", "Infectious and Immune-related Diseases",
  
  # Cancer (other than isolated cancers)
  "squamous cell lung carcinoma", "Cancer",
  "small cell lung carcinoma", "Cancer",
  "non-small cell lung carcinoma", "Cancer",
  "breast carcinoma", "Cancer",
  "breast cancer", "Cancer",
  "luminal B breast carcinoma", "Cancer",
  "luminal A breast carcinoma", "Cancer",
  "triple-negative breast carcinoma", "Cancer",
  "gastric cancer", "Cancer",
  "colorectal cancer", "Cancer",
  "colon sessile serrated adenoma/polyp", "Cancer",
  "follicular lymphoma", "Cancer",
  "B-cell acute lymphoblastic leukemia", "Cancer",
  "B-cell non-Hodgkin lymphoma", "Cancer",
  "acute myeloid leukemia", "Cancer",
  "acute promyelocytic leukemia", "Cancer",
  "plasma cell myeloma", "Cancer",
  "clear cell renal carcinoma", "Cancer",
  "nonpapillary renal cell carcinoma", "Cancer",
  "basal cell carcinoma", "Cancer",
  "colorectal neoplasm", "Cancer",
  "adenocarcinoma", "Cancer",
  "chromophobe renal cell carcinoma", "Cancer",
  "neuroendocrine carcinoma", "Cancer",
  "lung large cell carcinoma", "Cancer",
  "tongue cancer", "Cancer",
  "Wilms tumor", "Cancer",
  "pleomorphic carcinoma", "Cancer",
  "blastoma", "Cancer",
  
  # Neurodegenerative and Neurological Disorders
  "dementia", "Neurodegenerative and Neurological Disorders",
  "Alzheimer disease", "Neurodegenerative and Neurological Disorders",
  "Parkinson disease", "Neurodegenerative and Neurological Disorders",
  "amyotrophic lateral sclerosis", "Neurodegenerative and Neurological Disorders",
  "multiple sclerosis", "Neurodegenerative and Neurological Disorders",
  "Down syndrome", "Neurodegenerative and Neurological Disorders",
  "trisomy 18", "Neurodegenerative and Neurological Disorders",
  "frontotemporal dementia", "Neurodegenerative and Neurological Disorders",
  "temporal lobe epilepsy", "Neurodegenerative and Neurological Disorders",
  "Lewy body dementia", "Neurodegenerative and Neurological Disorders",
  "amyotrophic lateral sclerosis 26 with or without frontotemporal dementia", "Neurodegenerative and Neurological Disorders",
  
  # Respiratory Conditions
  "pulmonary fibrosis", "Respiratory Conditions",
  "respiratory system disorder", "Respiratory Conditions",
  "chronic obstructive pulmonary disease", "Respiratory Conditions",
  "cystic fibrosis", "Respiratory Conditions",
  "interstitial lung disease", "Respiratory Conditions",
  "hypersensitivity pneumonitis", "Respiratory Conditions",
  "non-specific interstitial pneumonia", "Respiratory Conditions",
  "aspiration pneumonia", "Respiratory Conditions",
  "pulmonary emphysema", "Respiratory Conditions",
  "pulmonary sarcoidosis", "Respiratory Conditions",
  
  # Cardiovascular Diseases
  "myocardial infarction", "Cardiovascular Diseases",
  "acute myocardial infarction", "Cardiovascular Diseases",
  "dilated cardiomyopathy", "Cardiovascular Diseases",
  "heart failure", "Cardiovascular Diseases",
  "arrhythmogenic right ventricular cardiomyopathy", "Cardiovascular Diseases",
  "congenital heart disease", "Cardiovascular Diseases",
  "non-compaction cardiomyopathy", "Cardiovascular Diseases",
  "cardiomyopathy", "Cardiovascular Diseases",
  "heart disorder", "Cardiovascular Diseases",
  
  # Metabolic and Other Disorders
  "type 2 diabetes mellitus", "Metabolic and Other Disorders",
  "chronic kidney disease", "Metabolic and Other Disorders",
  "digestive system disorder", "Metabolic and Other Disorders",
  "primary sclerosing cholangitis", "Metabolic and Other Disorders",
  "gastritis", "Metabolic and Other Disorders",
  "acute kidney failure", "Metabolic and Other Disorders",
  "tubular adenoma", "Metabolic and Other Disorders",
  "benign prostatic hyperplasia", "Metabolic and Other Disorders",
  "opiate dependence", "Metabolic and Other Disorders",
  "gingivitis", "Metabolic and Other Disorders",
  "hyperplastic polyp", "Metabolic and Other Disorders",
  "clonal hematopoiesis", "Metabolic and Other Disorders",
  "epilepsy", "Metabolic and Other Disorders",
  "age related macular degeneration 7", "Metabolic and Other Disorders",
  "kidney benign neoplasm", "Metabolic and Other Disorders",
  "malignant pancreatic neoplasm", "Metabolic and Other Disorders",
  "cataract", "Metabolic and Other Disorders",
  "macular degeneration", "Metabolic and Other Disorders",
  "hydrosalpinx", "Metabolic and Other Disorders",
  "tubulovillous adenoma", "Metabolic and Other Disorders",
  "gastric intestinal metaplasia", "Metabolic and Other Disorders",
  "Barrett esophagus", "Metabolic and Other Disorders",
  
  # Other Diseases
  "injury", "Other Diseases",
  "anencephaly", "Other Diseases",
  "primary biliary cholangitis", "Other Diseases",
  "keloid", "Other Diseases",
  "kidney oncocytoma", "Other Diseases",
  "respiratory failure", "Other Diseases",
  "pilocytic astrocytoma", "Other Diseases"
)

disease_data_grouped_coarse <- tribble(
  ~ disease_groups, ~disease_groups_coarse,
  "Other Diseases", "Other",
  "Glioblastoma", "Glioblastoma",
  "Cardiovascular Diseases", "Cardiovascular",
  "Normal", "Normal",
  "COVID-19 related", "COVID-19",
  "Systemic Lupus Erythematosus", "Lupus",
  "Cancer", "Cancer",
  "Metabolic and Other Disorders", "Metabolic/Other",
  "Lung Adenocarcinoma", "Lung AdenoCA",
  "Respiratory Conditions", "Respiratory",
  "Infectious and Immune-related Diseases", "Infectious/Immune",
  "Neurodegenerative and Neurological Disorders", "Neuro/Neurodeg"
)

shorten_technology_label <- tribble(
  ~ assay, ~ assay_shorten,
  "BD Rhapsody Targeted mRNA", "BD Rhapsody Panel",
  "inDrop", "inDrop",
  "GEXSCOPE technology", "GEXSCOPE",
  "BD Rhapsody Whole Transcriptome Analysis", "BD Rhapsody WTA",
  "10x 3' v2", "10x 3' v2",
  "MARS-seq", "MARS-seq",
  "microwell-seq", "Microwell",
  "Seq-Well", "Seq-Well",
  "Smart-seq2", "Smart-seq2",
  "STRT-seq", "STRT-seq",
  "10x 3' v3", "10x 3' v3",
  "ScaleBio single cell RNA sequencing", "ScaleBio",
  "TruDrop", "TruDrop",
  "sci-RNA-seq", "sci-RNA-seq",
  "Seq-Well S3", "Seq-Well S3",
  "CEL-seq2", "CEL-seq2",
  "10x 5' transcription profiling", "10x 5' TP",
  "SPLiT-seq", "SPLiT-seq",
  "10x 5' v1", "10x 5' v1",
  "Smart-seq v4", "Smart-seq4",
  "10x 5' v2", "10x 5' v2",
  "10x 3' v1", "10x 3' v1",
  "Drop-seq", "Drop-seq",
  "10x 3' transcription profiling", "10x 3' TP"
)
