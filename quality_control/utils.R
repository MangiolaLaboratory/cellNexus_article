get_tissue_grouped <- function(tissue) {

  # Sex-specific tissue entries (gonad differs by donor sex)
  sex_specific <- tibble::tribble(
    ~tissue,  ~sex,      ~tissue_groups,
    "gonad",  "male",    "male reproductive system (other)",
    "gonad",  "female",  "female reproductive system"
  )

  list(
    # Respiratory System
    "respiratory system" = c(
      "lung", "lung parenchyma", "alveolus of lung", "bronchus",
      "respiratory airway", "pleura", "pleural effusion", "middle lobe of right lung",
      "upper lobe of left lung", "lower lobe of left lung", "upper lobe of right lung",
      "lower lobe of right lung", "lingula of left lung", "right lung", "left lung"
    ),
    trachea = c("epithelium of trachea", "trachea"),

    # Cardiovascular System
    "cardiovascular system" = c(
      "heart", "heart left ventricle", "heart right ventricle", "cardiac ventricle",
      "cardiac atrium", "right cardiac atrium", "left cardiac atrium", "apex of heart",
      "aorta", "coronary artery",
      "venous blood", "anterior wall of left ventricle", "myocardium", "interventricular septum", "ventricular tissue", "basal zone of heart"
    ),
    vasculature = c("kidney blood vessel", "artery", "vein", "vasculature", "mesenteric artery"),
    # Umbilical Cord Blood
    "umbilical cord blood" = "umbilical cord blood",

    # Oesophagus
    "oesophagus" = c(
      "esophagus", "lower esophagus", "esophagus muscularis mucosa",
      "submucosal esophageal gland"
    ),

    # Stomach
    "stomach" = c(
      "stomach", "body of stomach", "cardia of stomach"
    ),

    # Small Intestine
    "small intestine" = c(
      "small intestine", "duodenum", "jejunum", "ileum"
    ),

    # Large Intestine
    "large intestine" = c(
      "large intestine", "colon", "left colon", "right colon",
      "sigmoid colon", "descending colon", "transverse colon",
      "ascending colon", "hepatic flexure of colon", "caecum",
      "rectum", "appendix", "vermiform appendix"
    ),

    # Digestive System (General)
    "digestive system (general)" = c(
      "intestine", "hindgut"
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
      "skin of pes", "axilla", "nose skin", "scalp"
    ),

    # Gastrointestinal Accessory Organs
    "gallbladder" = "gallbladder",

    # Gastrointestinal Accessory Organs
    "pancreas" = c("pancreas", "exocrine pancreas"),

    # Gastrointestinal Accessory Organs
    "liver" = c("liver", "caudate lobe of liver", "hepatic cecum"),

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
      "decidua", "decidua basalis", "isthmus of fallopian tube"
    ),
    "ovary" = "ovary",

    # Extraembryonic and Placental Structures
    "extraembryonic and placental structures" = c(
      "placenta", "yolk sac"
    ),

    # Male Reproductive System
    "male reproductive system (other)" = c(
      "testis"
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

    # Epithelium and Mucosal Tissues
    "epithelium and mucosal tissues" = c(
      "epithelium of small intestine", "epithelium of esophagus", "caecum epithelium",
      "jejunal epithelium", "ileal epithelium", "colonic epithelium",
      "submucosa of ascending colon", "submucosa of ileum", "lamina propria",
      "lamina propria of large intestine", "lamina propria of small intestine",
      "mucosa", "mucosa of colon", "lamina propria of mucosa of colon"
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
    enframe(name = "tissue_groups") |>
    distinct() |>
    unnest(value) |>
    dplyr::rename(tissue = value) |>
    mutate(sex = NA_character_) |>
    bind_rows(sex_specific) |>
    select(tissue, sex, tissue_groups)

}

tissue_group_conversion_tbl <- tibble::tibble(
  tissue_groups = c(
    "small intestine",
    "spleen",
    "large intestine",
    "sensory-related structures",
    "adipose tissue",
    "vasculature",
    "male reproductive system (other)",
    "umbilical cord blood",
    "stomach",
    "gallbladder",
    "liver",
    "pancreas",
    "extraembryonic and placental structures",
    "gastrointestinal accessory organs",
    "nasal, oral, and pharyngeal regions",
    "connective tissue",
    "muscular system (skeletal muscles)",
    "breast",
    "prostate",
    "digestive tract junctions and connections",
    "respiratory system",
    "lymphatic system",
    "epithelium and mucosal tissues",
    "oesophagus",
    "trachea",
    "ovary",
    "peritoneal and abdominal cavity structures",
    "blood",
    "brainstem and cerebellar structures",
    "integumentary system (skin)",
    "miscellaneous glands",
    "thymus",
    "limbic and basal systems",
    "endocrine system",
    "renal system",
    "general brain and major structures",
    "cerebral lobes and cortical areas",
    "bone marrow",
    "female reproductive system",
    "cardiovascular system",
    "digestive system (general)"
  ),
  tissue_groups_short = c(
    "Small Intestine",
    "Spleen",
    "Large Intestine",
    "Sensory",
    "Adipose",
    "Vasculature",
    "Male Repro.",
    "Cord Blood",
    "Stomach",
    "Gallbladder",
    "Liver",
    "Pancreas",
    "Extraembryonic/Placental",
    "GI Accessory",
    "Nasal/Oral/Pharyngeal",
    "Connective",
    "Muscle",
    "Breast",
    "Prostate",
    "GI Junctions",
    "Respiratory",
    "Lymphatic",
    "Epithelium/Mucosa",
    "Oesophagus",
    "Trachea",
    "Ovary",
    "Peritoneal/Abdominal",
    "Blood",
    "Brainstem/Cerebellum",
    "Skin",
    "Misc. Glands",
    "Thymus",
    "Limbic/Basal",
    "Endocrine",
    "Renal",
    "Brain (General)",
    "Cerebral Lobes",
    "Bone Marrow",
    "Female Repro.",
    "Cardiovascular",
    "Digestive"
  )
)