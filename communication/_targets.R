library(targets)
library(tidyverse)
library(targets)
library(tarchetypes)
library(glue)
library(qs)
library(crew)
library(crew.cluster)
new_elastic <- function(name, mem_gb, time_min, workers, crashes_max, 
    backup = NULL) {
    crew_controller_slurm(name = name, workers = workers, crashes_max = crashes_max, 
        seconds_idle = 30, options_cluster = crew_options_slurm(memory_gigabytes_required = mem_gb, 
            cpus_per_task = 8, time_minutes = time_min), backup = backup)
}
elastic_160 <- new_elastic("elastic_160", 160, 60 * 24, workers = 8, 
    crashes_max = 2)
elastic_80 <- new_elastic("elastic_80", 80, 60 * 4, workers = 16, 
    crashes_max = 1, backup = elastic_160)
elastic_40 <- new_elastic("elastic_40", 40, 60 * 4, workers = 24, 
    crashes_max = 1, backup = elastic_80)
elastic_20 <- new_elastic("elastic_20", 20, 60 * 4, workers = 32, 
    crashes_max = 1, backup = elastic_40)
elastic_10 <- new_elastic("elastic_10", 10, 60 * 4, workers = 48, 
    crashes_max = 1, backup = elastic_20)
elastic_5 <- new_elastic("elastic_5", 5, 60 * 4, workers = 64, 
    crashes_max = 6, backup = elastic_10)
controllers <- crew_controller_group(elastic_5, elastic_10, elastic_20, 
    elastic_40, elastic_80, elastic_160)
tar_option_set(memory = "transient", garbage_collection = 100, 
    storage = "worker", retrieval = "worker", error = "continue", 
    workspace_on_error = TRUE, format = "qs", workspaces = "estimates_chunk_633f6596569029e6", 
    debug = "estimates_chunk_8da8c24dbc0993bb", controller = controllers)
remove_unwanted_effect = function(fit, newdata, robust = FALSE, 
    correct_by_offset = T, re_formula = ~0) {
    fitted_residuals = residuals(fit, robust = robust, summary = FALSE)
    if (correct_by_offset) 
        fitted_residuals = sweep(fitted_residuals, 2, exp(fit$data$offset), 
            FUN = "/")
    fitted_values_ethnicity <- fitted(fit, newdata = newdata, 
        re_formula = re_formula, summary = FALSE, offset = 0)
    adjusted_counts = fitted_values_ethnicity + fitted_residuals
    fitted_residuals_tbl = as_tibble(posterior_summary(fitted_residuals, 
        robust = robust))
    colnames(fitted_residuals_tbl) = paste0("residuals___", colnames(fitted_residuals_tbl))
    fitted_values_ethnicity_tbl = as_tibble(posterior_summary(fitted_values_ethnicity, 
        robust = robust))
    colnames(fitted_values_ethnicity_tbl) = paste0("fitted___", 
        colnames(fitted_values_ethnicity_tbl))
    adjusted_counts_tbl = as_tibble(posterior_summary(adjusted_counts, 
        robust = robust))
    colnames(adjusted_counts_tbl) = paste0("adjusted___", colnames(adjusted_counts_tbl))
    bind_cols(bind_cols(adjusted_counts_tbl, fitted_residuals_tbl), 
        fitted_values_ethnicity_tbl)
}
remove_unwanted_effect_new = function(fit, newdata, robust = FALSE, 
    correct_by_offset = T, re_formula = ~0) {
    fitted_residuals = predictive_error(fit, robust = robust, 
        summary = FALSE, offset = 0)
    if (correct_by_offset) 
        fitted_residuals = sweep(fitted_residuals, 2, exp(fit$data$offset), 
            FUN = "/")
    fitted_values_ethnicity <- posterior_epred(fit, newdata = newdata, 
        re_formula = re_formula, offset = 0)
    adjusted_counts = fitted_values_ethnicity + fitted_residuals
    fitted_residuals_tbl = as_tibble(posterior_summary(fitted_residuals, 
        robust = robust))
    colnames(fitted_residuals_tbl) = paste0("residuals___", colnames(fitted_residuals_tbl))
    fitted_values_ethnicity_tbl = as_tibble(posterior_summary(fitted_values_ethnicity, 
        robust = robust))
    colnames(fitted_values_ethnicity_tbl) = paste0("fitted___", 
        colnames(fitted_values_ethnicity_tbl))
    adjusted_counts_tbl = as_tibble(posterior_summary(adjusted_counts, 
        robust = robust))
    colnames(adjusted_counts_tbl) = paste0("adjusted___", colnames(adjusted_counts_tbl))
    bind_cols(bind_cols(adjusted_counts_tbl, fitted_residuals_tbl), 
        fitted_values_ethnicity_tbl)
}
get_adjusted_matrix = function(summary_df, column_adjusted) {
    column_adjusted = enquo(column_adjusted)
    m = Matrix::Matrix(as(tidybulk:::as_matrix(pivot_wider(select(unnest(summary_df, 
        !!column_adjusted), .feature, adjusted___Estimate, sample_id), 
        names_from = sample_id, values_from = adjusted___Estimate), 
        rownames = ".feature"), "sparseMatrix"), sparse = T)
    max_rm_infinite = quantile(m[!is.infinite(m)], 0.999)
    m[m > max_rm_infinite] = max_rm_infinite
    m[m < 0] = 0
    return(m)
}
check_rclone_installation = function() {
    rclone_path <- path.expand("~/bin/rclone")
    if (!file.exists(rclone_path)) {
        stop("rclone was not found in the expected location '~/bin/rclone'.")
    }
}
source("https://gist.githubusercontent.com/stemangiola/8fe6c45b79dd95a200c0fd2314ec57d0/raw/14c753da54df886c29a02ac5fbf464684002956c/gistfile1.txt")
monotone_ageing_hypothesis_testing = function() {
    cell_type <- cell_type_df$cell_type[[1]] %>% make.names()
    run_path = glue("{target_path}/V1_{cell_type}/")
    pseudobulk_sample_id = qs_read(glue("{run_path}_targets/objects/pseudobulk_sample_id"))
    plan(multisession, workers = 16)
    age_sex_summary <- read_delim(glue("{run_path}_targets/meta/meta"), 
        delim = "|") %>% filter(name %>% str_starts("estimates_chunk")) %>% 
        filter(type == "branch" & is.na(error)) %>% pull(name) %>% 
        future_map_dfr(function(x, run_path) {
            tmp = qs_read(glue("{run_path}_targets/objects/{x}"))
            res = tmp$brms_fit[[1]]$fit %>% summary() %>% as.data.frame() %>% 
                filter(str_detect(rownames(.), "age_bin|sex") & 
                  str_detect(rownames(.), "^(b_|r_)")) %>% select(summary.mean) %>% 
                t %>% as.data.frame()
            rownames(res) = tmp$.feature
            return(res)
        }, run_path = run_path, .progress = TRUE)
    saveRDS(age_sex_summary, file = glue("{run_path}age_sex_summary.rds"))
}
build_contrast <- function(k, bins, minimum_difference = 0) {
    if (minimum_difference == 0) 
        sign = "="
    else if (minimum_difference < 0) 
        sign = "<"
    else if (minimum_difference > 0) 
        sign = ">"
    younger <- bins[1:k]
    older <- bins[(k + 1):length(bins)]
    glue::glue("({paste(older,   collapse = ' + ')})/{length(older)} - ", 
        "({paste(younger, collapse = ' + ')})/{length(younger)} {sign} {minimum_difference}")
}
hypothesis_abs <- function(fit, k, bins, abs_threshold = 0.2, 
    scope = "coef", group = "", re_formula = NA, resp = NULL) {
    h_pos <- build_contrast(k, bins, abs_threshold)
    h_neg <- stringr::str_replace(build_contrast(k, bins, -abs_threshold), 
        ">", "<")
    tbl_pos <- dplyr::mutate(tibble::as_tibble(brms::hypothesis(fit, 
        h_pos, scope = scope, group = group, re_formula = re_formula, 
        resp = resp)$hypothesis), direction = "positive", prob = Post.Prob)
    tbl_neg <- dplyr::mutate(tibble::as_tibble(brms::hypothesis(fit, 
        h_neg, scope = scope, group = group, re_formula = re_formula, 
        resp = resp)$hypothesis), direction = "negative", prob = Post.Prob)
    long_tbl <- dplyr::bind_rows(tbl_pos, tbl_neg)
    if (!"Group" %in% colnames(long_tbl)) 
        long_tbl = mutate(long_tbl, Group = "population")
    out <- dplyr::transmute(dplyr::ungroup(dplyr::slice_max(dplyr::mutate(dplyr::group_by(long_tbl, 
        Group), prob_abs = sum(prob), p_two_tail = 2 * min(prob), 
        p_two_tail = pmin(p_two_tail, 1)), prob, n = 1, with_ties = FALSE)), 
        Group, Hypothesis, Estimate, CI.Lower, CI.Upper, Post.Prob = p_two_tail, 
        direction)
    out
}
fit_to_age_monotonic_changes = function(fit) {
    if (is.null(fit)) 
        return(NULL)
    full_bins <- c("age_binInfancy", "age_binChildhood", "age_binAdolescence", 
        "age_binYoungAdulthood", "age_binMiddleAge", "age_binSenior_60", 
        "age_binSenior_70")
    vars <- brms::variables(fit)
    present_bins <- full_bins[paste0("b_", full_bins) %in% vars]
    if (length(present_bins) < 2) {
        return(tibble())
    }
    purrr::map_dfr(seq_len(length(present_bins) - 1), function(k) {
        h_txt <- build_contrast(k, present_bins)
        h_tot = hypothesis_abs(fit, k, present_bins, abs_threshold = 0.2, 
            scope = "coef", group = "tissue_groups")
        h_tot = h_tot %>% tibble::as_tibble() %>% dplyr::transmute(component = "total", 
            tissue = Group, split_after = present_bins[k], younger_bins = paste(present_bins[1:k], 
                collapse = ","), older_bins = paste(present_bins[(k + 
                1):length(present_bins)], collapse = ","), estimate = Estimate, 
            ci_lower = CI.Lower, ci_upper = CI.Upper, post_prob = Post.Prob, 
            Hypothesis = Hypothesis, direction = direction)
        h_fix <- hypothesis_abs(fit, k, present_bins, abs_threshold = 0.2, 
            scope = "standard")
        h_fix = h_fix %>% tibble::as_tibble() %>% dplyr::transmute(component = "fixed", 
            tissue = "population", split_after = present_bins[k], 
            younger_bins = paste(present_bins[1:k], collapse = ","), 
            older_bins = paste(present_bins[(k + 1):length(present_bins)], 
                collapse = ","), estimate = Estimate, ci_lower = CI.Lower, 
            ci_upper = CI.Upper, post_prob = Post.Prob, Hypothesis = Hypothesis, 
            direction = direction)
        dplyr::bind_rows(h_tot, h_fix)
    })
}
side_expr <- function(bins) paste(bins, collapse = " + ")
has_param <- function(expr) expr != "0"
safe_hyp <- purrr::possibly(brms::hypothesis, NULL)
prepare_database = function(tbl, ethnicity_imputed) {
    inner_join(distinct(tbl, sample_id, source, target, pathway_name, 
        interaction_weight, interaction_count, pathway_prob, 
        annotation), mutate(tidybulk:::.resolve_complete_confounders_of_non_interest_df(mutate(left_join(distinct(edit_covariates(distinct(get_metadata(), 
        sample_id, donor_id, dataset_id, tissue, age_days, sex, 
        self_reported_ethnicity, disease, assay, title, collection_id, 
        cell_type_unified_ensemble, cell_type, is_immune)), sample_id, 
        donor_id, dataset_id, tissue_groups, age_days, sex, ethnicity_groups, 
        disease_groups, assay_groups, age_decade, age_bin), ethnicity_imputed, 
        by = join_by(sample_id, ethnicity_groups)), ethnicity_groups_imputed = if_else(is.na(ethnicity_groups_imputed), 
        ethnicity_groups, ethnicity_groups_imputed)), dataset_id, 
        assay_groups, disease_groups), ethnicity_groups_imputed = fct_relevel(ethnicity_groups_imputed, 
        "European"), assay_groups___altered = fct_relevel(assay_groups___altered, 
        "10x Genomics 3"), disease_groups___altered = fct_relevel(disease_groups___altered, 
        "Normal"), age_bin = fct_relevel(age_bin, "Senior_50"), 
        age_decade = fct_relevel(age_decade, "5")), copy = TRUE)
}
get_pairs_to_consider = function(all_cell_types) {
    myeloid_lymphoid_pairs <- tribble(~source, ~target, ~mechanism, 
        ~doi, "cdc", "cd4 naive", "Antigen presentation + CD80/86–CD28 costimulation (priming)", 
        "10.1038/32588", "cdc", "cd4 tcm", "Recall priming / re-stimulation by cDC", 
        "10.1038/32588", "cdc", "cd4 tem", "Effector re-stimulation by cDC", 
        "10.1038/32588", "cdc", "cd4 fh em", "ICOSL/IL-6-driven Tfh differentiation by APCs incl. DCs", 
        "10.1016/j.immuni.2014.10.004", "cdc", "cd4 th1 em", 
        "IL-12 skews Th1 during DC-T priming", "10.1084/jem.184.2.741", 
        "cdc", "cd4 th2 em", "DC-guided Th2 under appropriate cues (e.g., OX40L/TSLP contexts)", 
        "10.1038/32588", "cdc", "cd4 th17 em", "IL-23/IL-6 axis promotes Th17 from DC priming", 
        "10.1084/jem.20041257", "cdc", "cd4 th1/th17 em", "Mixed Th1/Th17 polarisation from DC cytokine milieu", 
        "10.1084/jem.20041257", "cdc", "treg", "Specialised DCs (e.g., CD103+ intestinal DCs) induce peripheral FoxP3+ Tregs", 
        "10.1093/intimm/dxae042", "cdc", "cd8 naive", "Cross-presentation primes CD8+ T cells", 
        "10.1038/35100512", "cdc", "cd8 tcm", "Recall responses via cross-presenting DCs", 
        "10.1038/35100512", "cdc", "cd8 tem", "Effector re-stimulation via cross-presentation", 
        "10.1038/35100512", "cdc", "cytotoxic", "Licensing/priming of cytotoxic T cells by cross-presentation", 
        "10.1038/35100512", "cdc", "tgd", "DC antigen presentation & costimulation supports γδ T activation", 
        "10.1038/32588", "pdc", "nk", "Type I IFN from pDC activates NK cells", 
        "10.1126/science.284.5421.1835", "pdc", "b naive", "pDC (±viral trigger) induce B-cell activation/differentiation via IFN-I + IL-6", 
        "10.1016/S1074-7613(03)00208-5", "pdc", "b memory", "pDC (±viral trigger) induce B-cell activation/differentiation via IFN-I + IL-6", 
        "10.1016/S1074-7613(03)00208-5", "pdc", "plasma", "pDC drive plasmablast/plasma-cell differentiation", 
        "10.1016/S1074-7613(03)00208-5", "cd14 mono", "nk", "IL-12/IL-18 family from myeloid cells activates NK (myeloid–NK cytokine crosstalk)", 
        "10.3389/fimmu.2021.739220", "cd16 mono", "nk", "IL-12/IL-18 family from myeloid cells activates NK (myeloid–NK cytokine crosstalk)", 
        "10.3389/fimmu.2021.739220", "monocytic", "nk", "IL-12/IL-18 family from myeloid cells activates NK (myeloid–NK cytokine crosstalk)", 
        "10.3389/fimmu.2021.739220", "cd14 mono", "mait", "Cytokine-driven (IL-12/IL-18) MAIT activation by monocytes/APCs", 
        "10.3389/fimmu.2020.01014", "cd16 mono", "mait", "Cytokine-driven (IL-12/IL-18) MAIT activation by monocytes/APCs", 
        "10.3389/fimmu.2020.01014", "monocytic", "mait", "Cytokine-driven (IL-12/IL-18) MAIT activation by monocytes/APCs", 
        "10.3389/fimmu.2020.01014", "macrophage", "cd4 th1 em", 
        "Macrophage IL-12 promotes Th1 differentiation/maintenance", 
        "10.1084/jem.184.2.741", "macrophage", "cd4 th17 em", 
        "Macrophage/monocyte IL-23/IL-6 supports Th17", "10.1084/jem.20041257", 
        "macrophage", "cd8 naive", "Myeloid APCs (incl. macrophages) prime/boost CD8 under some settings", 
        "10.3389/fimmu.2013.00389", "macrophage", "cd8 tem", 
        "Myeloid APCs (incl. macrophages) prime/boost CD8 under some settings", 
        "10.3389/fimmu.2013.00389", "macrophage", "cytotoxic", 
        "Boosting/maintenance of cytotoxic CD8 via macrophage APC function", 
        "10.3389/fimmu.2013.00389", "macrophage", "nk", "Reciprocal cytokines (IL-12/IL-18 from myeloid ↔ IFN-γ/TNF-α from NK)", 
        "10.3389/fimmu.2021.739220", "macrophage", "nkt", "CD1d+ macrophages present glycolipids to iNKT (APC role)", 
        "10.1016/j.coi.2007.03.014", "macrophage", "mait", "MR1/cytokine-dependent MAIT activation by macrophages", 
        "10.3389/fimmu.2020.01014", "granulocyte", "cd4 tcm", 
        "Neutrophils can acquire APC features for memory CD4+", 
        "10.1182/blood-2016-10-744441", "granulocyte", "cd4 tem", 
        "Neutrophils can acquire APC features for memory CD4+", 
        "10.1182/blood-2016-10-744441", "granulocyte", "nk", 
        "Neutrophil–NK functional crosstalk in inflammation/cancer", 
        "10.3389/fimmu.2020.570380", "mast", "cd4 tem", "Mast cells enhance T-cell activation via TNF/costims incl. OX40L", 
        "10.3389/fimmu.2015.00394", "mast", "treg", "Bidirectional mast cell–Treg modulation (OX40–OX40L/cytokines)", 
        "10.3389/fimmu.2015.00394", "mast", "tgd", "Direct immune synapse with γδ T cells (antiviral context)", 
        "10.1172/JCI122530", "mast", "b naive", "Mast-cell mediators (e.g., IL-6) promote B-cell proliferation/differentiation", 
        "10.1182/blood-2009-10-250126", "mast", "b memory", "Mast-cell mediators (e.g., IL-6) promote B-cell proliferation/differentiation", 
        "10.1182/blood-2009-10-250126", "ilc", "macrophage", 
        "ILC2-derived IL-13 drives M2-like macrophage polarisation", 
        "10.1155/2020/5018975", "nk", "cdc", "DC↔NK licensing/activation loops (IL-12/IL-18, IL-15 trans-presentation)", 
        "10.3389/fimmu.2014.00159", "nk", "cd14 mono", "Myeloid–NK cytokine crosstalk (myeloid IL-12/18 → NK; NK IFN-γ/TNF-α → myeloid)", 
        "10.3389/fimmu.2021.739220", "nk", "macrophage", "Myeloid–NK cytokine crosstalk (myeloid IL-12/18 → NK; NK IFN-γ/TNF-α → myeloid)", 
        "10.3389/fimmu.2021.739220", "nk", "pdc", "pDC type I IFN promotes NK activation", 
        "10.1126/science.284.5421.1835", "mait", "cdc", "MAIT cells can mature DCs via CD40L/GM-CSF (feedback)", 
        "10.4049/jimmunol.1700615", "mait", "monocytic", "Cytokine-dependent activation by monocytes; MR1-dependent with macrophages", 
        "10.3389/fimmu.2020.01014", "nkt", "cdc", "iNKT rapidly license DCs and shape adaptive responses", 
        "10.1016/j.coi.2007.03.014", "nkt", "b naive", "iNKT provide cognate B-cell help (iNKTfh; IL-21-dependent)", 
        "10.1038/ni.2172", "nkt", "b memory", "iNKT provide cognate B-cell help (iNKTfh; IL-21-dependent)", 
        "10.1038/ni.2166", "nkt", "plasma", "iNKT help promotes antibody-secreting cell generation", 
        "10.1038/ni.2166")
    nonimmune <- c("endocrine", "endothelial", "epithelial", 
        "fat", "glial", "muscle", "myoepithelial", "neuron", 
        "pericyte", "pneumocyte", "progenitor", "renal", "secretory", 
        "stromal", "reproductive", "mesothelial", "lens", "sensory", 
        "epidermal", "cartilage", "liver", "bone")
    immune <- setdiff(all_cell_types, nonimmune)
    nonimmune_immune_pairs <- bind_rows(tidyr::crossing(source = factor(nonimmune, 
        levels = nonimmune), target = factor(immune, levels = immune)), 
        tidyr::crossing(target = factor(nonimmune, levels = nonimmune), 
            source = factor(immune, levels = immune)))
    pairs_to_consider = bind_rows(myeloid_lymphoid_pairs, nonimmune_immune_pairs)
    distinct(bind_rows(select(pairs_to_consider, source, target), 
        set_names(select(pairs_to_consider, source, target), 
            c("target", "source"))))
}
list(tar_target(cellchat_file, "/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/cellNexus_lr_signaling_pathway_strength.duckdb", 
    format = "file"), tar_target(ethnicity_imputed, {
    check_rclone_installation()
    temp_path = tempdir()
    system(glue("~/bin/rclone copy box_adelaide:/Mangiola_ImmuneAtlas/reports/ning/data/All_pseudobulk_1_0_6_ethnicity_imputed_colData.csv {temp_path}/"))
    mutate(select(read_csv(glue("{temp_path}/All_pseudobulk_1_0_6_ethnicity_imputed_colData.csv")), 
        sample_id, ethnicity_groups, ethnicity_groups_imputed = finalEthnicity_groups), 
        ethnicity_groups_imputed = str_replace(ethnicity_groups_imputed, 
            "_imp$", "_imputed"))
}, packages = c("glue", "readr", "dplyr", "stringr")), tar_target(cell_to_exclude, 
    c("immune", "blood", "t", "erythrocyte", "other", "cytotoxic", 
        "b", "monocytic", "dc", "t cd8", "t cd4")), tar_target(cell_pathway_combination, 
    {
        comb = as_tibble(distinct(filter(tbl(dbConnect(duckdb::duckdb(), 
            dbdir = cellchat_file), "lr_pathway_table"), !source %in% 
            cell_to_exclude & !target %in% cell_to_exclude), 
            source, target, pathway_name))
        tar_group(group_by(mutate(arrange(group_by(inner_join(comb, 
            distinct(select(get_pairs_to_consider(unique(pull(pivot_longer(select(comb, 
                source, target), everything()), value))), source, 
                target))), source, target), source, target), 
            .chunk = ceiling(row_number()/10)), source, target, 
            .chunk))
    }, iteration = "group", packages = c("dbplyr", "duckdb", 
        "cellNexus", "tarchetypes"), resources = tar_resources(crew = tar_resources_crew("elastic_5"))), 
    tar_target(estimates_chunk, mutate(cell_pathway_combination, 
        brms_fit = pmap(list(source, target, pathway_name), function(s, 
            t, p) {
            con = dbConnect(duckdb::duckdb(), dbdir = cellchat_file, 
                read_only = TRUE)
            data = as_tibble(filter(tbl(con, "lr_pathway_table"), 
                source == s, target == t, pathway_name == p))
            data = droplevels(as_tibble(filter(prepare_database(data, 
                ethnicity_imputed), sex != "unknown")))
            dbDisconnect(con, shutdown = TRUE)
            if (nrow(distinct(data, sample_id)) < 100) return(NULL)
            colnames(data) = stringr::str_replace_all(colnames(data), 
                "_+", "_")
            formula_chr = "log10(pathway_prob) ~ 1 + age_bin*sex + disease_groups_altered + ethnicity_groups_imputed + assay_groups_altered +\n              (1 | dataset_id_altered) +\n              (1 + age_bin*sex  | tissue_groups)"
            if (nrow(distinct(data, disease_groups_altered)) == 
                1) formula_chr = str_remove_all(formula_chr, 
                fixed("+ disease_groups_altered"))
            if (nrow(distinct(data, ethnicity_groups_imputed)) == 
                1) formula_chr = str_remove_all(formula_chr, 
                fixed("+ ethnicity_groups_imputed"))
            if (nrow(distinct(data, assay_groups_altered)) == 
                1) formula_chr = str_remove_all(formula_chr, 
                fixed("+ assay_groups_altered"))
            if (nrow(distinct(data, dataset_id_altered)) == 1) formula_chr = str_remove_all(formula_chr, 
                fixed("(1 | dataset_id_altered) +"))
            if (nrow(distinct(data, sex)) == 1) formula_chr = str_remove_all(formula_chr, 
                fixed("*sex"))
            formula <- bf(as.formula(formula_chr), sigma ~ 1)
            prior = eval(substitute(c(prior(student_t(3, i, s), 
                class = Intercept), prior(student_t(3, 0, s), 
                class = Intercept, dpar = "sigma"), prior(student_t(3, 
                0, 0.5), class = b)), env = list(i = mean(sqrt(data$interaction_weight)), 
                s = sd(sqrt(data$interaction_weight) * 2))))
            chains = 2
            brm(formula = formula, data = data, family = gaussian(), 
                prior = prior, chains = chains, cores = pmin(as.numeric(parallelly::availableCores()), 
                  chains), threads = threading(threads = floor((as.numeric(parallelly::availableCores())/chains))), 
                warmup = 500, refresh = 10, backend = "cmdstanr", 
                iter = 2400)
        }, .progress = TRUE)), pattern = map(cell_pathway_combination), 
        packages = c("brms", "glue", "stringr", "dplyr", "purrr"), 
        resources = tar_resources(crew = tar_resources_crew("elastic_5"))), 
    tar_target(hypothesis_age_monotonic, select(mutate(mutate(estimates_chunk, 
        hypothesis_age_monotonic = map(brms_fit, fit_to_age_monotonic_changes)), 
        summary = map(brms_fit, ~{
            if (is.null(.x)) return(NULL)
            rename(posterior::summarise_draws(.x), parameter = variable)
        })), -brms_fit), pattern = map(estimates_chunk), packages = c("brms", 
        "glue", "dplyr", "purrr", "rstan", "tibble", "purrr", 
        "posterior", "stringr"), resources = tar_resources(crew = tar_resources_crew("elastic_5"))))
