# VEGFA and PCNA expression in breast and lung cancer cell lines (DepMap 25Q2)
# Run from the project folder. Needs the three files below in ./data

library(tidyverse)

expr_file    <- "data/OmicsExpressionProteinCodingGenesTPMLogp1Stranded.csv"
model_file   <- "data/Model.csv"
hotspot_file <- "data/OmicsSomaticMutationsMatrixHotspot.csv"

my_lines <- c("MDAMB231", "MCF7", "NCIH358")

dir.create("figures", showWarnings = FALSE)

# 1. Read only the genes we need from the large expression file ----------------
get_genes <- function(file, genes) {
  hdr <- names(read_csv(file, n_max = 0, show_col_types = FALSE))
  pattern <- paste0("^(", paste(genes, collapse = "|"), ") \\(")
  idx <- c(1, grep(pattern, hdr))
  read_csv(file, col_select = all_of(idx), show_col_types = FALSE) %>%
    rename(ModelID = 1) %>%
    rename_with(~ sub(" \\(.*", "", .x), .cols = -1)
}

genes <- get_genes(expr_file, c("VEGFA", "PCNA", "ESR1", "PGR", "ERBB2"))

# 2. Join to cell line annotation, keep breast and lung -------------------------
model <- read_csv(model_file, show_col_types = FALSE)

d <- genes %>%
  left_join(model, by = "ModelID") %>%
  filter(OncotreeLineage %in% c("Breast", "Lung")) %>%
  mutate(mine = StrippedCellLineName %in% my_lines)

# 3. Plot helper ---------------------------------------------------------------
plot_genes <- function(df, group, title) {
  long <- df %>%
    pivot_longer(c(VEGFA, PCNA), names_to = "gene", values_to = "expression")
  ggplot(long, aes(.data[[group]], expression)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
    geom_point(data = filter(long, mine), colour = "red", size = 3) +
    geom_text(data = filter(long, mine), aes(label = StrippedCellLineName),
              colour = "red", hjust = -0.2, size = 3) +
    facet_wrap(~gene) +
    labs(x = NULL, y = "Expression (log2 TPM + 1)", title = title) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
}

# 4. Figure 1: breast vs lung --------------------------------------------------
p1 <- plot_genes(d, "OncotreeLineage", "Breast and lung cancer cell lines")
ggsave("figures/expression_breast_lung.png", p1, width = 8, height = 5, dpi = 300)

# 5. Figure 2: NSCLC by KRAS hotspot status ------------------------------------
hot <- read_csv(hotspot_file, show_col_types = FALSE)
kras_col <- grep("^KRAS \\(", names(hot), value = TRUE)
kras <- hot %>%
  select(ModelID = 1, all_of(kras_col)) %>%
  rename(KRAS_hotspot = 2)

nsclc <- d %>%
  filter(OncotreePrimaryDisease == "Non-Small Cell Lung Cancer") %>%
  left_join(kras, by = "ModelID") %>%
  filter(!is.na(KRAS_hotspot)) %>%
  mutate(KRAS = ifelse(KRAS_hotspot > 0, "KRAS hotspot mutant", "No KRAS hotspot"))

print(count(nsclc, KRAS))

p2 <- plot_genes(nsclc, "KRAS", "NSCLC cell lines")
ggsave("figures/expression_nsclc_kras.png", p2, width = 8, height = 5, dpi = 300)

# 6. Figure 3: breast by expression-based subtype ------------------------------
# DepMap does not label triple-negative breast cancer, so lines are grouped by
# ESR1 and ERBB2 expression. Cutoffs are rough and only checked on MDA-MB-231
# (should be TNBC-like) and MCF-7 (should be ER-positive-like).
breast <- d %>%
  filter(OncotreeLineage == "Breast") %>%
  mutate(type = case_when(
    ESR1 < 2 & ERBB2 < 7 ~ "TNBC-like",
    ESR1 >= 2            ~ "ER-positive-like",
    TRUE                 ~ "HER2-high, ER-low"))

print(count(breast, type))

p3 <- plot_genes(breast, "type", "Breast cancer cell lines")
ggsave("figures/expression_breast_subtype.png", p3, width = 8, height = 5, dpi = 300)

# 7. Results table for the three cell lines ------------------------------------
results <- bind_rows(
  breast %>% filter(mine) %>% mutate(group = type),
  nsclc  %>% filter(mine) %>% mutate(group = KRAS)
) %>%
  select(cell_line = StrippedCellLineName, group, VEGFA, PCNA)

print(results)
write_csv(results, "results_my_lines.csv")

# 8. Record the R session for reproducibility ----------------------------------
writeLines(capture.output(sessionInfo()), "session_info.txt")
