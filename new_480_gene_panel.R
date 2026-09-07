library(dplyr)   
library(readxl)  
library(purrr)      

PANEL_SIZE <- 480   # final panel size

nn_sheets <- c("HA ECM", "Calcification Mineralisation",   # sheets to pull non-negotiables
               "Inflammation Uraemia Fibrosis", "Therapeutic Targets")

non_negotiable <- map_dfr(nn_sheets, function(sh) {        # load + clean non-negotiables
  read_excel("SLM_shorter.xlsx", sheet = sh) %>%
    transmute(symbol = toupper(trimws(Gene)), non_negotiable_sheet = sh)
}) %>%
  distinct(symbol, .keep_all = TRUE)

consensus_wide <- read.csv("/Users/tyleradams/Desktop/DISS/consensus_panel_480.csv")   # load consensus
consensus_wide <- consensus_wide %>%
  mutate(symbol_upper = toupper(trimws(symbol)))            # uppercase key

nn_in_model <- consensus_wide %>%                           # non-negotiables present in model
  filter(symbol_upper %in% non_negotiable$symbol) %>%
  mutate(inclusion = "non_negotiable")

nn_missing <- non_negotiable %>%                            # non-negotiables missing from model
  filter(!symbol %in% consensus_wide$symbol_upper) %>%
  transmute(symbol, symbol_upper = symbol, inclusion = "non_negotiable")

n_forced <- nrow(nn_in_model) + nrow(nn_missing)            # forced count
n_fill   <- PANEL_SIZE - n_forced                           # remaining slots

if (n_fill < 0) {                                            # sanity check
  stop(sprintf("Non-negotiable genes alone (%d) exceed PANEL_SIZE (%d) - can't build a %d-gene list.",
               n_forced, PANEL_SIZE, PANEL_SIZE))
}

fill_genes <- consensus_wide %>%                            # ranked fill genes
  filter(!symbol_upper %in% non_negotiable$symbol) %>%
  slice(1:n_fill) %>%
  mutate(inclusion = "ranked_fill")

new_480_panel <- bind_rows(nn_in_model, nn_missing, fill_genes) %>%   # combine
  left_join(non_negotiable, by = "symbol") %>%
  distinct(symbol_upper, .keep_all = TRUE) %>%
  select(-symbol_upper)

write.csv(new_480_panel, "new_480_gene_panel.csv", row.names = FALSE)   # save panel

cat("Panel size:", nrow(new_480_panel), "(should equal", PANEL_SIZE, ")\n")
cat("  forced (non-negotiable):", n_forced, "\n")
cat("  ranked fill:", nrow(fill_genes), "\n")

consensus_full <- read.csv("/Users/tyleradams/Desktop/DISS/consensus_gene_table.csv",   # full consensus
                           stringsAsFactors = FALSE) %>%
  mutate(symbol_upper = toupper(trimws(symbol)))

wkru_consensus_overlap <- new_480_panel %>%                # WKRU + consensus membership
  mutate(in_consensus_ranking = symbol %in% consensus_full$symbol_upper) %>%
  select(symbol, inclusion, non_negotiable_sheet, in_consensus_ranking)

write.csv(wkru_consensus_overlap,
          "/Users/tyleradams/Desktop/DISS/WKRU_consensus_overlap.csv",
          row.names = FALSE)
