# ─────────────────────────────────────────────────────────────────────────────
# Load Required Packages
# ─────────────────────────────────────────────────────────────────────────────
req_pkgs <- c(
  "readxl", "dplyr", "ggplot2", "janitor",
  "lme4", "lmerTest", "broom.mixed", "tidyr", "stringr",
  "performance", 
)
new_pkgs <- setdiff(req_pkgs, rownames(installed.packages()))
if (length(new_pkgs)) install.packages(new_pkgs, dependencies = TRUE)
invisible(lapply(req_pkgs, library, character.only = TRUE))

# ─────────────────────────────────────────────────────────────────────────────
# Utility: Generate Fixed Effect Labels from Model
# ─────────────────────────────────────────────────────────────────────────────
make_fe_label <- function(mod, digits = 2) {
  tab <- broom.mixed::tidy(mod, effects = "fixed")
  
  if (!"p.value" %in% names(tab)) {
    mod <- lmerTest::lmer(formula(mod), data = model.frame(mod), REML = FALSE)
    tab <- broom.mixed::tidy(mod, effects = "fixed")
  }
  
  tab %>%
    mutate(
      estimate = signif(estimate, digits),
      p.value = signif(p.value, digits),
      term = dplyr::recode(term, `(Intercept)` = "Intercept")
    ) %>%
    transmute(lbl = sprintf("%s = %s (p = %s)", term, estimate, p.value)) %>%
    dplyr::pull(lbl) %>%
    paste(collapse = "\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# Utility: Convert p-values to Significance Stars
# ─────────────────────────────────────────────────────────────────────────────
pval_to_stars <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# Plot Configuration Parameters
# ─────────────────────────────────────────────────────────────────────────────
bottom_expand_mult <- 0.20  # padding below violin plot
upper_expand_mult  <- 0.25  # padding above violin plot

# ─────────────────────────────────────────────────────────────────────────────
# Load and Preprocess Dataset
# ─────────────────────────────────────────────────────────────────────────────
file_path <- "/Users/dhanyata.n/Documents/Stats tests/Round 12 - Reformatted R/Sound experiment results Reformatted.xlsx"
df <- readxl::read_excel(file_path) %>% janitor::clean_names()

# Convert pixels to microns (µm)
px_to_um <- 1 / 0.7752
len_vars  <- intersect(names(df), c("perimeter"))
area_vars <- intersect(names(df), c("area"))
if (length(len_vars))  df[len_vars]  <- df[len_vars]  * px_to_um
if (length(area_vars)) df[area_vars] <- df[area_vars] * (px_to_um^2)

# ─────────────────────────────────────────────────────────────────────────────
# Recode Identifiers and Set Factor Levels
# ─────────────────────────────────────────────────────────────────────────────
df <- df %>%
  mutate(
    timepoint = factor(timepoint, levels = c(
      "Pretreatment", "Exposure 1 Immediate Post",
      "Exposure 1 Delayed Post", "Exposure 2 Immediate Post",
      "Exposure 2 Delayed Post")),
    group               = factor(group, levels = c("Control", "Sound")),
    patient             = factor(patient),
    technical_replicate = factor(technical_replicate)
  )

# ─────────────────────────────────────────────────────────────────────────────
# Summary Table: Row Counts by Group and Timepoint
# ─────────────────────────────────────────────────────────────────────────────
cat("\nRow counts by group × timepoint:\n")
print(
  df %>%
    count(group, timepoint) %>%
    tidyr::pivot_wider(names_from = group, values_from = n, values_fill = 0)
)

# ─────────────────────────────────────────────────────────────────────────────
# Fit Mixed-Effects Models for All Numeric Features
# ─────────────────────────────────────────────────────────────────────────────
output_dir <- "~/Documents/Stats tests/Round 12 - Reformatted R/model_diagnostics"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

all_numeric <- df %>% select(where(is.numeric)) %>% names()
num_feats   <- setdiff(all_numeric, c("organoid_id", "technical_replicate"))
model_list  <- list()

# Add log1p-transformed columns for all numeric features
df <- df %>%
  mutate(across(all_of(num_feats), ~ log1p(.x), .names = "{.col}_log"))
num_feats_log <- paste0(num_feats, "_log")

# Fit models and save diagnostics
for (feat in num_feats_log) {
  dat <- df %>% drop_na(all_of(feat)) %>% droplevels()
  
  if (any(vapply(dat[c("group", "timepoint")],
                 \(x) nlevels(x) < 2, logical(1)))) {
    message("Skipping ", feat, ": <2 levels in a predictor.")
    next
  }
  
  mod <- lmerTest::lmer(
    paste0(feat, " ~ group * timepoint + (1|patient/technical_replicate)"),
    data = dat,
    REML = FALSE
  )
  
  cat("\n==== ", feat, " (log1p) ====\n")
  print(broom.mixed::tidy(mod, effects = "fixed"))
  
  diag_plot <- plot(performance::check_model(mod))
  ggsave(
    filename = file.path(output_dir, paste0(feat, "_model_diagnostics.png")),
    plot = diag_plot, width = 8, height = 6, dpi = 300)
  
  orig_feat <- sub("_log$", "", feat)
  model_list[[orig_feat]] <- mod
}

# ─────────────────────────────────────────────────────────────────────────────
# Generate Violin Plots for Each Numeric Feature with Significance Brackets
# ─────────────────────────────────────────────────────────────────────────────
metric_units <- c(
  area = "µm²", perimeter = "µm",
  solidity = "", roundness = "", aspect_ratio = "", circularity = ""
)

for (feat in num_feats) {
  this_mod <- model_list[[feat]]
  unit_lbl <- metric_units[[feat]] %||% ""
  clean_feat <- str_to_title(str_replace_all(feat, "_", " "))
  axis_lbl <- if (unit_lbl == "") clean_feat else paste0(clean_feat, " (", unit_lbl, ")")
  
  p <- ggplot(df, aes(group, !!sym(feat), fill = group)) +
    geom_violin(trim = FALSE, colour = NA, alpha = 0.7, scale = "width", adjust = 2) +
    geom_boxplot(width = .12, fill = "white", outlier.shape = NA) +
    facet_wrap(~ timepoint, nrow = 1, scales = "free_y", drop = FALSE) +
    labs(y = axis_lbl, x = "Group",
         title = paste("Distribution of", axis_lbl), fill = "Group") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(size = 9)) +
    scale_y_continuous(expand = expansion(mult = c(bottom_expand_mult, upper_expand_mult)))
  
  # Annotate plots with significance stars
  star_df <- broom.mixed::tidy(this_mod, effects = "fixed") %>%
    filter(term == "groupSound" | str_detect(term, "^groupSound:timepoint")) %>%
    mutate(
      timepoint = if_else(term == "groupSound", "Pretreatment",
                          str_remove(term, "^groupSound:timepoint")),
      stars = pval_to_stars(p.value)
    ) %>%
    select(timepoint, stars)
  
  y_limits <- df %>%
    group_by(timepoint) %>%
    summarise(y = max(.data[[feat]], na.rm = TRUE) * 1.05, .groups = "drop")
  
  star_df <- left_join(star_df, y_limits, by = "timepoint") %>%
    mutate(x = 1.5)
  
  bracket_df <- star_df %>%
    transmute(x_start = 1, x_end = 2, y = y, timepoint = timepoint)
  
  p <- p +
    geom_segment(data = bracket_df, aes(x = x_start, xend = x_end, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.6) +
    geom_segment(data = bracket_df, aes(x = x_start, xend = x_start,
                                        y = y, yend = y - 0.02 * diff(range(df[[feat]], na.rm = TRUE))),
                 inherit.aes = FALSE, linewidth = 0.6) +
    geom_segment(data = bracket_df, aes(x = x_end, xend = x_end,
                                        y = y, yend = y - 0.02 * diff(range(df[[feat]], na.rm = TRUE))),
                 inherit.aes = FALSE, linewidth = 0.6) +
    geom_text(data = star_df, aes(x = x, y = y + 0.08 * diff(range(df[[feat]], na.rm = TRUE)),
                                  label = stars), inherit.aes = FALSE, size = 5)
  
  print(p)
  ggsave(file.path(output_dir, paste0(feat, "_distribution_wFE.png")),
         p, width = 8, height = 6, dpi = 300)
}

# ─────────────────────────────────────────────────────────────────────────────
# End of Script
# ─────────────────────────────────────────────────────────────────────────────
