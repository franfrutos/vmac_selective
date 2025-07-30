# Set-up ----
if (!require(pacman)) {
  install.packages(pacman)
  library(pacman)
}

p_load(here, dplyr, readr, lme4, lmerTest, hypr, simr, pbapply, future, parallel, binom, ggplot2, marginaleffects, metafor, purrr, tidyr, patchwork)

source(here::here("Scripts/functions.R")) # Load data and all relevant functions

# Loading raw data from Garre-Frutos et al. (2024, 2025) from github repo
dlist <- list(
  "expM" = read.csv("https://raw.githubusercontent.com/franfrutos/VMAC_awareness/refs/heads/main/Input/data/raw_data_multi.csv")[,-1],
  "exp1" = read.csv("https://raw.githubusercontent.com/franfrutos/VMAC_awareness/refs/heads/main/Input/data/experiment_1.csv")[,-1],
  "exp2" = read.csv("https://raw.githubusercontent.com/franfrutos/VMAC_awareness/refs/heads/main/Input/data/experiment_2.csv")[,-1]
)

# Filter data from experiments in the original studies
dlist[["expM"]] <- filter_data(dlist[["expM"]], f.absent = F,
                               sd_filter = 2,
                               experiment = "eM") %>% select(Singleton, rt, ID, Block) %>%
  mutate(ID = paste0(ID, "expM"),
         group = "Instructions",
         Exp = "expM")

dlist[["exp1"]] <- filter_data(dlist[["exp1"]], f.absent = F,
                               sd_filter = 2,
                               experiment = "e1") %>% select(Singleton, rt, ID, Block) %>%
  mutate(ID = paste0(ID, "exp1"),
         group = "No instructions",
         Exp = "exp1")
dlist[["exp2"]] <- filter_data(dlist[["exp2"]], f.absent = F,
                               sd_filter = 2,
                               experiment = "e2") %>% select(Singleton, rt, ID, group, Block_num) %>%
  mutate(ID = paste0(ID, "exp2"),
         group = ifelse(group == "A", "Instructions", "No instructions"),
         Exp = "exp2",
         Block = Block_num
  )

# Loading data from current study: 

# Applying exclusion criteria (see: https://osf.io/f3bm8):
# Excluding participants with numerically low accuracy in the visual search task (< 0.7)
raw <- read.csv("Input/data.csv")

acc_ex <- raw %>%
  filter(Phase == "Rewarded") %>%
  summarise(ACC = mean(correct), .by = c(task, ID))%>%
  filter(ACC < .7) %>%
  pull(ID)

# Excluding participants not performing above chance in the report task
report_ex <-
  raw %>%
  filter(Phase %in% c("Rewarded", "Report")) %>%
  group_by(ID) %>%
  mutate(tmp = 1:n()) %>%
  ungroup() %>%
  filter(tmp > 30) %>%
  select(-tmp) %>%
  filter(Phase == "Report") %>%
  summarise(correct = sum(correct), n = n(), .by = c(task, ID))%>%
  mutate(
    ACC = correct / n,
    p_value = map2_dbl(correct, n, ~ binom.test(.x, .y)$p.value),
    cut_off = mean(ACC)-sd(ACC)*2.5
  ) %>%
  filter(p_value > .05 | cut_off > ACC)%>%
  pull(ID)

# And participants whose mean RTs could be considered outliers (+- 3SDs)
rt_ex <-
  raw %>%
  filter(Phase == "Rewarded", correct == 1, rt > 150, rt < 1800) %>%
  summarise(rt = mean(rt), .by = c(task, ID))%>%
  mutate(
    cut_off_low = mean(rt)-sd(rt)*3,
    cut_off_high = mean(rt)+sd(rt)*3
  ) %>%
  filter(cut_off_high < rt | cut_off_low > rt)%>%
  pull(ID)

# Final vector of all excluded participants based on above criteria
exclusions <- unique(c(report_ex, acc_ex, rt_ex))

# Filter data
d <- raw %>%
  filter(Phase == "Rewarded", rt > 150, rt < 1800,
         correct == 1, !ID %in% c(exclusions), Block_num > 2)

# Calculate effect sizes ----
effs <- rbind(
  dlist[["expM"]] %>%
    summarise(RT = mean(log(rt)), .by = c(ID, Singleton)) %>%
    filter(Singleton != "Absent") %>%
    pivot_wider(names_from = "Singleton", values_from = "RT") %>%
    summarise(
      mean_High = mean(High),
      mean_Low  = mean(Low),
      sd_High   = sd(High),
      sd_Low    = sd(Low),
      cor       = cor(High, Low),
      n         = n()
    ) %>%
    mutate(Study = "Garre-Frutos et al. (2024) (Instructions)",
           manipulation = 1),
  
  dlist[["exp1"]] %>%
    summarise(RT = mean(log(rt)), .by = c(ID, Singleton)) %>%
    filter(Singleton != "Absent") %>%
    pivot_wider(names_from = "Singleton", values_from = "RT") %>%
    summarise(
      mean_High = mean(High),
      mean_Low  = mean(Low),
      sd_High   = sd(High),
      sd_Low    = sd(Low),
      cor       = cor(High, Low),
      n         = n()
    ) %>%
    mutate(Study =   "Garre-Frutos, Lupiáñez et al. (2025) - Experiment 1 (No instructions)",
           manipulation = 0),
  
  dlist[["exp2"]] %>%
    summarise(RT = mean(log(rt)), .by = c(group, ID, Singleton)) %>%
    filter(Singleton != "Absent") %>%
    pivot_wider(names_from = "Singleton", values_from = "RT") %>%
    summarise(
      mean_High = mean(High),
      mean_Low  = mean(Low),
      sd_High   = sd(High),
      sd_Low    = sd(Low),
      cor       = cor(High, Low),
      n         = n(),
      .by      = group
    ) %>%
    mutate(
      Study = c(
        "Garre-Frutos, Lupiáñez et al. (2025) - Experiment 2 (Instructions)",
        "Garre-Frutos, Lupiáñez et al. (2025) - Experiment 2 (No instructions)"),
        manipulation = c(1, 0)
      ) %>%
    select(-group),
  
  d %>%
    filter(
      Phase == "Rewarded",
      rt > 150, rt < 1800,
      correct == 1,
      !ID %in% exclusions
    ) %>%
    summarise(RT = mean(log(rt)), .by = c(task, ID, Singleton)) %>%
    filter(Singleton != "Absent") %>%
    pivot_wider(names_from = "Singleton", values_from = "RT") %>%
    summarise(
      mean_High = mean(High),
      mean_Low  = mean(Low),
      sd_High   = sd(High),
      sd_Low    = sd(Low),
      cor       = cor(High, Low),
      n         = n(),
      .by      = task
    ) %>%
    mutate(
      Study = c("Present study - Location task", "Present study - Color task"),
      manipulation = c(0, 1)
    ) %>%
    select(-task)
)

SMCC <- escalc(
  measure = "SMCC",
  m1      = mean_High,
  m2      = mean_Low,
  sd1     = sd_High,
  sd2     = sd_Low,
  ni      = n,
  ri      = cor,
  data    = effs,
  slab    = Study
)

# Fit model ----
metaFit <- rma(yi = yi, vi = vi, data = SMCC, method = "REML",
               mods = ~manipulation)

# Prediction for each type of study
preds <- predict(metaFit, newmods = matrix(c(1, 0)))

# Plot meta-analysis ----

# format data
levels_chrono <- c(
  "Garre-Frutos et al. (2024) (Instructions)",
  "Garre-Frutos, Lupiáñez et al. (2025) - Experiment 1 (No instructions)",
  "Garre-Frutos, Lupiáñez et al. (2025) - Experiment 2 (Instructions)",
  "Garre-Frutos, Lupiáñez et al. (2025) - Experiment 2 (No instructions)",
  "Present study - Location task",
  "Present study - Color task"
)

df_studies <- tibble(
  study  = SMCC$Study,
  est    = SMCC$yi,
  vi     = SMCC$vi,
  weight = 1/SMCC$vi
) %>%
  mutate(
    lo    = est - qnorm(0.975)*sqrt(vi),
    hi    = est + qnorm(0.975)*sqrt(vi),
    label = sprintf("%.2f [%.2f, %.2f]", est, lo, hi),
    group = case_when(
      grepl("Instructions", study) | grepl("Color task", study) ~ "Instr",
      TRUE                                                     ~ "NoInstr"
    ),
    study = factor(study, levels = rev(levels_chrono))
  )

df_summary <- tibble(
  study  = "Meta-estimate",
  est    = preds$pred,
  lo     = preds$ci.lb,
  hi     = preds$ci.ub,
  group = c("Instr", "NoInstr"),
  weight = NA_real_,
  label  = sprintf("%.2f [%.2f, %.2f]", preds$pred, preds$ci.lb, preds$ci.ub)
)

# color palette
cols  <- c(Instr = "#258da5", NoInstr = "#C15F1E")
fills <- c(Instr = "#7CBAC9", NoInstr = "#E0A67A")

# plot
# individual studies plot
p1 <- ggplot(df_studies, aes(y = study, x = est)) +
  geom_errorbarh(aes(xmin = lo, xmax = hi, color = group),
                 height = 0, size = 0.8) +
  geom_point(aes(size = weight, fill = group, color = group),
             shape = 22) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values  = fills) +
  scale_size_area(max_size = 6) +
  geom_text(aes(x = .6, label = label),
            hjust = 0, nudge_y = 0.15, size = 4) +
  scale_x_continuous(name   = "Standardized Mean Change",
                     limits = c(-0.2,0.8),
                     breaks = seq(-0.2,1,0.2)) +
  geom_vline(xintercept = 0,    linetype = "dashed", color = "gray70") +
  geom_vline(xintercept = 0.35, linetype = "dashed", color = cols["Instr"], size = 1) +
  coord_cartesian(xlim = c(-0.2, 0.8)) +
  labs(y = "Study") +
  theme_Publication(base_size = 12, text_size = 12) +
  theme(
    legend.position     = "none",
    axis.title.x        = element_blank(),
    axis.text.x         = element_blank(),
    axis.ticks.x        = element_blank(),
    axis.line.x         = element_blank(),
    panel.grid.major.y  = element_blank(),
    plot.margin         = margin(t = 5, r = 10, b = 0, l = 10)
  )

# p2: Summary
p2 <- ggplot(df_summary, aes(y = study, x = est, color = group, fill = group)) +
  geom_errorbarh(aes(xmin = lo, xmax = hi),
                 height = 0, size = 0.8) +
  geom_point(shape = 21, size = 6) +
  geom_text(aes(label = label),
            hjust = 0, nudge_y = 0.2, size = 4) +
  geom_vline(xintercept = 0,    linetype = "dashed", color = "gray70") +
  geom_vline(xintercept = 0.35, linetype = "dashed", color = cols["Instr"], size = 1) +
  scale_x_continuous(name   = "Standardized Mean Difference",
                     limits = c(-0.2,0.8),
                     breaks = seq(-0.2,1,0.2)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values  = fills) +
  theme_Publication(base_size = 12, text_size = 12) +
  theme(
    axis.title.y        = element_blank(),
    legend.position     = "none",
    panel.grid.major.y  = element_blank(),
    plot.margin         = margin(t = 0, r = 10, b = 5, l = 10)
  )

# save plot ----
final_plot <- p1 / p2 + plot_layout(heights = c(5, 1))

# Save figure with customs dimensions:
if(!file.exists("Output/plots/fig3.pdf")) {
  ggsave(plot=final_plot, device = "pdf", filename=here::here("Output/plots/fig3.pdf"),
         units = "cm", height = 17, width = 27)
}

if(!file.exists("Output/plots/fig3.png")) {
  ggsave(plot=final_plot, filename=here::here("Output/plots/fig3.png"),
         units = "cm", height = 17, width = 27)
}

# Interaction: ----
# 
effs2 <- rbind(
  dlist[["exp2"]] %>%
    summarise(RT = mean(log(rt)), .by = c(group, ID, Singleton)) %>%
    filter(Singleton != "Absent") %>%
    mutate(group = ifelse(group == "Instructions", 1, 0)) %>%
    pivot_wider(names_from = "Singleton", values_from = "RT") %>%
    summarise(
      VMAC = mean(High) - mean(Low),
      .by      = c(ID, group)
    ) %>%
    summarise(
      mean_VMAC = mean(VMAC),
      sd_VMAC = sd(VMAC),
      n = n(),
      .by = group
    ) %>%
    pivot_wider(names_from = "group", values_from = c("mean_VMAC", "sd_VMAC", "n")) %>%
    mutate(
      Study = c("Garre-Frutos, Lupiáñez et al. (2025) - Experiment 2"),
    ),
  
  d %>%
    filter(
      Phase == "Rewarded",
      rt > 150, rt < 1800,
      correct == 1,
      !ID %in% exclusions
    ) %>%
    summarise(RT = mean(log(rt)), .by = c(task, ID, Singleton)) %>%
    filter(Singleton != "Absent") %>%
    mutate(group = ifelse(task == "C", 1, 0)) %>%
    pivot_wider(names_from = "Singleton", values_from = "RT") %>%
    summarise(
      VMAC = mean(High) - mean(Low),
      .by      = c(ID, group)
    ) %>%
    summarise(
      mean_VMAC = mean(VMAC),
      sd_VMAC = sd(VMAC),
      n = n(),
      .by = group
    ) %>%
    pivot_wider(names_from = "group", values_from = c("mean_VMAC", "sd_VMAC", "n")) %>%
    mutate(
      Study = c("Present study"),
    )
)

SMD <- escalc(
  measure = "SMD",
  m1      = mean_VMAC_1,
  m2      = mean_VMAC_0,
  sd1     = sd_VMAC_1,
  sd2     = sd_VMAC_0,
  n1i     = n_1,
  n2i     = n_0,
  data    = effs2,
  slab    = Study
)

meta2 <- rma(yi = yi, vi = vi, method = "REML", data = SMD)

# Construct df
df_studies <- tibble(
  study  = SMD$Study,
  est    = SMD$yi,
  vi     = SMD$vi,
  weight = 1/SMD$vi
) %>%
  mutate(
    lo    = est - qnorm(0.975)*sqrt(vi),
    hi    = est + qnorm(0.975)*sqrt(vi),
    label = sprintf("%.2f [%.2f, %.2f]", est, lo, hi),
    study = factor(study, levels = c(
      c(
        "Garre-Frutos, Lupiáñez et al. (2025) - Experiment 2",
        "Present study"
      )
    ))
  )

df_summary <- tibble(
  study  = "Meta-estimate",
  est    = meta2$b[1,1],
  lo     = meta2$ci.lb,
  hi     = meta2$ci.ub,
  weight = NA_real_,
  label  = sprintf("%.2f [%.2f, %.2f]", est, lo, hi)
)

# Plot
p1 <- ggplot(df_studies, aes(y = study, x = est)) +
  geom_errorbarh(
    data = df_studies[!is.infinite(df_studies$hi),],
    aes(xmin = lo, xmax = hi), height = 0, size = 0.8,
    color = "#6C3483") +
  geom_segment(data = filter(df_studies, is.infinite(hi)),
               aes(x = lo, xend = est, y = study, yend = study),
               size = 0.8, color = "#6C3483") +
  geom_segment(data = filter(df_studies, is.infinite(hi)),
               aes(x = est, y = study, xend = 0.8, yend = study),
               arrow = arrow(type="closed", length=unit(0.15,"cm")),
               size = 0.8, color = "#6C3483") +
  geom_point(shape = 22, fill = "#D7BDE2", color = "#6C3483", aes(size = weight)) +
  geom_text(aes(x = 0.8, label = label), hjust = 1, size = 4, nudge_y = 0.15) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
  geom_vline(xintercept = 0.349, linetype = "dashed", color = "#6C3483", size = 1) +
  scale_x_continuous(limits = c(0,0.8), breaks = seq(-0.2,1,0.2)) +
  theme_Publication(base_size = 12, text_size = 12) +
  labs(y = "Study") +
  scale_size(
    range = c(4, 5)    # tamaño máximo que quieras
  ) +
  theme(
    legend.position = "none",
    axis.title.x     = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.line.x      = element_blank(),  # <— quita la barrita Y
    panel.grid.major.y = element_blank(),
    plot.margin      = margin(t = 5, r = 10, b = 0, l = 10)
  )

p2 <- ggplot(df_summary, aes(y = study, x = est)) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, size = 0.8, color = "#6C3483") +
  geom_point(shape = 21, size = 6, fill = "#D7BDE2", color = "#6C3483") +
  geom_vline(xintercept = 0.349, linetype = "dashed", color = "#6C3483", size = 1) +
  geom_text(aes(x = 0.8, label = label), hjust = 1, size = 4, nudge_y = 0.1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
  scale_x_continuous(name = "Standarize Mean Difference",
                     limits = c(0,0.8),
                     breaks = seq(-0.2,1,0.2)) +
  theme_Publication(base_size = 12, text_size = 12) +
  theme(
    axis.title.y      = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.margin      = margin(t = 0, r = 10, b = 5, l = 10)
  )

(final_plot <- p1 / p2 + plot_layout(heights = c(5,1)))

# Save figure with customs dimensions:
if(!file.exists("Output/plots/fig4.pdf")) {
  ggsave(plot=final_plot, device = "pdf", filename=here::here("Output/plots/fig4.pdf"),
         units = "cm", height = 8, width = 27)
}

if(!file.exists("Output/plots/fig4.png")) {
  ggsave(plot=final_plot, filename=here::here("Output/plots/fig4.png"),
         units = "cm", height = 8, width = 27)
}
