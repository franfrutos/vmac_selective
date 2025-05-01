# Author: Francisco Garre-Frutos

# Loading packages ----
if (!require(pacman)) install.packages(pacman)

p_load(tidyr, Rmisc, here, dplyr, stringr,
       ggplot2, betareg, hypr, marginaleffects, lmerTest,
       purrr, ggpubr)

summarise <- dplyr::summarise
mutate <- dplyr::mutate

# Load raw data ----
raw <- read.csv("Input/data.csv")

# Applying exclusion criteria (see: https://osf.io/f3bm8):
# Excluding participants with numerically low accuracy in the visual search task (< 0.7)
acc_ex <-
  raw %>%
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

# Pre-registered analysis ----
# Fitted model list to use in reproducible manuscript
fit_list <- list()

# Visual search task ----
# RTs:
HcRep <- hypr(
  VMAC = High ~ Low,
  AC = Low ~ Absent,
  levels = c("High", "Low", "Absent")
) # Contrast for VMAC and AC effects

# Filter data
d <- raw %>%
  filter(Phase == "Rewarded", rt > 150, rt < 1800,
         correct == 1, !ID %in% c(exclusions), Block_num > 2)

# Set contrasts for model fitting
d$Singleton <- factor(d$Singleton, levels = c("High", "Low", "Absent"))
d$Task <- factor(d$task, levels = c("C", "L"))
contrasts(d$Task) <- contr.sum(2)
colnames(contrasts(d$Task)) <- c("c_vs_L")
d$log_RT <- log(d$rt)
contrasts(d$Singleton) <- contr.hypothesis(HcRep)

fit <-
  lmer(
    log_RT ~ Singleton * Task + (Singleton | ID),
    control = lmerControl(optimizer = 'bobyqa'),
    data = d
  )

# Model (raw) summary
summary(fit)

# one-tailed p-value for VMAC:Task
pt(summary(fit)$coefficients[5, 4], summary(fit)$coefficients[5, 3], lower.tail = F) 

# Save for later
fit_list[["RTs_VST"]] <- fit

# Reliability of VMAC effect:
rel <- multi.s::splith(
  outcome = "rt",
  data = d[d$Singleton != "Absent",],
  permutations = 5e3,
  variable = "Singleton",
  subject = "ID",
  include_block = F,
  block = "Epoch",
  average = "mean"
)$final_estimates

# Plot
dodge <- position_dodge(width = 0.4)

d_within <- summarySEwithin(
  data = d %>%
    summarise(rt = mean(rt), .by = c(task, ID, Singleton)) %>%
    mutate(Singleton = factor(Singleton, levels = c("High", "Low", "Absent"))),               # The data frame containing your data
  measurevar = "rt",            # The dependent variable (e.g., response time)
  withinvars = c("Singleton"),  # The within-subject factors (e.g., distractor type)
  betweenvars = c("task"),      # The between-subject factor (e.g., task type)
  idvar = "ID"                  # The subject identifier (e.g., participant ID)
)

fig2a <- ggplot(data = d_within, aes(x = Singleton, y = rt, color = task,
                            fill = task)) +
  geom_col(position = position_dodge(width = 0.4), width = 0.3, alpha = 0.3, size = .5) +
  coord_cartesian(ylim = c(650, 735)) +
  facet_wrap(
    ~ task,
    labeller = labeller(task = c(C = "Color task", L = "Location task"))
  ) +
  geom_errorbar(data = d_within, position = dodge, aes(ymin = rt-ci, ymax = rt+ci), width = 0, size = .5) +
  geom_line(position = dodge, aes(group = task, linetype = task), size = .5) +
  geom_point(data = d_within, position = dodge, size = 2, shape = 21) +
  scale_color_manual(values = c("C" = "#258da5", "L" = "#C15F1E"), # Customize colors
                     labels = c("C" = "Color", "L" = "Location")) + # Custom labels
  scale_fill_manual(
    values = c(
      "C" = "#7CBAC9",  # ~40% mix with white
      "L" = "#E0A67A"   # ~40% mix with white
    ),
    labels = c("C" = "Color", "L" = "Location")
  ) +
  scale_linetype_manual(
    values = c(
      "C" = "solid",  # ~40% mix with white
      "L" = "dashed"   # ~40% mix with white
    ),
    labels = c("C" = "Color", "L" = "Location")
  ) +
  labs(x = "Distractor type", y = "Response times (ms)") + # Change legend title
  theme_Publication(base_size = 10, text_size = 10) +
  theme(
    legend.position = "none"
    )

# Accuracy:
d <- raw %>%
  filter(Phase == "Rewarded", rt > 150, rt < 1800, !ID %in% c(exclusions), Block_num > 2) %>%
  dplyr::summarise(ACC = mean(correct), n = n(), .by = c(task, ID, Singleton))

# Set hypothesis matrix:
d$Singleton <-
  factor(d$Singleton, levels = c("High", "Low", "Absent"))
d$Task <- ifelse(d$task == "L", -.05, .05)
contrasts(d$Singleton) <- contr.hypothesis(HcRep)

# Fit model:
fitACC <-
  glmer(
    ACC ~ Singleton*task + (1 | ID),
    control = glmerControl(optimizer = 'bobyqa'),
    data = d,
    weights = n,
    family = binomial()
  )

# Model summary
summary(fitACC)

# Save for later
fit_list[["ACC_VST"]] <- fitACC

# Report task accuracy ----
# Load data:
d_Rep <-
  raw %>%
  filter(Phase %in% c("Rewarded", "Report"), !ID %in% exclusions) %>%
  group_by(ID) %>%
  dplyr::mutate(tmp = 1:n()) %>%
  ungroup() %>%
  filter(tmp > 30) %>%
  select(-tmp) %>%
  filter(Phase == "Report") %>%
  summarise(correct = sum(correct), n = n(), .by = c(task, ID))%>%
  mutate(
    ACC = correct / n,
  )

# Set contrast
d_Rep$Task <- ifelse(d_Rep$task == "L", -.05, .05)

# Fit model
fit_report <- glmer(
  ACC ~ task + (1 | ID),
  control = glmerControl(optimizer = 'bobyqa'),
  data = d_Rep,
  weights = n,
  family = binomial()
)

# Model summary
summary(fit_report)

# Save for later
fit_list[["REPORT"]] <- fit_report


# Awareness tests ----
# Load data:
awareness <- raw %>% filter(Phase == "Awareness1" |
                                                Phase == "Confidence1" |
                                                Phase == "Confidence2" |
                                                Phase == "Awareness2", !ID %in% c(exclusions)) %>% select(ID, task, Phase, response) %>%
  dplyr::mutate(response = as.numeric(response)/100,
                n = n(),
                response = transform_data(response, n=n)) %>%
  summarise(response = max(response),
            n = max(n),
            .by = c(ID, Phase, task)) %>%
  tidyr::pivot_wider(names_from = "Phase", values_from = "response")


awareness$Task <- factor(awareness$task, levels = c("C", "L"))
contrasts(awareness$Task) <- contr.sum(2)
colnames(contrasts(awareness$Task)) <- c("c_vs_L")

# Between-group differences
# Contingency Belief
summary(betareg(Awareness1~Task, data = awareness)) # What is the percentage of bonus trials you earned with each color distractor?

# Contingency awareness
summary(betareg(Awareness2~Task, data = awareness)) # What is the percentage of bonus trials you earned with each color distractor?

# Validity of the measures:
# 1: Is response confidence associated with awareness
# Color:
# Contingency Belief
summary(betareg(Awareness1~scale(Confidence1)*Task, data = awareness)) # No significant association

# Contingency awareness
summary(betareg(Awareness2~scale(Confidence2)*Task, data = awareness)) # Positively related

# 2: Are both measures of awareness associated with each other?
summary(betareg(Awareness2~scale(Awareness1)*Task, data = awareness)) # Positively related

# Plot:
fig2b <- marginaleffects::avg_predictions(betareg(Awareness2~scale(Awareness1)*Task, data = awareness), condition = c("Awareness1", "Task"),
                                                    by =c("Awareness1", "Task")) %>%
  ggplot(aes(Awareness1, estimate, color = Task, fill = Task))+
  geom_line(aes(linetype = Task), size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2, color = NA)+
  geom_point(data = awareness, aes(x = Awareness1, y = Awareness2, color = task), alpha = .3, size = 1)+
  scale_color_manual(values = c("C" = "#258da5", "L" = "#C15F1E"), # Customize colors
                     labels = c("C" = "Color", "L" = "Location")) + # Custom labels
  scale_fill_manual(
    values = c(
      "C" = "#7CBAC9",  # ~40% mix with white
      "L" = "#E0A67A"   # ~40% mix with white
    ),
    labels = c("C" = "Color", "L" = "Location")
  ) +
  scale_linetype_manual(
    values = c(
      "C" = "solid",  # ~40% mix with white
      "L" = "dashed"   # ~40% mix with white
    ),
    labels = c("C" = "Color", "L" = "Location")
  ) +
  labs(y = "Comntingency Awareness",
       x = "Contingency Belief") +
  geom_hline(yintercept = c(0.5), linetype = "dashed") +
  geom_vline(xintercept = c(0.5), linetype = "dashed") +
  scale_y_continuous(limits = c(0,1)) +
  theme_Publication(text_size = 10) +
  theme(legend.position = "none")

# Is the previous interaction explained away by individual differences in awareness?
# Only adding scale(Awareness2) as a covariate

# Load data
d <- raw %>%
  filter(Phase == "Rewarded", rt > 150, rt < 1800, correct == 1, !ID %in% c(exclusions), Block_num > 2)

d <- left_join(d, awareness[,c("Awareness2", "ID")])

# Set contrasts
d$Singleton <-
  factor(d$Singleton, levels = c("High", "Low", "Absent"))
d$Task <- factor(d$task, levels = c("C", "L"))
contrasts(d$Task) <- contr.sum(2)
colnames(contrasts(d$Task)) <- c("c_vs_L")
d$log_RT <- log(d$rt)
contrasts(d$Singleton) <- contr.hypothesis(HcRep)

# Fit model
fit <-
  lmer(
    log_RT ~ Singleton * Task + scale(Awareness2) * Singleton + (Singleton | ID),
    control = lmerControl(optimizer = 'bobyqa'),
    data = d
  )

# Raw summary
summary(fit)

# one-tailed p-value for VMAC:Task
pt(summary(fit)$coefficients[6, 4], summary(fit)$coefficients[6, 3], lower.tail = F) 

# Save for later
fit_list[["fit_aw"]] <- fit

# Plot results:
comps<-comparisons(fit, variables = "Singleton", newdata = datagrid(Task = unique,
                                                                Singleton = c("High", "Low"),
                                                                ID = NA,
                                                                Awareness2 = seq(0, 1, .01)),
               comparison = function(hi, lo) exp(hi + (sigma(fit)^2)/2) - exp(lo + (sigma(fit)^2)/2),
               re.form = NA) %>%
  filter(contrast == "Low, High") %>%
  mutate(estimate = -estimate,
         conf.high = -conf.high,
         conf.low = - conf.low)

raw_data <- d %>%
  summarise(rt = mean(rt), .by = c("Task", "ID", "Singleton")) %>%
  spread(Singleton, rt) %>%
  mutate(VMAC = High - Low) %>%
  left_join(., awareness, by = "ID") %>%
  mutate(Task = task)


vmac_sum <- summarySE(
  data = raw_data,
  measurevar = "VMAC",
  groupvars = "Task"
) %>%
  mutate(CI_VMAC = ci) %>%
  select(Task, VMAC, CI_VMAC)

awareness_sum <- summarySE(
  data = raw_data,
  measurevar = "Awareness2",
  groupvars = "Task"
) %>%
  mutate(CI_aw = ci) %>%
  select(Task, Awareness2, CI_aw)

sum_data <- left_join(vmac_sum, awareness_sum, by = "Task")

fig2c <- ggplot(comps, aes(Awareness2, estimate, color = Task)) +
  geom_errorbar(data = sum_data, aes(y = VMAC, ymin = VMAC-CI_VMAC, ymax = VMAC+CI_VMAC), width = 0) +
  geom_errorbarh(data = sum_data, aes(y = VMAC,xmin = Awareness2-CI_aw, xmax = Awareness2+CI_aw), height = 0) +
  geom_line(aes(linetype = Task), size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Task), color = NA, alpha = .2) +
  geom_point(data = raw_data, aes(y = VMAC), alpha = .3, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(data = sum_data, aes(y = VMAC, x = Awareness2, fill = Task), alpha = 1, size = 2.2, shape = 21) +
  scale_color_manual(values = c("C" = "#258da5", "L" = "#C15F1E"), # Customize colors
                     labels = c("C" = "Color", "L" = "Location")) + # Custom labels
  scale_fill_manual(
    values = c(
      "C" = "#7CBAC9",  # ~40% mix with white
      "L" = "#E0A67A"   # ~40% mix with white
    ),
    labels = c("C" = "Color", "L" = "Location")
  ) +
  scale_linetype_manual(
    values = c(
      "C" = "solid",  # ~40% mix with white
      "L" = "dashed"   # ~40% mix with white
    ),
    labels = c("C" = "Color", "L" = "Location")
  ) +
  scale_y_continuous(breaks = seq(-75, 125, 25), limits = c(-75, 125)) +
  scale_x_continuous(breaks = seq(0, 1, .1), limits = c(0.0, 1)) +
  labs(y = "VMAC effect (ms)", x = "Contingency Awareness") +
  theme_Publication(text_size = 10) +
  theme(legend.position = "none",
        plot.margin = margin(t = .1,  # Top margin
                             r = 2,  # Right margin
                             b = .1,  # Bottom margin
                             l = 2,
                             unit = "cm"))

# Combine all plots in one figure:
fig2 <- ggarrange(
  ggarrange(
    fig2a,
    fig2b,
    nrow = 1,
    ncol = 2,
    labels = c("A", "B")
    ),
  fig2c,
  nrow = 2,
  labels = c("", "C")
  )

# Save figure with customs dimensions:
if(!file.exists("Output/plots/fig2.pdf")) {
  ggsave(plot=fig2, device = "pdf", filename=here("Output/plots/fig2.pdf"),
         units = "cm", height = 14, width = 17)
}
