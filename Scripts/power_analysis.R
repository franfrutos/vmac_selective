# Author: Francisco Garre-Frutos
# Date: 20/06/2024

# Power analysis from Garre-Frutos et al. (2024)

# Set random seed:
set.seed(666)

# Set-up ----
if (!require(pacman)) {
  install.packages(pacman)
  library(pacman)
}

p_load(here, dplyr, readr, lme4, lmerTest, hypr, simr, pbapply, future, parallel, binom, ggplot2, marginaleffects)

source(here("Scripts/functions.R")) # Load data and all relevant functions

# Loading raw data from Garre-Frutos et al. (2024b) from github repo
dlist <- list(
  "exp2" = read_csv("https://raw.githubusercontent.com/franfrutos/VMAC_awareness/main/Input/experiment_2.csv")[,-1]
)

d <- filter_data(dlist[["exp2"]], f.absent = F,
                               sd_filter = 2,
                               experiment = "e2") %>% select(Singleton, rt, ID, group, Block_num) %>%
  mutate(ID = paste0(ID, "exp2"),
         group = ifelse(group == "A", "Instructions", "No instructions"),
         Exp = "exp2",
         Block = Block_num
         )

d <- d %>% filter(Singleton != "Absent") %>%
  mutate(Singleton = factor(Singleton, levels = c("High", "Low")),
         group = factor(group, levels = c("Instructions", "No instructions")),
         Exp = factor(Exp, levels = c("expM", "exp1", "exp2"))
         )

d$Singleton <- ifelse(d$Singleton == "High", 0.5, -0.5)
d$Group <- ifelse(d$group == "Instructions", 0.5, -0.5)


fit <- lmer(log(rt)~Singleton*Group +(Singleton|ID),control=lmerControl(optimizer='bobyqa'),
                 data = d[d$Block > 2,])

summary(fit) # Model used for simulation


# Power-analysis based on previous study ----
# VMAC main effect and VMAC x Block x Group interaction
# If the file with the power analysis result is detected, do not run power analysis
if (!file.exists(here("Output/powerVMAC.rds"))) {
  # Setting simulation params
  ncls <- 10 # Clusters for parallel computing
  steps <- seq(20*2, 100*2, 20) # Steps for power curve
  nsim <- 1000 # number of simulations
  seed <- 123 # seed for reproducible output
  li <- list() # store results
  for (i in 1:length(steps)) {
    print(paste("N =", steps[i], "Participants.", steps[i]/2, "participants per group."))
    set.seed(seed)
    model <- makeModel(steps[i], mod.meta = fit, abs = F) # Create a model with the specified structure
    dat <- getData(model) # Get data from the model to simulate
    cl <- makeCluster(ncls) # Set clusters
    clusterExport(cl, c("sim_lme4", "doSim", "lmer", "which", "dat",
                        "model", "fit", "lmerControl", "isSingular")) # Export objects and functions
    power <- pbreplicate(n = nsim, expr = sim_lme4(dat, model), simplify = T, cl = cl)
    stopCluster(cl)
    li[[paste(i)]]$Int <- binom::binom.confint(sum(unlist(power)), nsim, method = "exact")
    print(li[[i]])
  }
  
  # Get power df
  power <- do.call(rbind,do.call(rbind, li))
  power$N <- rep(steps/2, length.out = nrow(power))
  power$Effect <- rep(c("Planned-contrast (one-tailed)"), each = length(steps))
  power$Effect <- factor(power$Effect, levels = c("Planned-contrast (one-tailed)"))
  saveRDS(power, file = "Output/powerVMAC.rds")
  
} else {
  power <- readRDS(here("Output/powerVMAC.rds"))
}

(plot_power <- 
  ggplot(data=power %>% filter(Effect == "Planned-contrast (one-tailed)", N < 101),
         aes(x = N, y = mean, color = Effect, fill = Effect)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = unique(power$N), labels = unique(power$N)) +
  coord_cartesian(ylim = c(0, 1))+
  scale_y_continuous(labels = paste0(seq(0, 1, .1)*100, "%"), breaks = seq(0, 1, .1)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 1/3, color = NA)+
  geom_hline(yintercept = .80, linetype = "dashed") +
  labs(y = "Power", x = "Number of participants per group")+
  scale_fill_manual(values=c("darkblue")) +
  scale_color_manual(values=c("darkblue")) +
  theme_Publication(base_size = 12) +
  theme(plot.margin = margin(0, 0, 10, 10),
        legend.position = "none"))

ggsave(plot = plot_power,
       filename = here("Output/plots/powerCurve.png"),
       height = 12,
       width = 17,
       units = "cm",
       dpi = 900)

# Predictions from the simulated model (is the model doing what is supposed to do?) ----
# This code is just to check if the above simulation power analysis produce the expected output.
# Run the model with more participants to have more precise estimates.
simM <- makeModel(80*2, mod.meta = fit) # Fitted model with the hypothetical structure for 80 participants per group (using simR::makeLmer)
simD <- getData(simM) # Data from the hypothetical model
simD$rt <- doSim(simM) # Simulate response variable

simF <- lmer(rt~Singleton*Group +(Singleton|ID),control=lmerControl(optimizer='bobyqa'),
            data = simD) # Fit new model

summary(simF) # Coefficients with custom contrasts in the Singleton*Group interaction (two-tailed by default).

# Model predictions in the response scale (ms) by Group
avg_comparisons(simF, variables = "Singleton", re.form = NULL, newdata = simD, by = "Group",
                comparison = \(hi, lo) exp(hi + (sigma(simF) ^ 2) / 2) - exp(lo + (sigma(simF) ^ 2) / 2))

# Raw data visualization
(pred_plot <- simD %>%
  # Raw data are in the log-scale. We need to back-transform mean RTs assuming a lognormal likelihood
  summarise(RT = mean(exp(rt+(sigma(simF)^2)/2)),
            .by = c("Group", "ID", "Singleton")) %>%
  mutate(Group = case_when(
    Group == "A"~"Distractor color",
    Group == "B"~"Distracor location",
  ),
  Group = factor(Group, levels = c("Distractor color", "Distracor location")),
  Singleton = ifelse(Singleton == .5, "High", "Low")) %>%
  tidyr::spread(Singleton, RT) %>%
  mutate(VMAC = High - Low) %>%
  ggplot(aes(Group, VMAC)) +
  geom_point(position = position_jitter(width = .1), alpha = .1) +
  stat_summary(size = .3, fun.data = "mean_cl_normal") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(breaks = seq(-100, 100, 40)) +
  labs(y = "VMAC effect (ms)") +
  theme_Publication(base_size = 12))

ggsave(filename = here("Output/plots/modelPreds.png"),
       height = 15,
       width = 17,
       units = "cm",
       dpi = 900)
