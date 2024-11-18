# Author: Francisco Garre-Frutos

# Script with helper functions

#source("scripts/load_data.R") # Automatically load data

# Loading relevant packages ----

if (!require(pacman)) {
  install.packages(pacman)
  library(pacman)
}

p_load(dplyr, Hmisc, grid, ggthemes, simr, car)

# Defining functions ----

# Function to filter data. The decisions are based on Garre-Frutos et al. (2024):
filter_data <- function(data,
                        f.absent = T,
                        acc = T,
                        sd_filter = NULL,
                        fixed = T,
                        phase = "Rewarded",
                        experiment = "e1") {

  if (experiment == "eM") {
    data$correct <- data$Accuracy
    data$rt <- data$RT 
    }
  if (experiment == "e1") data$correct <- data$acc
  
  acc_ex <- data %>%
    filter(Phase %in% c("Rewarded")) %>%
    group_by(ID) %>%
    dplyr::summarise(mean_ACC = mean(correct)) %>% filter(mean_ACC < .7) %>% pull(ID)
  
  out <- data %>%
    filter(Phase %in% phase, !ID %in% acc_ex)
  
  if (f.absent)
    out <- out[which(out$Singleton != "Absent"), ]
  
  if (fixed)
    out <- out[which(out$rt > 150 & out$rt < 1800), ]
  
  if (is.numeric(sd_filter) &
      ifelse(is.null(sd_filter), F, sd_filter %in% 1:3)) {
    out <- out %>%
      group_by(ID) %>%
      dplyr::mutate(
        high_rt = mean(rt, na.rm = T) + sd(rt, na.rm = T) * sd_filter,
        low_rt = mean(rt, na.rm = T) - sd(rt, na.rm = T) * sd_filter
      ) %>%
      ungroup()  %>%
      filter(rt > low_rt, rt < high_rt) %>%
      select(-c(high_rt, low_rt))
  }
  
  if (acc)
    out <- out[which(out$correct == 1), ]
  
  return(out)
}


# Funtions to create epochs of blocks:
create_epochs <- function(blocks, epoch = 2) {
  vapply(blocks, function(x, e = epoch) {
    ceiling(x / e)
  }, FUN.VALUE = numeric(1))
  
}

# Theme used in the plots
theme_Publication <-
  function(base_size = 12,
           base_family = "sans",
           text_size = 11,
           lengend_pos) {
    (
      theme_foundation(base_size = base_size, base_family = base_family)
      + theme(
        plot.title = element_text(
          face = "bold",
          size = rel(1.2),
          #hjust = 0.5
        ),
        text = element_text(size = text_size),
        panel.background = element_rect(colour = NA),
        plot.background = element_rect(colour = NA),
        panel.border = element_rect(colour = NA),
        axis.title = element_text(face = "bold", size = rel(1)),
        axis.title.y = element_text(angle = 90, vjust = 2),
        axis.title.x = element_text(vjust = -0.2),
        axis.text = element_text(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(),
        panel.grid.major = element_line(colour = "#f0f0f0"),
        panel.grid.minor = element_blank(),
        legend.key = element_rect(colour = NA),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.size = unit(0.2, "cm"),
        #legend.margin = margin(0, "cm"),
        #legend.title = element_text(face="italic"),
        plot.margin = unit(c(10, 5, 5, 5), "mm"),
        strip.background = element_rect(colour = NA, fill = NA),
        strip.text = element_text(face = "bold")
      )
    )
  }

rm_legend <- function(p){
  p + theme(legend.position = "none")
}

# Make new data. Argument I represents the number of participants. 
# abs indicates whether to include absent singleton trials and mod is a lmer model fit.
makeModel <- function(I, mod.meta, form = formula(rt ~ Singleton * Group + (Singleton | ID)), abs = F) {
  # Function to produce a vector of shuffled Singleton conditions
  makeSingleton <- function(vec, I) {
    v <- vector()
    for (i in 1:I) {
      for (j in 1:10) {
        v <- c(v, sample(vec))
      }
    }
    return(v)
  }
  
  # Update Singleton condition
  singleton <- makeSingleton(c(rep(c("High", "Low"), each = 10), rep("Absent", length.out = 4)), I) # Sample singleton conditions
  N <- 24
  d <- data.frame(
    Block = rep(rep(1:10, each = N), length.out = I * 10 * N), # 10 blocks of N trials and I participants
    Singleton = factor(singleton, levels = c("High", "Low", "Absent")),
    Group = factor(rep(c("A", "B"), 
                       each = (I / 2) * 10 * N)), # Adjusting for two levels
    ID = as.factor(rep(paste0("ID", 1:I), each = 10 * N)) # Participant ID
  )
  
  # Setting repeated contrasts with hypr package
  if (abs) {
    HcRep <- hypr(
      VMAC = High ~ Low,
      AC = Low ~ Absent,
      levels = c("High","Low", "Absent")
    )
  } else {
    d <- d[which(d$Singleton != "Absent"),] # Filter out 'Absent'
    d$Singleton <- ifelse(d$Singleton == "High", .5, -.5) # Contrast for Singleton
  }
  HcRep <- hypr(
    AB = A ~ B,
    levels = c("A", "B")
  )
  contrasts(d$Group) <- HcRep # Set repeated contrasts for Group variable
  
  # Modify fixef for repeated contrasts in a 2x2 design
  old_fixef <- fixef(mod.meta)
  
  # Create a new vector of fixef
  new_fixef <- old_fixef[1]
  
  # Singleton VMAC adjustment
  new_fixef["Singleton"] <- old_fixef[2]  # Adjust as necessary
  new_fixef["GroupAB"] <- 0  # No difference between A and B
  
  # Set the interaction
  new_fixef["Singleton:GroupAB"] <- old_fixef["Singleton:Group"]  # Adjusted interaction
  
  # Simulate the model using the updated fixed effects
  m <- makeLmer(form, fixef = new_fixef, 
                VarCorr = matrix(VarCorr(mod.meta)$`ID`, ncol = 2),
                sigma = sigma(mod.meta), data = d)
  
  return(m)
}

# Function to perform power analysis
sim_lme4 <- function(dat, model) {
  isSing <- T
  while (isSing) { # Fit the model again if is singular
    dat$rt <- doSim(model) # Simulate data
    f <- model@call[["formula"]]
    tmp <- lmer(
      formula = f,
      data = dat,
      control = lmerControl(optimizer = "bobyqa")
    )
    isSing <- isSingular(tmp) # Check singularity
  }

  return(list(
    inter = pt(summary(tmp)$coefficients[4, 4], summary(tmp)$coefficients[4, 3], lower.tail = F) < .05 # one-sided test only for VMAC x Group
  )) # VMAC x Group interaction
}

# Fuction to get the sample VMAC effect:
get_effects <- function(d, exp = "e2") {
  if (exp == "e2") {
    return(
      d %>%
        group_by(group, ID, Singleton) %>%
        dplyr::summarise(RT = mean(rt)) %>%
        spread(Singleton, RT) %>%
        mutate(VMAC = High - Low) %>%
        ungroup() %>%
        drop_na()
    )
  } else {
    return(
      d %>%
        group_by(ID, Singleton) %>%
        dplyr::summarise(RT = mean(rt)) %>%
        spread(Singleton, RT) %>%
        mutate(VMAC = High - Low) %>%
        ungroup() %>%
        drop_na()
    )
  }
}

# Functions to plot:

# This function gets predictions for both models presented in the main text to plot
# d indicate raw data fed in a model for each phase. m indicate models for each phase
# Epoch is set to TRUE, which use epochs of 2 for raw data
# returns a list with 1) averaged raw data by Phase, Block (or Epoch) and Singleton
# and 2) model predictions for both models

get_predictions <- function(d, m, epoch = T, exp = "e1") {
  library(marginaleffects)
  library(dplyr)
  
  if (exp == "e1") {
    fPreds <- function(exp, model) {
      return(datagrid(
        Singleton = unique,
        Block = seq(1, 12, .01),
        ID = NA,
        model = model
      ))
    }
    rawW <- c("Phase", "Block", "Singleton")
    rawWse <- c("Phase", "Block", "Singleton")
    rawB = NULL
    avarageG <- c("Singleton", "Block")
    if (epoch) {
      d$Block <- create_epochs(d$Block)
    }
  } else {
    fPreds <- function(exp, model) {
      return(datagrid(
        Singleton = unique,
        Group = unique,
        Block_num = seq(1, 12, .01),
        ID = NA,
        model = model
      ))
    }
    rawW <- c("Group", "Phase", "Block", "Singleton")
    rawWse <- rawW[-1]
    rawB <- rawW[1]
    avarageG <- c("Group", "Singleton", "Block_num")
    if (epoch) {
      d$Block <- create_epochs(d$Block_num)
    }
  }
  
  raw_data <- d %>%
    dplyr::summarise(estimate = mean(rt), .by = c("ID", rawW)) %>%
    Rmisc::summarySEwithin(
      data = .,
      betweenvars = rawB,
      measurevar = "estimate",
      withinvars = rawWse,
      idvar = "ID"
    ) %>% mutate(Phase = "Rewarded")
  
  raw_data$Singleton <-
    factor(raw_data$Singleton, levels = c("High", "Low", "Absent"))
  raw_data$Block <- as.numeric(raw_data$Block) * 2 - .5
  
  model_preds <- predictions(
    m,
    newdata = fPreds(exp = exp, m),
    re.form = NA,
    transform = \(x) exp(x + (sigma(m) ^ 2) / 2),
    by = avarageG
  )
  
  model_preds$Singleton <-
    factor(model_preds$Singleton, levels = c("High", "Low", "Absent"))
  
  if (exp == "e2") {
    model_preds$Block <- model_preds$Block_num
    return(list(
      raw = raw_data %>% mutate(Group = ifelse(
        Group == .5, "Instructions", "No instructions"
      )),
      mod = model_preds %>% mutate(Group = ifelse(
        Group == .5, "Instructions", "No instructions"
      ))
    ))
    
  }
  
  return(list(raw = raw_data, mod = model_preds))
}

# Intermediate function to get VMAC and AC effects in get_comparisons()
get_raw_effect <- function(d, epoch = T, exp = "e1") {
  if (exp == "e1") {
    rawW <- c("Phase", "Block", "Singleton")
    rawWse <- c("Phase", "Block")
    rawB = NULL
    if (epoch) {
      d$Block <- create_epochs(d$Block)
    }
  } else {
    rawW <- c("Group", "Phase", "Block", "Singleton")
    rawWse <- rawW[c(-1, -2, -4)]
    rawB <- rawW[1]
    if (epoch) {
      d$Block <- create_epochs(d$Block_num)
    }
  }
  
  effs <- d %>%
    dplyr::summarise(RT = mean(rt), .by = c("ID", rawW)) %>%
    spread(Singleton, RT) %>%
    mutate(VMAC = High - Low, AC = Low - Absent) %>%
    ungroup() %>%
    drop_na()
  
  return(
    rbind(
      summarySEwithin2(
        data = effs,
        measurevar = "VMAC",
        betweenvars = rawB,
        withinvars = rawWse,
        idvar = "ID"
      ) %>%
        dplyr::rename("estimate" = "VMAC") %>% mutate(Block = as.numeric(Block), Effect = "VMAC") %>%
        .[, -5],
      summarySEwithin2(
        data = effs,
        measurevar = "AC",
        betweenvars = rawB,
        withinvars = rawWse,
        idvar = "ID"
      ) %>% dplyr::rename("estimate" = "AC") %>%
        mutate(Block = as.numeric(Block), Effect = "AC") %>% .[, -5]
    )
  )
  
}

# This function gets conditional effects for each contrasts in both models presented in the main text to plot
# d indicate raw data fed in a model for each phase. m indicate models for each phase
# Epoch is set to TRUE, which use epochs of 2 for raw data
# returns a list with 1) averaged raw data by Phase, Block (or Epoch) and Singleton
# and 2) model predictions for both models

get_comparisons <- function(d, m, epoch = T, exp = "e1") {
  if (exp == "e1") {
    fPreds <- function(exp, model) {
      return(datagrid(
        Block = seq(1, 12, .01),
        ID = NA,
        model = model
      ))
    }
    if (epoch) {
      d$Block <- create_epochs(d$Block)
    }
  } else {
    fPreds <- function(exp, model) {
      return(datagrid(
        Group = unique,
        Block_num = seq(1, 12, .01),
        ID = NA,
        model = model
      ))
    }
  }
  
  raw_data <- get_raw_effect(d, epoch, exp)
  raw_data$Effect <-
    factor(raw_data$Effect, levels = c("VMAC", "AC"))
  raw_data$Block <- as.numeric(raw_data$Block) * 2 - .5
  
  
  model_comps <- comparisons(
    m,
    variables = list(Singleton = "revsequential"),
    newdata = fPreds(exp = exp, m),
    re.form = NA,
    transform_pre = \(hi, lo) exp(hi + (sigma(m) ^ 2) / 2) - exp(lo +
                                                                   (sigma(m) ^ 2) / 2),
  ) %>% mutate(Effect = rep(c("VMAC", "AC"), each = nrow(.) / 2), Phase = "Rewarded")
  
  model_comps$Effect <-
    factor(model_comps$Effect, levels = c("VMAC", "AC"))
  
  return(list(raw = raw_data, mod = model_comps))
}

# Functions to get average and SE for raw data in mixed designs. The source code is in: :
summarySEwithin2 <- function (data = NULL,
                              measurevar,
                              betweenvars = NULL,
                              withinvars = NULL,
                              idvar = NULL,
                              na.rm = FALSE,
                              conf.interval = 0.95,
                              .drop = TRUE)
{
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop = FALSE], FUN = is.factor, FUN.VALUE = logical(1))
  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message(
      "Automatically converting the following non-factors to factors: ",
      paste(nonfactorvars, collapse = ", ")
    )
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }
  datac <- summarySE(
    data,
    measurevar,
    groupvars = c(betweenvars, withinvars),
    na.rm = na.rm,
    conf.interval = conf.interval,
    .drop = .drop
  )
  
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop = .drop)
  measurevar_n <- paste(measurevar, "_norm", sep = "")
  ndatac <- summarySE(
    ndata,
    measurevar_n,
    groupvars = c(betweenvars, withinvars),
    na.rm = na.rm,
    conf.interval = conf.interval,
    .drop = .drop
  )
  nWithinGroups <- prod(vapply(ndatac[, withinvars, drop = FALSE], FUN = nlevels, FUN.VALUE = numeric(1)))
  correctionFactor <- sqrt(nWithinGroups / (nWithinGroups - 1))
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor
  merge(datac, ndatac)
}
normDataWithin <- function(data = NULL,
                           idvar,
                           measurevar,
                           betweenvars = NULL,
                           na.rm = FALSE,
                           .drop = TRUE) {
  # Measure var on left, idvar + between vars on right of formula.
  data.subjMean <- plyr::ddply(
    data,
    c(idvar, betweenvars),
    .drop = .drop,
    .fun = function(xx, col, na.rm) {
      c(subjMean = mean(xx[, col], na.rm = na.rm))
    },
    measurevar,
    na.rm
  )
  
  # Put the subject means with original data
  data <- merge(data, data.subjMean)
  
  # Get the normalized data in a new column
  measureNormedVar <- paste(measurevar, "_norm", sep = "")
  data[, measureNormedVar] <- data[, measurevar] - data[, "subjMean"] +
    mean(data[, measurevar], na.rm = na.rm)
  
  # Remove this subject mean column
  data$subjMean <- NULL
  
  return(data)
}

# Plot results from multiverse analysis. "d" indicates multiverse estimates to show in the top row, and
# "spec" are the specifications to show in the bottom row. "ord" refers is a character to order estimates
# and specifications.

multiverse_plot <- function(d,
                            spec,
                            ord,
                            y_name = "Split-half reliability",
                            minimun_t = T,
                            spec_seg = NA,
                            cor = F) {
  # Ordering estimates and decisions:
  raw_sp <-
    d$spearmanbrown
  decisions_raw <-
    spec[, 2] %>% unnest(cols = c("data")) # unorded decisions
  decisions <-
    decisions_raw %>% mutate(tmp = raw_sp) %>% arrange(tmp) %>% select(-tmp)
  
  est_ord <- d %>%
    arrange(spearmanbrown) %>% mutate(Universe = 1:nrow(decisions))
  est_ord <- cbind(est_ord, decisions)
  
  if (cor) {
    est_ord$spearmanbrown <- est_ord$cor
    est_ord$SB_low <- est_ord$lwr.ci
    est_ord$SB_high <- est_ord$ipr.ci
  }
  
  # Creating formatted data for dashboard:
  final_d <- decisions %>% mutate(Universe = 1:nrow(decisions)) %>%
    dplyr::rename(
      "Fixed" = "RT_fix_cutof",
      "Relative" = "RT_sd_cutoff",
      "Averaging" = "averaging_method",
      "Length" = "NBlocks",
      "First_trials" = "Two_Trials"
    ) %>%
    mutate(
      Averaging = ifelse(Averaging == "mean", "Mean", "Median"),
      Relative = paste("Filter", Relative, "SD"),
      logtransform = ifelse(logtransform == "yes", "Log", "Raw"),
      First_trials = ifelse(First_trials == "TRUE", "Remove first", "Don't remove"),
      Length = paste(Length, "Blocks"),
      Fixed = ifelse(Fixed == "TRUE", "Fixed cut-off", "No cut-off")
    ) %>%
    mutate(Relative = ifelse(Relative == "Filter 0 SD", "No filter", Relative)) %>%
    rstatix::convert_as_factor(Relative, Fixed, Length, First_trials)
  
  final_d$Length <-
    factor(final_d$Length, levels = c("6 Blocks", "12 Blocks"))
  
  final_d %>%
    pivot_longer(!Universe, names_to = "Bigdecision", values_to = "Decision") -> final_d
  
  final_d$Decision <-
    factor(
      final_d$Decision,
      levels = c(
        "Median",
        "Mean",
        "Don't remove",
        "Remove first",
        "Fixed cut-off",
        "No cut-off",
        "6 Blocks",
        "12 Blocks",
        "Raw",
        "Log",
        "No filter",
        "Filter 2 SD",
        "Filter 2.5 SD",
        "Filter 3 SD"
      )
    )
  
  # Top plot
  plot_est <- est_ord %>%
    ggplot(aes(x = Universe, y = spearmanbrown, )) +
    geom_point(size = 1) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    geom_ribbon(aes(ymin = SB_low, ymax = SB_high),
                alpha = .2,
                color = NaN) +
    labs(y = y_name) +
    scale_x_continuous(breaks = c(1, seq(0, 128, 10)[-1], nrow(decisions)),
                       limits = c(1, 128)) +
    theme_Publication(text_size = 11) +
    theme(
      axis.title.x = element_blank(),
      legend.position = 'none',
      plot.margin = margin(1, 1, 0, 1, unit = "cm"),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x = element_blank()
    )
  
  if (minimun_t) {
    plot_est <- plot_est +
      coord_cartesian(ylim = c(0, .9)) +
      annotate(
        "text",
        x = 20,
        y = .8,
        label = "Minimum threshold",
        size = 3.5
      ) +
      geom_hline(yintercept = .7, linetype = "dashed")
  }
  
  if (!cor) {
    plot_est <- plot_est + scale_y_continuous(breaks = seq(0, 1, .2))
  } else {
    plot_est <- plot_est + geom_hline(yintercept = 0, linetype = "dashed")
  }
  
  # Bottom plot
  dashboard <- ggplot(data = final_d, aes(x = Universe, y = Decision)) +
    facet_grid(Bigdecision ~ ., scales = "free_y") +
    #scale_colour_brewer(palette = 'Set1') +
    geom_point(shape = 108, size = 3) +
    #annotate("rect",xmin = 273-1, xmax = 273 + 1, ymin = -Inf, ymax = Inf,
    #alpha = .3, color = "gray")+
    labs(x = "Specification") +
    scale_y_discrete(position = "left") +
    theme_Publication(text_size = 11) +
    scale_x_continuous(breaks = c(1, seq(0, 128, 10)[-1], nrow(decisions)),
                       limits = c(1, 128)) +
    theme(
      legend.position = 'none',
      strip.background = element_blank(),
      strip.text.y = element_blank(),
      plot.margin = margin(0.3, 1, 1, 1, unit = "cm"),
      text = element_text(size = 11)
    )
  
  # Joint plot
  if (!is.na(spec_seg)) {
    plot_est <- plot_est +
      annotate(
        "segment",
        xend = spec_seg,
        x = spec_seg,
        y = Inf,
        yend = -Inf,
        colour = "blue",
        alpha = .3,
        linewidth = 1.5
      )
    dashboard <- dashboard +
      annotate(
        "segment",
        x = spec_seg,
        xend = spec_seg,
        y = -Inf,
        yend = Inf,
        colour = "blue",
        alpha = .3,
        linewidth = 1.5
      )
  }
  
  return((plot_est / dashboard + plot_layout(heights = c(4, 6))))
}


# Get effects for participants at different quantiles in the distribution of contingency
# to plot the results of the contingencty awaraness analysis:

get_effectCont <- function(d,
                           probs = c(.60),
                           epoch = T) {
  if ("Block" %in% colnames(d)) {
    d$Block_num <- d$Block
  }
  if (epoch)
    d$Block <- create_epochs(d$Block_num, epoch = 2)
  
  effs <- d %>%
    mutate(Contingency = case_when(
      Contingency < quantile(.$Contingency, probs)[1] ~ "Low Contingency",
      T ~ "High Contingency"
    )) %>%
    group_by(Contingency, ID, Block, Singleton) %>%
    dplyr::summarise(RT = mean(rt)) %>%
    spread(Singleton, RT) %>%
    mutate(VMAC = High - Low, AC = Low - Absent) %>%
    ungroup() %>%
    drop_na()
  
  Quants <- bind_rows(
    Rmisc::summarySEwithin(
      data = effs,
      measurevar = "VMAC",
      withinvars = c("Block"),
      betweenvars = c("Contingency"),
      idvar = "ID"
    )  %>%
      mutate(
        Block = as.numeric(Block),
        Effect = "VMAC",
        estimate = VMAC
      ) %>% select(-VMAC),
    Rmisc::summarySEwithin(
      data = effs,
      measurevar = "AC",
      withinvars = c("Block"),
      betweenvars = c("Contingency"),
      idvar = "ID"
    )  %>%
      mutate(
        Block = as.numeric(Block),
        Effect = "AC",
        estimate = AC
      ) %>% select(-AC)
  ) %>% mutate(Block = Block * 2 - .5)
  return(Quants)
}