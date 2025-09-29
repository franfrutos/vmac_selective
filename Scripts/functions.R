# Author: Francisco Garre-Frutos

# Script with helper functions

# Loading relevant packages ----

if (!require(pacman)) {
  install.packages(pacman)
  library(pacman)
}

p_load(dplyr, Hmisc, grid, ggthemes, simr, car)

# Defining functions ----

# Function to filter data in Garre-Frutos et al. (2025)
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

# Transform data to avoid exact 0s and 1s
transform_data <- function(y, epsilon = 0.5, n) {
  (y * (n - 1) + epsilon) / n
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


format_interaction <- function(vars) {
  # vars: vector of characters, e.g. c("VMAC") o c("VMAC","Group")
  # format each element as \mathrm{…}
  formatted <- sprintf("\\mathrm{%s}", vars)
  # if more than one element add " \\times " in between
  if (length(formatted) == 1) {
    return(formatted)
  } else {
    return(paste(formatted, collapse = " \\times "))
  }
}

# Function to report predictors from (G)LMMs:
report_coef <- function(model, term, term_name, one.sided = F, digits = 3) {
  model_sum <- summary(model)
  if (!term %in% rownames(model_sum$coefficients)) {
    stop(sprintf("The term '%s' cannot be found in the model.", term))
  }
  coef_est <- model_sum$coefficients[term, "Estimate"]
  p_colname <- grep("^Pr\\(>\\|.*\\)$", colnames(model_sum$coefficients), value = TRUE)
  p_value  <- model_sum$coefficients[term, p_colname]
  p_output <-  if (p_value < .001) "p < 0.001" else paste0("p = ", format(round(p_value, digits), nsmall = digits))



  if ("t value" %in% colnames(model_sum$coefficients)) {
    test_stat <- model_sum$coefficients[term, "t value"]
    df_value  <- model_sum$coefficients[term, "df"]
    if (one.sided) {
      p_value <- pt(test_stat, df_value, lower.tail = F)
      p_output <-  if (p_value < .001) "p < 0.001" else paste0("p = ", format(round(p_value, digits), nsmall = digits))
    }
    test_name <- "t"
    name <- format_interaction(term_name)
    latex_code <- sprintf("$\\beta_{%s} = %s, \\; %s_{%s} = %s, \\; %s$",
                          name,
                          format(round(coef_est, digits), nsmall = digits),
                          test_name,
                          format(round(df_value, 0)),
                          format(round(test_stat, 2), nsmall = 2),
                          p_output)
  } else if ("z value" %in% colnames(model_sum$coefficients)) {
    test_stat <- model_sum$coefficients[term, "z value"]
    test_name <- "z"
    latex_code <- sprintf("$\\beta_{%s} = %s, \\; %s = %s, \\; %s$",
                          name,
                          format(round(coef_est, digits), nsmall = digits),
                          test_name,
                          format(round(test_stat, 2), nsmall = 2),
                          format(round(p_value, digits), nsmall = digits))
  } else {
    stop("There is no t or z value column in the model.")
  }
  return(latex_code)
}

# Function to select p-values
get_least_non_significant <- function(model) {
  coefs <- summary(model)$coefficients


  # Extract p-values and their names
  pvals <- coefs[, 4]
  names(pvals) <- rownames(coefs)

  # Filter for p > 0.05
  nonsig <- pvals[pvals > 0.05]

  if (length(nonsig) == 0) {
    message("All fixed effects are significant (p <= 0.05).")
    return(NULL)
  }

  # Select the one with the smallest p-value
  selected <- which.min(nonsig)

  return(nonsig[selected])
}
