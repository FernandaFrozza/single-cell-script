library(survival)
library(broom)
library(dplyr)
library(ggplot2)

#Hazard ratios (HRs) and 95% confidence intervals were estimated using multivariable Cox proportional hazards regression. Results are displayed as a forest plot with HRs on a logarithmic scale. P-values are categorized into predefined significance levels to facilitate clinical interpretation.

# Multivariable Cox proportional hazards model
cox_multi <- coxph(
  Surv(Survival.days, Survival) ~ 
    Age +
    Sex +
    Stage +
    Autoimmunity +
    Hypothyroidism +
    Smoking.Status +
    Race +
    Histological.Subtype,
  data = matriz_cox,
  x = TRUE
)

cox_df <- tidy(
  cox_multi,
  exponentiate = TRUE,
  conf.int = TRUE
)

cox_df_plot <- cox_df %>%
  
  # Categorize p-values for visualization
  mutate(
    p_cat = case_when(
      p.value < 0.01 ~ "< 0.01",
      p.value < 0.05 ~ "0.01 – 0.05",
      TRUE ~ "≥ 0.05"
    )
  )

ggplot(cox_df_plot, aes(x = estimate, y = term)) +
  
  # Reference line
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  
  # Confidence intervals
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.25,
    color = "black"
  ) +
  
  # Hazard ratio points
  geom_point(
    aes(color = p_cat),
    size = 6
  ) +
  
  # Log-scale HR
  scale_x_log10() +
  
  # Discrete p-value colors
  scale_color_manual(
    values = c(
      "< 0.01" = "#B2182B",
      "0.01 – 0.05" = "#EF8A62",
      "≥ 0.05" = "#4D4D4D"
    ),
    name = "p-value"
  ) +
  
  # Clean theme
  theme_classic(base_family = "Arial") +
  
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  ) +
  
  labs(
    x = "Hazard Ratio (log scale)",
    y = NULL,
    title = "Multivariate Cox Regression – Overall Survival"
  )