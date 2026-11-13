# ANOVA + Tukey + plot with letters (BLUP_Ocurrence_upper ~ twist)
library(ggplot2)
library(dplyr)
library(emmeans)     # pairwise Tukey + letters
library(multcomp)    # needed by emmeans::cld()
library(rlang)     # !!sym



df0<-read.csv("Q2_TraitVariation_Plasticity_Ancestry/trait_pi.csv")

twist<-read.csv("Q3_GWAS/TWIST_Ancestry.csv")

df<-merge(df0, twist, by= "ID", all.y=TRUE )


df$twist<-df$Ocurrence_upper_Potri.015G022300
df$twist<-as.factor(df$twist)
cols_sp<- c("dodgerblue2","gray80","darkolivegreen2")


# One-way ANOVA + Tukey HSD for multiple traits in a loop
# Ensure factor
df$twist <- as.factor(df$twist)

# Columns to analyze/plot
y_vars <- c(
  "Vermont_Ocurrence_upper",
  "Virginia_Ocurrence_upper",
  "Plasticity_Ocurrence_upper"
)

# Containers
plots <- list()
aov_summaries <- list()
tukey_results <- list()
letters_tables <- list()

for (v in y_vars) {
  # 1) Fit ANOVA for the current variable
  fml <- as.formula(paste(v, "~ twist"))
  fit <- aov(fml, data = df)
  aov_summaries[[v]] <- summary(fit)
  
  # 2) Tukey pairwise
  emm <- emmeans(fit, ~ twist)
  tuk <- pairs(emm, adjust = "tukey")
  tukey_results[[v]] <- as.data.frame(tuk)
  
  # 3) Compact Letter Display
  lett <- cld(emm, adjust = "tukey", Letters = letters)
  lett <- as.data.frame(lett) |> mutate(.group = gsub("\\s+", "", .group))
  
  # y positions per group (a bit above group max)
  ypos <- df |>
    group_by(twist) |>
    summarise(y = max(.data[[v]], na.rm = TRUE), .groups = "drop")
  rng <- range(df[[v]], na.rm = TRUE)
  offset <- 0.05 * (rng[2] - rng[1])
  labs_df <- left_join(lett, ypos, by = "twist") |>
    mutate(y = y + offset)
  
  # 4) Plot for this variable
  p <- ggplot(df, aes(x = twist, y = !!sym(v))) +
    geom_boxplot(outlier.shape = NA, fill="white") +
    geom_jitter(aes(fill = twist), color="black", shape=21, width = 0.2, alpha = 0.8, size = 3, show.legend = FALSE) +
    geom_text(data = labs_df, aes(x = twist, y = y, label = .group),
              vjust = 0, fontface = "bold", inherit.aes = FALSE) +
    labs(x = "twist", y = v, title = paste0("ANOVA with Tukey HSD — ", v)) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(color = "black", size = 10),
          axis.text.y = element_text(color = "black", size = 10),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 11),
          plot.title = element_text(size = 14),
          plot.subtitle = element_text(size = 14, color = "black"),
          panel.grid = element_blank())+
    scale_fill_manual(values = cols_sp) +
    scale_color_manual(values = cols_sp)
  
    
  
  plots[[v]] <- p
  letters_tables[[v]] <- labs_df
}

# Show plots in the viewer
for (p in plots) print(p)

# (Optional) Save all plots to a multi-page PDF
pdf("Q3_GWAS/anova_tukey_by_trait.pdf", width = 4, height = 4)
for (p in plots) print(p)
dev.off()






