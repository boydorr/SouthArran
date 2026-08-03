#### Creating Models
library(lme4)
library(lmtest)
library(vegan)

# getting PCA results
pca_result <- rda(results_nmds$comm_used)  
pc_scores <- scores(pca_result, display = "sites", choices = 1:2) %>%
  as.data.frame() %>%
  rename(PC1 = PC1, PC2 = PC2)
pc_scores$sample_id <- results_nmdss$meta_used$sample_id


df_modeling <- results_nmds$meta_used %>%
  left_join(pc_scores %>% select(sample_id, PC1), by = "sample_id") %>%
  rename(protection = protection_level,
         phi        = Mean_phi,
         depth      = Depth,
         station    = GrabSite_station) %>% 
  mutate(period     = factor(period, levels = c("Before", "After")),
         protection = factor(protection),
         station    = factor(station))

m1 <- lmer(PC1 ~ period * protection + phi + depth + (1 | station), data = df_modeling, REML = FALSE)
m2 <- lmer(PC1 ~ period + protection + phi + depth + (1 | station), data = df_modeling, REML = FALSE)
m3 <- lmer(PC1 ~ period + phi + depth + (1 | station), data = df_modeling, REML = FALSE)
m4 <- lmer(PC1 ~ protection + phi + depth + (1 | station), data = df_modeling, REML = FALSE)
m5 <- lmer(PC1 ~ phi + depth + (1 | station), data = df_modeling, REML = FALSE)
#LRT
lrtest(m1, m2)   # interaction needed?
lrtest(m2, m3)   # protection needed, given period + phi + depth?
lrtest(m2, m4)   # period needed, given protection + phi + depth?
lrtest(m4, m5) 
lrtest(m3, m5) 



# LRT on T survey sites 
# getting PCA results
pca_result <- rda(results_nmds_period_Tsites$comm_used)
pc_scores <- scores(pca_result, display = "sites", choices = 1:2) %>%
  as.data.frame() %>%
  rename(PC1 = PC1, PC2 = PC2)
pc_scores$sample_id <- results_nmds_period_Tsites$meta_used$sample_id

df_modeling <- results_nmds_period_Tsites$meta_used %>%
  left_join(pc_scores %>% select(sample_id, PC1), by = "sample_id") %>%
  rename(protection = protection_level,
         phi        = Mean_phi,
         depth      = Depth,
         station    = GrabSite_station) %>%
  mutate(period     = factor(period, levels = c("Before", "After")),
         protection = factor(protection),
         station    = factor(station))

m1_T <- lmer(PC1 ~ period * protection + phi + depth + (1 | station), data = df_modeling, REML = FALSE)
m2_T <- lmer(PC1 ~ period + protection + phi + depth + (1 | station), data = df_modeling, REML = FALSE)
m3_T <- lmer(PC1 ~ period + phi + depth + (1 | station), data = df_modeling, REML = FALSE)
m4_T <- lmer(PC1 ~ protection + phi + depth + (1 | station), data = df_modeling, REML = FALSE)
m5_T <- lmer(PC1 ~ phi + depth + (1 | station), data = df_modeling, REML = FALSE)

#LRT
lrtest(m1_T, m2_T)   # interaction needed?
lrtest(m2_T, m3_T)   # protection needed, given period + phi + depth?
lrtest(m2_T, m4_T)   # period needed, given protection + phi + depth?
lrtest(m4_T, m5_T)
lrtest(m3_T, m5_T)



# testing 
pc_scores_traits <- results_PCA_abundance_protection_Tsites$pca_out$scores

df_modeling_traits <- spatial_result$meta %>%
  filter(grepl("^T", GrabSite_station)) %>%  
  left_join(
    pc_scores_traits %>% select(sample_id, PC1),
    by = "sample_id"
  ) %>%
  rename(
    protection = protection_level,
    phi = Mean_phi,
    depth = Depth,
    station = GrabSite_station
  ) %>%
  mutate(
    protection = factor(protection),
    station = factor(station))

m1_T_traits <- lmer(PC1 ~ period * protection + phi + depth + (1 | station), data = df_modeling_traits, REML = FALSE)
m2_T_traits <- lmer(PC1 ~ period + protection + phi + depth + (1 | station), data = df_modeling_traits, REML = FALSE)
m3_T_traits <- lmer(PC1 ~ period + phi + depth + (1 | station), data = df_modeling_traits, REML = FALSE)
m4_T_traits <- lmer(PC1 ~ protection + phi + depth + (1 | station), data = df_modeling_traits, REML = FALSE)
m5_T_traits <- lmer(PC1 ~ phi + depth + (1 | station), data = df_modeling_traits, REML = FALSE)

#LRT
lrtest(m1_T_traits, m2_T_traits)   # interaction needed?
lrtest(m2_T_traits, m3_T_traits)   # protection needed, given period + phi + depth?
lrtest(m2_T_traits, m4_T_traits)   # period needed, given protection + phi + depth?
lrtest(m4_T_traits, m5_T_traits)
lrtest(m3_T_traits, m5_T_traits)