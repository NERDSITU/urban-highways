###  Urban highways are barriers to social ties ###
###  regression on number of social connections crossing highways ###
### by sandorjuhasz ###


###### structure ######
# (1) data preparation -- for multiple models of the paper and SI
# (2) Table 1 regression for the main text
# (3) Table SI3 regressions for the Supplementary Information
# (4) Table SI4 regressions for the Supplementary Information
# (5) Fig. SI7 and related regressions for the Supplementary Information
# (6) Fig. SI6 regressions to collect coefficients and p-values
# (7) rebuttal -- regressions that reflect the setting of Fig. 2 main text
# (8) gravity model verification -- census tract based -- for rebuttal
# (9) gravity model verification -- cell based -- for rebuttal
# (10) distance of tract behind model settings -- data export for rebuttal figure
# (11) null model ties as dependent variable
######



### packages
library(arrow)
library(data.table)
library(tidyr)
library(dplyr)
library(MASS)
library(pscl)
library(stargazer)
library(gravity)
library(interplot)
library(cowplot)
library(fixest)



###### (1) data preparation -- for multiple models of the paper and SI ######

### NOTE : aim is to construct the dataframe 'df' with all information for all possible combinations of tracts
### 'df' is filtered directly for each model specification at place

# parquet file from Luca -- census tract-census tract real/null model connections
df1_obs <- open_dataset("../data/regressions/df_pairwise_averaged_part.parquet/") %>%
  collect() %>%
  data.table()

# census tract - cbsacode - cbsa short name
metro_names <- fread("../data/regressions/cbsacode_shortname_tracts.csv", drop = 1)

# add replaced variable on num_motorways
df3_num_motorway <- open_dataset("../data/regressions/df_pairwise_num_motorways.parquet/") %>%
  collect() %>%
  data.table()

df <- merge(
  dplyr::select(df1_obs, -num_motorways),
  df3_num_motorway,
  by = c("cbsacode", "tract_user_id1", "tract_user_id2"),
  all.x = TRUE,
  all.y = FALSE
)



# add missing rows for all possible tract-tract models
fn <- list.files("../data/regressions/missing_pairs_nomotorway_links_part/")
missings <- list()
for(f in 1:length(fn)){
  # print(fn[f])
  missings[[f]] <- open_dataset(paste0("../data/regressions/missing_pairs_nomotorway_links_part/", fn[f])) %>%
    collect() %>%
    data.table()
}

# combine and save
missings <- rbindlist(missings)


# missing edges from Luca
df2_missings <- missings %>%
  rename(
    tract_user_id1 = tract1,
    tract_user_id2 = tract2
    #motorway_real = num_motorways
  ) %>%
  data.table()
df2_missings$count_real <- 0
df2_missings$count_null <- 0

# select key columns
key_columns <- c(
  "tract_user_id1",
  "tract_user_id2",
  "cbsacode",
  "count_real",
  "count_null",
  "distance",
  "num_motorways",
  "population2014_tract1",
  "population2014_tract2",
  "income2014_tract1",
  "income2014_tract2",
  "white2014_tract1",
  "black2014_tract1",
  "asian2014_tract1",
  "white2014_tract2",
  "black2014_tract2",
  "asian2014_tract2"
)
df <- df[, ..key_columns]
df2_missings <- df2_missings[, ..key_columns]

# combine to one large dataframe
df <- rbind(df, df2_missings)
nrow(df)



# drop AB-BA ties
df <- subset(df, tract_user_id1 < tract_user_id2)



# ADD other barriers
df4_ob <- open_dataset("../data/regressions/df_pairwise_num_barriers_all_part.parquet/") %>%
  collect() %>%
  data.table()

df <- merge(
  df,
  df4_ob,
  by.x = c("cbsacode", "tract_user_id1", "tract_user_id2"),
  by.y = c("cbsacode", "tract1", "tract2"),
  all.x = TRUE,
  all.y = FALSE
)



# variable manipulation -- social connections
df$social <- df$count_real
df$social_null <- df$count_null
df$log_social <- log10(df$count_real)
df$log_social_p1 <- log10(1 + df$count_real)
df$social01 <- ifelse(df$count_real > 0, 1, 0)
df$log_social_null <- log10(1 + df$social_null)
df$social_null01 <- ifelse(df$count_null > 0, 1, 0)

# variable manipulation -- distance
df$log_dist <- log10((df$distance / 1000) + 1)

# variable manipulation -- motorway
df$log_motorways <- log10(df$num_motorways + 1)
df$highway01 <- ifelse(df$num_motorways > 0, 1, 0)
table(df$social01, df$highway01)
df$log_num_barriers <- log10(df$num_barriers + 1)

# variable manipulation -- segregation
df <- df %>% 
  mutate_at(c("income2014_tract1", "income2014_tract2"), ~replace_na(.,0)) %>%
  mutate(income_abs_diff = log10(abs((income2014_tract1 / 1000) - (income2014_tract2 / 1000)))) %>%
  data.table()
df$income_abs_diff[is.infinite(df$income_abs_diff)==1] <- 0

# variable manipulation -- racial composition
df[, "highest_share1"] <- apply(dplyr::select(df, white2014_tract1, black2014_tract1, asian2014_tract1), 1, max)
df$highest_share1[is.na(df$highest_share1)==TRUE] <- 0
df$white2014_tract1[is.na(df$white2014_tract1)==TRUE] <- 0
df[, "highest_share2"] <- apply(dplyr::select(df, white2014_tract2, black2014_tract2, asian2014_tract2), 1, max)
df$highest_share2[is.na(df$highest_share2)==TRUE] <- 0
df$white2014_tract2[is.na(df$white2014_tract2)==TRUE] <- 0
df$dominant_white1 <- ifelse(df$highest_share1 == df$white2014_tract1, 1, 0)
df$dominant_white2 <- ifelse(df$highest_share2 == df$white2014_tract2, 1, 0)
df$same_race <- ifelse(df$dominant_white1 == df$dominant_white2, 1, 0)

# variable manipulation -- population
df$pop1 <- df$population2014_tract1 / 1000
df$pop2 <- df$population2014_tract2 / 1000
df$log_pop1 <- log10(df$population2014_tract1 / 1000)
df$log_pop2 <- log10(df$population2014_tract2 / 1000)
df$pp <- log10(df$pop1 * df$pop2 + 1)
df$pop_sum <- df$pop1 + df$pop2






###### (2) Table 1 regression for the main text ######

### OLS -- sample of observed AND null model connected tracts ###
reg_df <- subset(df, (social > 0) | (social_null > 0))
table(reg_df$social01, reg_df$social_null01)

summary(ols_p1 <- lm(log_social_p1 ~ log_dist + pp + as.factor(cbsacode), data = reg_df))
summary(ols_p2 <- lm(log_social_p1 ~ log_motorways + log_dist + pp + as.factor(cbsacode), data = reg_df))
summary(ols_p3 <- lm(log_social_p1 ~ log_num_barriers + log_dist + pp + as.factor(cbsacode), data = reg_df))
summary(ols_p4 <- lm(log_social_p1 ~ income_abs_diff + log_dist + pp + as.factor(cbsacode), data = reg_df))
summary(ols_p5 <- lm(log_social_p1 ~ same_race + log_dist + pp + as.factor(cbsacode), data = reg_df))
summary(ols_p6 <- lm(log_social_p1 ~ log_motorways + log_num_barriers + income_abs_diff + same_race + log_dist + pp + as.factor(cbsacode), data = reg_df))


stargazer(ols_p1,
          ols_p2,
          ols_p3,
          ols_p4,
          ols_p5,
          ols_p6,
          dep.var.labels = "Number of social connections (log)",
          dep.var.caption = "",
          covariate.labels = c("Nr highways (log)", "Nr barriers (log)", "Income abs difference", "Racial homophily", "Distance (log)", "Population (log product)"),
          omit = c("cbsacode"),
          add.lines=list(c("Metro FE", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes")),
          omit.stat = c("f", "ser", "aic"),
          #out="../outputs/paper_ols_counterfact.tex"
          out="../outputs/si_all_barriers_ols_counterfact.html"
)


### same table in a fixest setting for robustness checks
reg_df <- subset(df, (social > 0) | (social_null > 0))
table(reg_df$social01, reg_df$social_null01)

summary(ols_p1 <- feols(log_social_p1 ~ log_dist + pp | cbsacode, vcov = "HC1", data = reg_df))
summary(ols_p2 <- feols(log_social_p1 ~ log_motorways + log_dist + pp | cbsacode, vcov = "HC1", data = reg_df))
summary(ols_p3 <- feols(log_social_p1 ~ log_num_barriers + log_dist + pp | cbsacode, vcov = "HC1", data = reg_df))
summary(ols_p4 <- feols(log_social_p1 ~ income_abs_diff + log_dist + pp | cbsacode, vcov = "HC1", data = reg_df))
summary(ols_p5 <- feols(log_social_p1 ~ same_race + log_dist + pp | cbsacode, vcov = "HC1", data = reg_df))
summary(ols_p6 <- feols(log_social_p1 ~ log_motorways + log_num_barriers + income_abs_diff + same_race + log_dist + pp | cbsacode, vcov = "HC1", data = reg_df))
etable(ols_p1, ols_p2, ols_p3, ols_p4, ols_p5, ols_p6)

summary(ols_p1 <- feols(log_social_p1 ~ log_dist + pp | cbsacode, data = reg_df))
summary(ols_p2 <- feols(log_social_p1 ~ log_motorways + log_dist + pp | cbsacode, data = reg_df))
summary(ols_p3 <- feols(log_social_p1 ~ log_num_barriers + log_dist + pp | cbsacode, data = reg_df))
summary(ols_p4 <- feols(log_social_p1 ~ income_abs_diff + log_dist + pp | cbsacode, data = reg_df))
summary(ols_p5 <- feols(log_social_p1 ~ same_race + log_dist + pp | cbsacode, data = reg_df))
summary(ols_p6 <- feols(log_social_p1 ~ log_motorways + log_num_barriers + income_abs_diff + same_race + log_dist + pp | cbsacode, data = reg_df))
etable(ols_p1, ols_p2, ols_p3, ols_p4, ols_p5, ols_p6)





###### (3) Table SI3 regressions for the Supplementary Information ######

### OLS -- log(t_ij) -- sample of tracts with ONLY observed social connections ###
reg_df <- subset(df, social > 0)
reg_df <- subset(reg_df, is.infinite(log_social)==FALSE)
summary(si_m2 <- lm(log_social ~ log_motorways + log_num_barriers + income_abs_diff + same_race + log_dist + pp + as.factor(cbsacode), data = reg_df))



### OLS -- log(1 + t_ij) -- sample of all possible census tract pairs ###
reg_df <- df
ols_all_poss <- feols(log_social_p1 ~ log_motorways + log_num_barriers + income_abs_diff + same_race + log_dist + pp | cbsacode, data = reg_df)
print(ols_all_poss)
### manual transfer to stargazer-made latex table ###



### PPML -- t_ij -- sample of observed AND null model connected tracts ###
reg_df <- subset(df, (social > 0) | (social_null > 0))
reg_df$cbsacode_char <- as.character(reg_df$cbsacode)
table(reg_df$social01, reg_df$social_null01)

summary(pm1 <- ppml(
  dependent_variable = "social",
  distance = "log_dist",
  additional_regressors = c("log_motorways", "log_num_barriers", "income_abs_diff", "same_race", "pp", "cbsacode_char"),
  data = reg_df
))
### manual transfer to stargazer-made latex table ###



### PPML -- t_ij > 0 -- sample of tracts with ONLY observed social connections ###
reg_df <- subset(df, social > 0)
reg_df$cbsacode_char <- as.character(reg_df$cbsacode)

summary(pm2 <- ppml(
  dependent_variable = "social",
  distance = "log_dist",
  additional_regressors = c("log_motorways", "log_num_barriers", "income_abs_diff", "same_race", "pp", "cbsacode_char"),
  data = reg_df
))
### manual transfer to stargazer-made latex table ###



### collector stargazer table for the SI -- needs manual preparation because of the multiple package not compatible with stargazer ###
stargazer(ols_p6,
          si_m2,
          si_m2,
          si_m2,
          si_m2,
          covariate.labels = c("Nr highways (log)", "Nr barriers (log)", "Income abs difference", "Racial homophily", "Distance (log)", "Population (log product)"),
          omit = c("cbsacode"),
          add.lines=list(c("Metro FE", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes")),
          omit.stat = c("f", "ser", "aic"),
          out="../outputs/si_for_manual_update.html"
          #out="../outputs/si_for_manual_update.tex"
)






###### (4) Table SI4 regressions for the Supplementary Information ######

### OLS -- sample of observed AND null model connected tracts ###
reg_df <- subset(df, (social > 0) | (social_null > 0))
table(reg_df$social01, reg_df$social_null01)

summary(ols_interp1_base <- lm(log_social_p1 ~ log_dist + log_motorways + log_num_barriers + income_abs_diff + same_race + pp + as.factor(cbsacode), data = reg_df))
summary(ols_interp1 <- lm(log_social_p1 ~ log_dist * log_motorways + log_num_barriers + income_abs_diff + same_race + pp + as.factor(cbsacode), data = reg_df))



### OLS -- log(t_ij) -- sample of tracts with ONLY observed social connections ###
reg_df <- subset(df, social > 0)
reg_df <- subset(reg_df, is.infinite(log_social)==FALSE)

summary(ols_interp2_base <- lm(log_social ~ log_dist + log_motorways + log_num_barriers + income_abs_diff + same_race + pp + as.factor(cbsacode), data = reg_df))
summary(ols_interp2 <- lm(log_social ~ log_dist * log_motorways + log_num_barriers + income_abs_diff + same_race + pp + as.factor(cbsacode), data = reg_df))



### collector stargazer for the SI -- needs manual preparation ###
stargazer(ols_interp1_base,
          ols_interp1,
          ols_interp2_base,
          ols_interp2,
          #covariate.labels = c("Distance (log)", "Nr highways (log)", "Income abs difference", "Racial homophily", "Population (log product)", "Distance X Nr highways"),
          omit = c("cbsacode"),
          add.lines=list(c("Metro FE", "Yes", "Yes", "Yes", "Yes")),
          omit.stat = c("f", "ser", "aic"),
          out="../outputs/si_interactions.html"
          #out="../outputs/si_for_manual_update.tex"
)
### manual transfer ###






###### (5) Fig. SI7 and related regressions for the Supplementary Information ######

### OLS -- sample of observed AND null model connected tracts ###
reg_df <- subset(df, (social > 0) | (social_null > 0))
table(reg_df$social01, reg_df$social_null01)

summary(ols_interp1_base <- lm(log_social_p1 ~ log_dist + log_motorways + log_num_barriers + income_abs_diff + same_race + pp + as.factor(cbsacode), data = reg_df))
summary(ols_interp1 <- lm(log_social_p1 ~ log_dist * log_motorways + log_num_barriers + income_abs_diff + same_race + pp + as.factor(cbsacode), data = reg_df))

# interplot
title <- "interplot_ols_counterfactual"
file_name <- paste0("../plots/", title, ".png")
png(file_name, width=600, height=400, units = 'px')

ip1 <- interplot(m = ols_interp1, var1 = "log_motorways", var2 = "log_dist", size=3, rfill = "#A30000") +
  aes(color = "") +
  xlab("Distance") +
  ylab("Estimated coefficient for\nnumber of highways crossed") +
  geom_hline(yintercept = 0, linetype = "dashed", size=1.5) +
  ylim(-0.3, 0.2) +
  theme_cowplot(12) +
  theme(axis.text = element_text(size=15), axis.title=element_text(size=25)) +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5), labels = c(expression(10^0), "", expression(10^1), "", expression(10^2), "")) +
  theme(plot.margin = unit(c(2, 1, 0, 1), "cm")) +
  guides(fill="none") +
  theme(legend.position="none")
dev.off()



### OLS -- log(t_ij) -- sample of tracts with ONLY observed social connections ###
reg_df <- subset(df, social > 0)
reg_df <- subset(reg_df, is.infinite(log_social)==FALSE)

summary(ols_interp2_base <- lm(log_social ~ log_dist + log_motorways + log_num_barriers + income_abs_diff + same_race + pp + as.factor(cbsacode), data = reg_df))
summary(ols_interp2 <- lm(log_social ~ log_dist * log_motorways + log_num_barriers + income_abs_diff + same_race + pp + as.factor(cbsacode), data = reg_df))


# interplot
title <- "interplot_ols_non0"
file_name <- paste0("../plots/", title, ".png")
png(file_name, width=600, height=400, units = 'px')

ip2 <- interplot(m = ols_interp2, var1 = "log_motorways", var2 = "log_dist", size=3, rfill = "#A30000") +
  aes(color = "") +
  xlab("Distance") +
  ylab("Estimated coefficient for\nnumber of highways crossed") +
  geom_hline(yintercept = 0, linetype = "dashed", size=1.5) +
  ylim(-0.3, 0.2) +
  theme_cowplot(12) +
  theme(axis.text = element_text(size=15), axis.title=element_text(size=25)) +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5), labels = c(expression(10^0), "", expression(10^1), "", expression(10^2), "")) + 
  theme(plot.margin = unit(c(2, 1, 0, 1), "cm")) +
  guides(fill="none") +
  theme(legend.position="none")
dev.off()


# combined version
title <- "combined_interplots"
file_name <- paste0("../plots/", title, ".png")
png(file_name, width=950, height=400, units = 'px')
plot_grid(ip1, ip2, labels = c('A', 'B'), label_size = 30)
dev.off()



### collector stargazer for the SI -- needs manual preparation ###
stargazer(ols_interp1_base,
          ols_interp1,
          ols_interp2_base,
          ols_interp2,
          #covariate.labels = c("Distance (log)", "Nr highways (log)", "Income abs difference", "Racial homophily", "Population (log product)", "Distance X Nr highways"),
          omit = c("cbsacode"),
          add.lines=list(c("Metro FE", "Yes", "Yes", "Yes", "Yes")),
          omit.stat = c("f", "ser", "aic"),
          out="../outputs/si_interactions.html"
          #out="../outputs/si_for_manual_update.tex"
)




##### spagetti plot version for SI ######

### OLS -- sample of observed AND null model connected tracts ###
reg_df <- subset(df, (social > 0) | (social_null > 0))
table(reg_df$social01, reg_df$social_null01)

summary(ols_interp1 <- lm(log_social_p1 ~ log_dist * log_motorways + log_num_barriers + income_abs_diff + same_race + pp + as.factor(cbsacode), data = reg_df))
interp_main1 <- interplot(m = ols_interp1, var1 = "log_motorways", var2 = "log_dist", size=5, rfill = "#A30000")
x_coords_main <- ggplot_build(interp_main1)$data[[1]]$x
y_coords_main <- ggplot_build(interp_main1)$data[[1]]$y



# loop for the single-city-lines
cbs <- unique(reg_df$cbsacode)
ips <- list()
x_coords <- list()
y_coords <- list()

for(r in 1:length(cbs)){
  print(cbs[r])
  test_set <- subset(reg_df, cbsacode == cbs[r])
  ols_interp_piece <- lm(log_social_p1 ~ log_dist * log_motorways, data = test_set)
  
  ips[[r]] <- interplot(m = ols_interp_piece, var1 = "log_motorways", var2 = "log_dist", size=3, rfill = "Black") +
    aes(color = "") +
    xlab("Distance") +
    ylab("Coefficient for\nnr. highways crossed") +
    geom_hline(yintercept = 0, linetype = "dashed", size=1.5) +
    #xlim(0, 2) +
    #ylim(-1, 1) +
    theme_cowplot(12) +
    theme(axis.text = element_text(size=15), axis.title=element_text(size=25)) +
    #scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5), labels = c(expression(10^0), "", expression(10^1), "", expression(10^2), "")) + 
    theme(plot.margin = unit(c(2, 1, 0, 1), "cm")) +
    guides(fill="none") +
    theme(legend.position="none")
  
  x_coords[[r]] <- ggplot_build(ips[[r]])$data[[1]]$x
  y_coords[[r]] <- ggplot_build(ips[[r]])$data[[1]]$y
}


selected_version <- ggplot() +
  geom_line(data = data.frame(x = x_coords[[1]], y = y_coords[[1]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[2]], y = y_coords[[2]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[3]], y = y_coords[[3]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[4]], y = y_coords[[4]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[5]], y = y_coords[[5]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[6]], y = y_coords[[6]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[7]], y = y_coords[[7]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[8]], y = y_coords[[8]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[9]], y = y_coords[[9]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[10]], y = y_coords[[10]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[11]], y = y_coords[[11]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[12]], y = y_coords[[12]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[13]], y = y_coords[[13]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[14]], y = y_coords[[14]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[15]], y = y_coords[[15]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[16]], y = y_coords[[16]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[17]], y = y_coords[[17]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[18]], y = y_coords[[18]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[19]], y = y_coords[[19]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[20]], y = y_coords[[20]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[21]], y = y_coords[[21]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[22]], y = y_coords[[22]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[23]], y = y_coords[[23]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[24]], y = y_coords[[24]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[25]], y = y_coords[[25]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[26]], y = y_coords[[26]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[27]], y = y_coords[[27]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[28]], y = y_coords[[28]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[29]], y = y_coords[[29]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[30]], y = y_coords[[30]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[31]], y = y_coords[[31]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[32]], y = y_coords[[32]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[33]], y = y_coords[[33]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[34]], y = y_coords[[34]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[35]], y = y_coords[[35]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[36]], y = y_coords[[36]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[37]], y = y_coords[[37]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[38]], y = y_coords[[38]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[39]], y = y_coords[[39]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[40]], y = y_coords[[40]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[41]], y = y_coords[[41]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[42]], y = y_coords[[42]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[43]], y = y_coords[[43]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[44]], y = y_coords[[44]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[45]], y = y_coords[[45]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[46]], y = y_coords[[46]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[47]], y = y_coords[[47]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[48]], y = y_coords[[48]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[49]], y = y_coords[[49]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[50]], y = y_coords[[50]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords_main, y = y_coords_main), aes(x, y), size = 4, color = "#A30000") +
  xlab("Distance in km (log)") +
  ylab("Coefficient for\nnr. highways crossed") +
  geom_hline(yintercept = 0, linetype = "dashed", size=1.5) +
  ylim(-0.5, 0.5) +
  theme_cowplot(12) +
  theme(axis.text = element_text(size=15), axis.title=element_text(size=25)) +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5), labels = c(expression(10^0), "", expression(10^1), "", expression(10^2), "")) +
  theme(plot.margin = unit(c(2, 1, 0, 1), "cm")) +
  guides(fill="none") +
  theme(legend.position="none")



### OLS -- log(t_ij) -- sample of tracts with ONLY observed social connections ###
reg_df <- subset(df, social > 0)
reg_df <- subset(reg_df, is.infinite(log_social)==FALSE)

summary(ols_interp2 <- lm(log_social ~ log_dist * log_motorways + log_num_barriers + income_abs_diff + same_race + pp + as.factor(cbsacode), data = reg_df))
interp_main2 <- interplot(m = ols_interp2, var1 = "log_motorways", var2 = "log_dist", size=5, rfill = "#A30000")
x_coords_main <- ggplot_build(interp_main2)$data[[1]]$x
y_coords_main <- ggplot_build(interp_main2)$data[[1]]$y



# loop for the single-city-lines
cbs <- unique(reg_df$cbsacode)
ips <- list()
x_coords <- list()
y_coords <- list()

for(r in 1:length(cbs)){
  print(cbs[r])
  test_set <- subset(reg_df, cbsacode == cbs[r])
  ols_interp_piece <- lm(log_social_p1 ~ log_dist * log_motorways, data = test_set)
  
  ips[[r]] <- interplot(m = ols_interp_piece, var1 = "log_motorways", var2 = "log_dist", size=3, rfill = "Black") +
    aes(color = "") +
    xlab("Distance") +
    ylab("Coefficient for\nnr. highways crossed") +
    geom_hline(yintercept = 0, linetype = "dashed", size=1.5) +
    theme_cowplot(12) +
    theme(axis.text = element_text(size=15), axis.title=element_text(size=25)) +
    theme(plot.margin = unit(c(2, 1, 0, 1), "cm")) +
    guides(fill="none") +
    theme(legend.position="none")
  
  x_coords[[r]] <- ggplot_build(ips[[r]])$data[[1]]$x
  y_coords[[r]] <- ggplot_build(ips[[r]])$data[[1]]$y
}


observed_version <- ggplot() +
  geom_line(data = data.frame(x = x_coords[[1]], y = y_coords[[1]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[2]], y = y_coords[[2]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[3]], y = y_coords[[3]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[4]], y = y_coords[[4]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[5]], y = y_coords[[5]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[6]], y = y_coords[[6]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[7]], y = y_coords[[7]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[8]], y = y_coords[[8]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[9]], y = y_coords[[9]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[10]], y = y_coords[[10]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[11]], y = y_coords[[11]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[12]], y = y_coords[[12]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[13]], y = y_coords[[13]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[14]], y = y_coords[[14]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[15]], y = y_coords[[15]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[16]], y = y_coords[[16]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[17]], y = y_coords[[17]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[18]], y = y_coords[[18]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[19]], y = y_coords[[19]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[20]], y = y_coords[[20]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[21]], y = y_coords[[21]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[22]], y = y_coords[[22]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[23]], y = y_coords[[23]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[24]], y = y_coords[[24]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[25]], y = y_coords[[25]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[26]], y = y_coords[[26]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[27]], y = y_coords[[27]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[28]], y = y_coords[[28]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[29]], y = y_coords[[29]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[30]], y = y_coords[[30]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[31]], y = y_coords[[31]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[32]], y = y_coords[[32]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[33]], y = y_coords[[33]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[34]], y = y_coords[[34]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[35]], y = y_coords[[35]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[36]], y = y_coords[[36]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[37]], y = y_coords[[37]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[38]], y = y_coords[[38]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[39]], y = y_coords[[39]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[40]], y = y_coords[[40]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[41]], y = y_coords[[41]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[42]], y = y_coords[[42]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[43]], y = y_coords[[43]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[44]], y = y_coords[[44]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[45]], y = y_coords[[45]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[46]], y = y_coords[[46]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[47]], y = y_coords[[47]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[48]], y = y_coords[[48]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[49]], y = y_coords[[49]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords[[50]], y = y_coords[[50]]), aes(x, y), color = "grey") +
  geom_line(data = data.frame(x = x_coords_main, y = y_coords_main), aes(x, y), size = 4, color = "#A30000") +
  xlab("Distance in km (log)") +
  ylab("Coefficient for\nnr. highways crossed") +
  geom_hline(yintercept = 0, linetype = "dashed", size=1.5) +
  ylim(-0.5, 0.5) +
  theme_cowplot(12) +
  theme(axis.text = element_text(size=15), axis.title=element_text(size=25)) +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5), labels = c(expression(10^0), "", expression(10^1), "", expression(10^2), "")) +
  theme(plot.margin = unit(c(2, 1, 0, 1), "cm")) +
  guides(fill="none") +
  theme(legend.position="none")



# combined version
title <- "combined_interplots_extra"
file_name <- paste0("../plots/", title, ".pdf")
pdf(file_name, width=12, height=5)  # Adjust width and height as needed
plot_grid(selected_version, observed_version, labels = c('A', 'B'), label_size = 30)
dev.off()






###### (6) Fig. SI6 regressions to collect coefficients and p-values ######

### OLS -- sample of observed AND null model connected tracts ###
reg_df <- subset(df, (social > 0) | (social_null > 0))
table(reg_df$social01, reg_df$social_null01)

cbsacodes <- unique(df$cbsacode)
coeff_tables <- list()

for(c in 1:length(cbsacodes)){
  # cbsa focus
  fdf <- subset(reg_df, cbsacode == cbsacodes[c])
  
  # running the model
  summary(model_stats <- lm(log_social_p1 ~ log_dist + log_motorways + log_num_barriers + income_abs_diff + same_race + pp, data = fdf))
  
  # put the coeffs into a dataframe
  coeff_tables[[c]] <- data.table(
    c("log_dist", "log_motorways", "log_num_barriers", "income_abs_diff", "same_race", "pp"),
    model_stats$coefficients[2:6],
    round(sqrt(diag(vcov(summary(model_stats)))), 3)[2:6],
    round(summary(model_stats)$coefficients[,4], 3)[2:6]
  )
  coeff_tables[[c]]$cbsacode <- cbsacodes[c]
}

# combine coeff dataframes and export
combined <- rbindlist(coeff_tables, use.names = TRUE, fill = TRUE)
colnames(combined) <- c("variable", "coeff", "se", "pvalue", "cbsacode")

combined <- merge(
  combined,
  unique(dplyr::select(metro_names, cbsacode, short_name)),
  by = "cbsacode",
  all.x = TRUE,
  all.y = FALSE
)

write.table(
  combined,
  "../outputs/si_coeff_table_main_models_ols_counterfactual_other_barriers.csv",
  sep=";",
  row.names = FALSE
)







###### (7) rebuttal -- regressions that reflect the setting of Fig. 2 main text ######
reg_df <- subset(df, (social > 0) & (num_motorways > 0))

summary(ols_p1 <- feols(log_social_p1 ~ log_dist + pp | cbsacode, vcov = "HC1", data = reg_df))
summary(ols_p2 <- feols(log_social_p1 ~ log_motorways + log_dist + pp | cbsacode, vcov = "HC1", data = reg_df))
summary(ols_p3 <- feols(log_social_p1 ~ log_num_barriers + log_dist + pp | cbsacode, vcov = "HC1", data = reg_df))
summary(ols_p4 <- feols(log_social_p1 ~ income_abs_diff + log_dist + pp | cbsacode, vcov = "HC1", data = reg_df))
summary(ols_p5 <- feols(log_social_p1 ~ same_race + log_dist + pp | cbsacode, vcov = "HC1", data = reg_df))
summary(ols_p6 <- feols(log_social_p1 ~ log_motorways + log_num_barriers + income_abs_diff + same_race + log_dist + pp | cbsacode, vcov = "HC1", data = reg_df))
etable(ols_p1, ols_p2, ols_p3, ols_p4, ols_p5, ols_p6)
### manual transfer ###






###### (8) gravity model verification -- census tract based -- for rebuttal ######

### OLS -- sample of observed AND null model connected tracts ###
t_cf <- subset(df, (social > 0) | (social_null > 0))
t_cf$log_pop1 <- log10(t_cf$population2014_tract1 + 1)
t_cf$log_pop2 <- log10(t_cf$population2014_tract2 + 1)

summary(
  grav_m1 <- feols(
    log_social_p1 ~ log_dist + log_pop1 + log_pop2 | cbsacode,
    vcov = "HC1",
    data = t_cf
  )
)



### OLS -- log(t_ij) -- sample of tracts with ONLY observed social connections ###
t_observed <- subset(df, social > 0)
t_observed <- subset(t_observed, is.infinite(log_social)==FALSE)
t_observed$log_pop1 <- log10(t_observed$population2014_tract1 + 1)
t_observed$log_pop2 <- log10(t_observed$population2014_tract2 + 1)

summary(
  grav_m2 <- feols(
    log_social ~ log_dist + log_pop1 + log_pop2 | cbsacode,
    vcov = "HC1",
    data = t_observed
  )
)


### OLS -- log(1 + t_ij) -- sample of all possible census tract pairs ###
t_all <- df
t_all$log_pop1 <- log10(t_all$population2014_tract1 + 1)
t_all$log_pop2 <- log10(t_all$population2014_tract2 + 1)

summary(
  grav_m3 <- feols(
    log_social_p1 ~ log_dist + log_pop1 + log_pop2 | cbsacode,
    vcov = "HC1",
    data = t_all
  )
)

etable(grav_m1, grav_m2, grav_m3)


### PPML -- t_ij -- sample of observed AND null model connected tracts ###
t_cf <- subset(df, (social > 0) | (social_null > 0))
t_cf$log_pop1 <- log10(t_cf$population2014_tract1 + 1)
t_cf$log_pop2 <- log10(t_cf$population2014_tract2 + 1)
t_cf$cbsacode_char <- as.character(t_cf$cbsacode)
#table(t_cf$social01, t_cf$social_null01)

summary(grav_ppml1 <- ppml(
  dependent_variable = "social",
  distance = "log_dist",
  additional_regressors = c("log_pop1", "log_pop2", "cbsacode_char"),
  data = t_cf
))



### PPML -- t_ij > 0 -- sample of tracts with ONLY observed social connections ###
t_observed <- subset(df, social > 0)
t_observed$log_pop1 <- log10(t_observed$population2014_tract1 + 1)
t_observed$log_pop2 <- log10(t_observed$population2014_tract2 + 1)
t_observed$cbsacode_char <- as.character(t_observed$cbsacode)

summary(grav_ppml2 <- ppml(
  dependent_variable = "social",
  distance = "log_dist",
  additional_regressors = c("log_pop1", "log_pop2", "cbsacode_char"),
  data = t_observed
))
### manual transfer ###






###### (9) gravity model verification -- cell based -- for rebuttal ######

# data from Luca 
cell_counts <- data.table(read_parquet("../data/regressions/gravity_all_cell_counts.parquet"))
cell_links <- open_dataset("../data/regressions/gravity_all_cell_links_part.parquet/") %>%
  collect() %>%
  data.table()


# combine for gravity regression based on cells at different sizes
cell_df <- merge(
  cell_links,
  cell_counts,
  by.x = c("cbsacode", "cell_size", "cell_id1"),
  by.y = c("cbsacode", "cell_size", "cell_id"),
  all.x = TRUE,
  all.y = FALSE
)
cell_df <- merge(
  cell_df,
  cell_counts,
  by.x = c("cbsacode", "cell_size", "cell_id2"),
  by.y = c("cbsacode", "cell_size", "cell_id"),
  all.x = TRUE,
  all.y = FALSE,
  suffixes = c("1", "2")
)


# undirected
cell_df <- subset(cell_df, cell_id1 < cell_id2)


# variable manipulation
cell_df$log_ties <- log10(cell_df$weight)
cell_df$log_dist <- log10(round(cell_df$distance_beeline_m / 1000) + 1)
cell_df$log_users1 <- log10(cell_df$num_users1)
cell_df$log_users2 <- log10(cell_df$num_users2)
cell_df$cbsacode_char <- as.character(cell_df$cbsacode)


# OLS
ds <- unique(cell_df$cell_size)
ds <- c(1000, 5000, 10000, 20000)
d_ols <- list()
d_ppml <- list()
for(d in 1:length(ds)){
  print(ds[d])
  temp <- subset(cell_df, cell_size == ds[d])
  
  d_ols[[d]] <- feols(log_ties ~ log_dist + log_users1 + log_users2 | cbsacode, data = temp)
  print(summary(d_ols[[d]]))
  
  d_ppml[[d]] <- summary(
    ppml(
      dependent_variable = "weight",
      distance = "log_dist",
      additional_regressors = c("log_users1", "log_users2", "cbsacode_char"),
      data = temp
    )
  )
  print(d_ppml[[d]])
  print(nrow(temp))
}



### all possible combinations w/ zeros -- work in progress
cell_counts <- data.table(read_parquet("../data/regressions/gravity_all_cell_counts.parquet"))
all_cell_links <- open_dataset("../data/regressions/gravity_all_cell_links_part.parquet/") %>%
  collect() %>%
  data.table()


# combine for gravity regression based on cells at different sizes
all_cell_df <- merge(
  all_cell_links,
  cell_counts,
  by.x = c("cbsacode", "cell_size", "cell_id1"),
  by.y = c("cbsacode", "cell_size", "cell_id"),
  all.x = TRUE,
  all.y = FALSE
)
all_cell_df <- merge(
  all_cell_df,
  cell_counts,
  by.x = c("cbsacode", "cell_size", "cell_id2"),
  by.y = c("cbsacode", "cell_size", "cell_id"),
  all.x = TRUE,
  all.y = FALSE,
  suffixes = c("1", "2")
)


# undirected
all_cell_df <- subset(all_cell_df, cell_id1 < cell_id2)


# variable manipulation
all_cell_df$log_ties <- log10(all_cell_df$weight + 1)
all_cell_df$log_dist <- log10(round(all_cell_df$distance_beeline_m / 1000) + 1)
all_cell_df$log_users1 <- log10(all_cell_df$num_users1)
all_cell_df$log_users2 <- log10(all_cell_df$num_users2)
all_cell_df$cbsacode_char <- as.character(all_cell_df$cbsacode)


# OLS
ds <- unique(all_cell_df$cell_size)
ds <- c(5000, 10000, 20000)
d_ols <- list()
d_ppml <- list()
for(d in 1:length(ds)){
  print(ds[d])
  temp <- subset(all_cell_df, cell_size == ds[d])
  
  d_ols[[d]] <- feols(log_ties ~ log_dist + log_users1 + log_users2 | cbsacode, data = temp)
  print(summary(d_ols[[d]]))
  
  
  
  d_ppml[[d]] <- summary(
    ppml(
      dependent_variable = "weight",
      distance = "log_dist",
      additional_regressors = c("log_users1", "log_users2", "cbsacode_char"),
      data = temp
    )
  )
  print(d_ppml[[d]])
  print(nrow(temp))
}



###### (10) distance of tract behind model settings -- data export for rebuttal figure ######
plot_df <- df
plot_df$eid <- seq(1, nrow(plot_df), 1)
plot_df <- plot_df %>%
  dplyr::select(eid, cbsacode, social, social_null, distance, num_motorways) %>%
  data.table()

write.table(
  plot_df,
  "../outputs/distance_across_model_settings.csv",
  sep=";",
  row.names = FALSE
)



###### (11) null model ties as dependent variable ######

# observed +null model ties
reg_df <- subset(df, (social > 0) | (social_null > 0))
table(reg_df$social01, reg_df$social_null01)

summary(ols_n1 <- feols(log_social_null ~ log_dist + pp | cbsacode, vcov = "HC1", data = reg_df))
summary(ols_n2 <- feols(log_social_null ~ log_motorways + log_dist + pp | cbsacode, vcov = "HC1", data = reg_df))
summary(ols_n3 <- feols(log_social_null ~ log_num_barriers + log_dist + pp | cbsacode, vcov = "HC1", data = reg_df))
summary(ols_n4 <- feols(log_social_null ~ income_abs_diff + log_dist + pp | cbsacode, vcov = "HC1", data = reg_df))
summary(ols_n5 <- feols(log_social_null ~ same_race + log_dist + pp | cbsacode, vcov = "HC1", data = reg_df))
summary(ols_n6 <- feols(log_social_null ~ log_motorways + log_num_barriers + income_abs_diff + same_race + log_dist + pp | cbsacode, vcov = "HC1", data = reg_df))
etable(ols_n1, ols_n2, ols_n3, ols_n4, ols_n5, ols_n6)


# only observed null model ties
reg_df <- subset(df, social_null > 0)
table(reg_df$social01, reg_df$social_null01)

summary(ols_n1 <- feols(log_social_null ~ log_dist + pp | cbsacode, vcov = "HC1", data = reg_df))
summary(ols_n2 <- feols(log_social_null ~ log_motorways + log_dist + pp | cbsacode, vcov = "HC1", data = reg_df))
summary(ols_n3 <- feols(log_social_null ~ log_num_barriers + log_dist + pp | cbsacode, vcov = "HC1", data = reg_df))
summary(ols_n4 <- feols(log_social_null ~ income_abs_diff + log_dist + pp | cbsacode, vcov = "HC1", data = reg_df))
summary(ols_n5 <- feols(log_social_null ~ same_race + log_dist + pp | cbsacode, vcov = "HC1", data = reg_df))
summary(ols_n6 <- feols(log_social_null ~ log_motorways + log_num_barriers + income_abs_diff + same_race + log_dist + pp | cbsacode, vcov = "HC1", data = reg_df))
etable(ols_n1, ols_n2, ols_n3, ols_n4, ols_n5, ols_n6)


# all possible connections
reg_df <- subset(df)

summary(ols_n1 <- feols(log_social_null ~ log_dist + pp | cbsacode, vcov = "HC1", data = reg_df))
summary(ols_n2 <- feols(log_social_null ~ log_motorways + log_dist + pp | cbsacode, vcov = "HC1", data = reg_df))
summary(ols_n3 <- feols(log_social_null ~ log_num_barriers + log_dist + pp | cbsacode, vcov = "HC1", data = reg_df))
summary(ols_n4 <- feols(log_social_null ~ income_abs_diff + log_dist + pp | cbsacode, vcov = "HC1", data = reg_df))
summary(ols_n5 <- feols(log_social_null ~ same_race + log_dist + pp | cbsacode, vcov = "HC1", data = reg_df))
summary(ols_n6 <- feols(log_social_null ~ log_motorways + log_num_barriers + income_abs_diff + same_race + log_dist + pp | cbsacode, vcov = "HC1", data = reg_df))
etable(ols_n1, ols_n2, ols_n3, ols_n4, ols_n5, ols_n6)







