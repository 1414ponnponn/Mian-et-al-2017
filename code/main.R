#####
# set up
library(tidyverse)
library(haven)
library(vars)
library(psych)
library(stargazer)
library(patchwork)
library(fixest)
library(modelsummary)
library(sandwich)
raw_data1 <- haven::read_dta("data/masterdata_annual.dta")
data2 <- read.csv("data/VAR3_levels_fixedeffects_IVprctiles.csv", header = FALSE)
source("code/functions.R")

#####
# table 1

data1 <- raw_data1 %>% 
  tidyr::drop_na(c(lnGDP_Fd3, L1D3HHD_GDP, L1D3NFD_GDP)) %>% 
  dplyr::select(
    c,
    lnGDP_d1,                #\Delta y                            : log real GDP
    lnGDP_Fd3,               #\Delta_{3}y                         : 
    D1PD_GDP,                #\Delta d^{Private}                  : private nonfinancial debt to GDP
    D3PD_GDP,                #\Delta_{3}d^{Private}               : 
    D1HHD_GDP,               #\Delta d^{HH}                       : household debt to GDP
    D3HHD_GDP,               #\Delta_{3}d^{HH}                    : 
    D1NFD_GDP,               #\Delta d^{F}                        : nonfinancial firm debt to GDP
    D3NFD_GDP,               #\Delta_{3}d^{F}                     : 
    L1D3GD_GDP_Sol2,         #\Delta_{3}d^{Gov}                   : government debt to GDP
    L1D3NFL_GDP,             #\Delta_{3}d^{Netforeign}            : the change in net foreign liabilities (sum of current account deficits to GDP)
    lnhhC_realLCU_d1,        #\Delta c                            : log real consumption
    lnhhCdur_realLCU_Fd1,    #\Delta c^{dur}                      : log real durable consumption
    lnhhCnondur_realLCU_Fd1, #\Delta c^{nondur}                   : log real nondurable consumption
    hhC_GDP_d1,              #\Delta\frac{C}{Y}                   : consumption to GDP
    lnI_realLCU_d1,          #\Delta i                            : log real investment
    lnGC_realLCU_d1,         #\Delta g                            : log real government consumption
    lnX_ifs_d1,              #\Delta x                            : log nominal exports
    lnIm_ifs_d1,             #\Delta m                            : log nominal imports
    NX_ifs_GDP_d1,           #\Delta\frac{NX}{Y}                  : net exports to GDP
    CA_GDP_d1,               #\Delta\frac{CA}{Y}                  : current account to GDP
    s_Xc_d1,                 #\Delta s^{XC}                       : the share of consumption exports to total exports
    s_Imc_d1,                #\Delta s^{MC}                       : the share of consumption imports to total imports
    D1lnREER_BIS_m12,        #\Delta reer                         : log real effective exchange rate
    unrate_Fd1,              #\Delta u                            : the unemployment rate
    unrate_Fd3,              #\Delta_{3}u                         : 
    Fd3_lnGDP_WEO_F,         #\Delta_{3}y_{t+3|t}^{WEO}           : the IMF Fall World Economic Outlook time t forecast of growth from t to t+3
    Fd3_WEOerror_F,          #\Delta_{3}(y_{t+3}-y_{t+3|t}^{WEO}) : 
    L1D3lnRHPIq4,            #\Delta_{3}\ln(P^{Housing})          : the log real house price index
    spr_r_q4,                #spr^{real}                          : the real 10-year government bond yield spread with respect to the United States
    ih_ig10,                 #spr^{MS}                            : mortgage-sovereign spread
    ic_ig10, 
    icorp_ig10,              #spr^{corp}                          : the corporate-sovereign spread
    L1toL3hys,
  )

table1 <- psych::describe(data1[,c(2:31, 33)]) %>% 
  data.frame() %>% 
  dplyr::select(n, mean, median, sd) %>% 
  dplyr::rename(N = n,
                Mean = mean,
                Median = median,
                `Std. dev.` = sd) 
sd_deltay <- table1$`Std. dev.`[1] 
table1 <- table1 %>% 
  dplyr::mutate(`Std. dev./Std. dev. (Delta y)` = `Std. dev.` / sd_deltay)


#####
# figure 1

p <- 5
lag_max <-  10
index_y <- c(7, 8, 9)
index_x <- 10:ncol(data2)
nboot <- 2000
clevel <- 95

var_model <- main_VAR(data2, index_y, index_x)
var_bias_corrected <- correct_bias(var_model)
var_bc <- chol_VAR(var_bias_corrected)
var_bs_bc <- wild_bootstrap(var_bc, nboot = nboot, clevel = clevel)

g1 <- plot_irf(var_bc$cholIRF[,1,1], var_bs_bc$CholirsH1[,1,1], var_bs_bc$CholirsL1[,1,1]) +
  ggplot2::ylim(c(-0.2, 1.8)) +
  ggplot2::labs(title = expression(paste(d^{HH}, "->", d^{HH})))

g2 <- plot_irf(var_bc$cholIRF[,3,1], var_bs_bc$CholirsH1[,3,1], var_bs_bc$CholirsL1[,3,1]) +
  ggplot2::ylim(c(-0.8, 0.6)) +
  ggplot2::labs(title = expression(paste(d^{HH}, "->", y)))

g3 <- plot_irf(var_bc$cholIRF[,3,2], var_bs_bc$CholirsH1[,3,2], var_bs_bc$CholirsL1[,3,2]) +
  ggplot2::ylim(c(-0.5, 0.4)) +
  ggplot2::labs(title = expression(paste(d^{F}, "->", y)))

ggsave("figure/Figure1_1.png", g1, width = 5, height = 6)
ggsave("figure/Figure1_2.png", g2, width = 5, height = 6)
ggsave("figure/Figure1_3.png", g3, width = 5, height = 6)

#####
# table 2
result2_1 <- regress_lagged_y(raw_data1, 0)
result2_2 <- regress_lagged_y(raw_data1, 1)
result2_3 <- regress_lagged_y(raw_data1, 2)
result2_4 <- regress_lagged_y(raw_data1, 3)
result2_5 <- regress_lagged_y(raw_data1, 4)
result2_6 <- regress_lagged_y(raw_data1, 5)
result2_7 <- regress_lagged_y(raw_data1, 6)

rows <- dplyr::tribble(
  ~term, ~result2_1, ~result2_2, ~result2_3, ~result2_4, ~result2_5, ~result2_6, ~result2_7,
  "Country fixed effects", "✓", "✓", "✓", "✓", "✓", "✓", "✓"
)
attr(rows, "position") <- 5

table2 <- modelsummary:::modelsummary(
  list(
    result2_1, result2_2, result2_3, result2_4, result2_5, result2_6, result2_7
  ),
  stars = TRUE,
  gof_map = c("r2.within", "nobs", "FE:c", "FE:year"),
  add_rows = rows,
  vcov = ~ year + CountryCode
)

#####
# table 3
result3_1 <- fixest::feols(lnGDP_Fd3 ~ L1D3PD_GDP | CountryCode, 
                           data = raw_data1 %>% 
                             tidyr::drop_na(L1D3NFD_GDP) %>% 
                             dplyr::filter(mainsmp1 == 1))
result3_2 <- fixest::feols(lnGDP_Fd3 ~ L1D3HHD_GDP | CountryCode, 
                           data = raw_data1 %>% 
                             tidyr::drop_na(L1D3NFD_GDP) %>% 
                             dplyr::filter(mainsmp1 == 1))
result3_3 <- fixest::feols(lnGDP_Fd3 ~ L1D3NFD_GDP | CountryCode, 
                           data = raw_data1 %>% 
                             tidyr::drop_na(L1D3NFD_GDP) %>% 
                             dplyr::filter(mainsmp1 == 1))
result3_4 <- fixest::feols(lnGDP_Fd3 ~ L1D3HHD_GDP + L1D3NFD_GDP | CountryCode, 
                           data = raw_data1 %>% 
                             tidyr::drop_na(L1D3NFD_GDP) %>% 
                             dplyr::filter(mainsmp1 == 1))
result3_5 <- fixest::feols(lnGDP_Fd3 ~ L1lnGDP_d1 + L2lnGDP_d1 + L3lnGDP_d1 + L1D3HHD_GDP + L1D3NFD_GDP |
                             CountryCode, 
                           data = raw_data1 %>% 
                             tidyr::drop_na(L1D3NFD_GDP) %>% 
                             dplyr::filter(mainsmp1 == 1))
result3_6 <- fixest::feols(lnGDP_Fd3 ~ L1lnGDP_d1 + L2lnGDP_d1 + L3lnGDP_d1 + L1D3HHD_GDP + L1D3NFD_GDP + L1D3GD_GDP_Sol2 |
                             CountryCode, 
                           data = raw_data1 %>% 
                             tidyr::drop_na(L1D3NFD_GDP) %>% 
                             dplyr::filter(mainsmp1 == 1))
result3_7 <- fixest::feols(lnGDP_Fd3 ~ L1lnGDP_d1 + L2lnGDP_d1 + L3lnGDP_d1 + L1D3HHD_GDP + L1D3NFD_GDP + L1D3NFL_GDP |
                             CountryCode, 
                           data = raw_data1 %>% 
                             dplyr::filter(mainsmp1 == 1))
result3_8 <- fixest::feols(lnGDP_Fd3 ~  L1lnGDP_d1 + L2lnGDP_d1 + L3lnGDP_d1 + L1D3HHD_GDP + L1D3NFD_GDP + L1D3NFLpos + L1D3NFLposXL1D3HHD_GDP|
                             CountryCode,
                           data = raw_data1 %>% 
                             dplyr::filter(mainsmp1 == 1) %>% 
                             mutate(L1toL3CA = dplyr::lag(CA_USD) + dplyr::lag(CA_USD, 2) + dplyr::lag(CA_USD, 3),
                                    L1D3NFLpos = (L1toL3CA < 0),
                                    L1D3NFLpos = if_else(is.na(L1toL3CA), NA, L1D3NFLpos),
                                    L1D3NFLposXL1D3HHD_GDP = L1D3NFLpos * L1D3HHD_GDP))

rows <- dplyr::tribble(
  ~term, ~result3_1, ~result3_2, ~result3_3, ~result3_4, ~result3_5, ~result3_6, ~result3_7, ~result3_8,
  "Country fixed effects", "✓", "✓", "✓", "✓", "✓", "✓", "✓", "✓",
  "Distributed lag in diff log GDP", "", "", "", "", "✓", "✓", "✓", "✓"
)
attr(rows, "position") <- c(15, 16)

table3 <- modelsummary:::modelsummary(
  list(
    result3_1, result3_2, result3_3, result3_4, result3_5, result3_6, result3_7, result3_8
  ),
  stars = TRUE,
  coef_omit = ".*lnGDP_d1",
  gof_map = c("r2.within", "nobs"),
  add_rows = rows,
  vcov = ~ year + CountryCode
)

#####
# table 4a
result4a_1 <- fixest::feols(lnGDP_Fd3 ~ L1lnGDP_d3 + L1D3HHD_GDP + L1D3NFD_GDP |
                             CountryCode, 
                           data = raw_data1 %>% 
                             tidyr::drop_na(lnGDP_Fd3, L1D3NFD_GDP) %>% 
                             dplyr::filter(mainsmp1 == 1,
                                           year %% 3 == 1))
# result4a_2 <- fixest::feols(lnGDP_Fd3 ~ lnGDP_d3 |
#                              L1D3HHD_GDP + L1D3NFD_GDP ~ L4D3HHD_GDP + L4D3NFD_GDP,
#                            data = raw_data1 %>% 
#                              dplyr::mutate(L4D3HHD_GDP = dplyr::lag(L1D3HHD_GDP, 3),
#                                            L4D3NFD_GDP = dplyr::lag(L1D3NFD_GDP, 3)) %>% 
#                              tidyr::drop_na(lnGDP_Fd3, L1D3NFD_GDP) %>% 
#                              dplyr::filter(mainsmp1 == 1,
#                                            year %% 3 == 1)) #incorrect result
result4a_3 <- fixest::feols(lnGDP_Fd3 ~ L1lnGDP_d1 + L2lnGDP_d1 + L3lnGDP_d1 + L1D3HHD_GDP + L1D3NFD_GDP, 
                           data = raw_data1 %>% 
                             tidyr::drop_na(lnGDP_Fd3, L1D3NFD_GDP) %>% 
                             dplyr::filter(mainsmp1 == 1))
# skip result4a_4 
result4a_5 <- fixest::feols(lnGDP_Fd3 ~ L1lnGDP_d1 + L2lnGDP_d1 + L3lnGDP_d1 + L1D3HHD_GDP + L1D3NFD_GDP + year |
                             CountryCode, 
                           data = raw_data1 %>% 
                             tidyr::drop_na(lnGDP_Fd3, L1D3NFD_GDP) %>% 
                             dplyr::filter(mainsmp1 == 1))
result4a_6 <- fixest::feols(lnGDP_Fd3 ~ L1lnGDP_d1 + L2lnGDP_d1 + L3lnGDP_d1 + L1D3HHD_GDP + L1D3NFD_GDP |
                             year + CountryCode, 
                           data = raw_data1 %>% 
                             tidyr::drop_na(lnGDP_Fd3, L1D3NFD_GDP) %>% 
                             dplyr::filter(mainsmp1 == 1))
result4a_7 <- fixest::feols(lnGDP_Fd3 ~ L1lnGDP_d1 + L2lnGDP_d1 + L3lnGDP_d1 + L1D3HHD_L3GDP + L1D3NFD_L3GDP |
                             CountryCode, 
                           data = raw_data1 %>% 
                             tidyr::drop_na(lnGDP_Fd3, L1D3NFD_GDP) %>% 
                             dplyr::filter(mainsmp1 == 1))

rows <- dplyr::tribble(
  ~term, ~result4a_1, ~result4a_3, ~result4a_5, ~result4a_6, ~result4a_7,
  "Country fixed effects", "✓", "", "✓", "✓", "✓", 
  "Distributed lag in diff log GDP", "✓", "✓", "✓", "✓", "✓", 
  "Year fixes effects", "", "", "", "✓", ""
)
attr(rows, "position") <- c(11:13)

table4a <- modelsummary:::modelsummary(
  list(
    result4a_1, result4a_3, result4a_5, result4a_6, result4a_7
  ),
  stars = TRUE,
  coef_omit = c(1, 4:7),
  gof_map = c("r2.within", "nobs"),
  add_rows = rows,
  vcov = ~ year + CountryCode
)

# table 4b
result4b_1 <- fixest::feols(lnGDP_Fd3 ~ L1lnGDP_d1 + L2lnGDP_d1 + L3lnGDP_d1 + L1D3HHD_GDP + L1D3NFD_GDP |
                              CountryCode, 
                            data = raw_data1 %>% 
                              tidyr::drop_na(lnGDP_Fd3, L1D3NFD_GDP) %>% 
                              dplyr::filter(mainsmp1 == 1,
                                            emergingdummy == 0))
result4b_2 <- fixest::feols(lnGDP_Fd3 ~ L1lnGDP_d1 + L2lnGDP_d1 + L3lnGDP_d1 + L1D3HHD_GDP + L1D3NFD_GDP |
                              CountryCode, 
                            data = raw_data1 %>% 
                              tidyr::drop_na(lnGDP_Fd3, L1D3NFD_GDP) %>% 
                              dplyr::filter(mainsmp1 == 1,
                                            emergingdummy == 1))
result4b_3 <- fixest::feols(lnGDP_Fd3 ~ L1lnGDP_d1 + L2lnGDP_d1 + L3lnGDP_d1 + L1D3HHD_GDP + L1D3NFD_GDP |
                              CountryCode, 
                            data = raw_data1 %>% 
                              tidyr::drop_na(lnGDP_Fd3, L1D3NFD_GDP) %>% 
                              dplyr::filter(mainsmp1 == 1,
                                            year <= 1992))
result4b_4 <- fixest::feols(lnGDP_Fd3 ~ L1lnGDP_d1 + L2lnGDP_d1 + L3lnGDP_d1 + L1D3HHD_GDP + L1D3NFD_GDP |
                              CountryCode, 
                            data = raw_data1 %>% 
                              tidyr::drop_na(lnGDP_Fd3, L1D3NFD_GDP) %>% 
                              dplyr::filter(mainsmp1 == 1,
                                            year <= 2003))
result4b_5 <- fixest::feols(lnGDP_Fd3 ~ L1lnGDP_d1 + L2lnGDP_d1 + L3lnGDP_d1 + L1D3HHD_GDP + L1D3NFD_GDP + year |
                              CountryCode, 
                            data = raw_data1 %>% 
                              tidyr::drop_na(lnGDP_Fd3, L1D3NFD_GDP) %>% 
                              dplyr::filter(mainsmp1 == 1,
                                            year <= 2003))
result4b_6 <- fixest::feols(lnGDP_Fd3 ~ L1lnGDP_d1 + L2lnGDP_d1 + L3lnGDP_d1 + L1D3HHD_GDP + L1D3NFD_GDP |
                              CountryCode + year, 
                            data = raw_data1 %>% 
                              tidyr::drop_na(lnGDP_Fd3, L1D3NFD_GDP) %>% 
                              dplyr::filter(mainsmp1 == 1,
                                            year <= 2003))

rows <- dplyr::tribble(
  ~term, ~result4b_1, ~result4b_2, ~result4b_3, ~result4b_4, ~result4b_5, ~result4b_6,
  "Country fixed effects", "✓", "✓", "✓", "✓", "✓", "✓", 
  "Distributed lag in diff log GDP", "✓", "✓", "✓", "✓", "✓", "✓", 
  "Year fixed effects", "", "", "", "", "", "✓", 
)
attr(rows, "position") <- c(7:9)

table4b <- modelsummary:::modelsummary(
  list(
    result4b_1, result4b_2, result4b_3, result4b_4, result4b_5, result4b_6
  ),
  stars = TRUE,
  coef_omit = ".*lnGDP_d1",
  gof_map = c("r2.within", "nobs"),
  add_rows = rows,
  vcov = ~ year + CountryCode
)

#####  
# table 5
result5_1 <- fixest::feols(hhC_GDP_d1 ~ D1HHD_GDP + D1NFD_GDP |
                             CountryCode,
                           data = raw_data1 %>% 
                             tidyr::drop_na(lnGDP_Fd3, L1D3NFD_GDP) %>% 
                             filter(mainsmp1 == 1))
result5_2 <- fixest::feols(hhCnondur_GDP_d1 ~ D1HHD_GDP + D1NFD_GDP |
                             CountryCode,
                           data = raw_data1 %>% 
                             tidyr::drop_na(lnGDP_Fd3, L1D3NFD_GDP) %>% 
                             filter(mainsmp1 == 1))
result5_3 <- fixest::feols(hhCdur_GDP_d1 ~ D1HHD_GDP + D1NFD_GDP |
                             CountryCode,
                           data = raw_data1 %>% 
                             tidyr::drop_na(lnGDP_Fd3, L1D3NFD_GDP) %>% 
                             filter(mainsmp1 == 1))
result5_4 <- fixest::feols(hhCservices_GDP_d1 ~ D1HHD_GDP + D1NFD_GDP |
                             CountryCode,
                           data = raw_data1 %>% 
                             tidyr::drop_na(lnGDP_Fd3, L1D3NFD_GDP) %>% 
                             filter(mainsmp1 == 1))
result5_5 <- fixest::feols(I_GDP_d1 ~ D1HHD_GDP + D1NFD_GDP |
                             CountryCode,
                           data = raw_data1 %>% 
                             tidyr::drop_na(lnGDP_Fd3, L1D3NFD_GDP) %>% 
                             filter(mainsmp1 == 1))
result5_6 <- fixest::feols(NX_ifs_GDP_d1 ~ D1HHD_GDP + D1NFD_GDP |
                             CountryCode,
                           data = raw_data1 %>% 
                             tidyr::drop_na(lnGDP_Fd3, L1D3NFD_GDP) %>% 
                             filter(mainsmp1 == 1))
result5_7 <- fixest::feols(CA_GDP_d1 ~ D1HHD_GDP + D1NFD_GDP |
                             CountryCode,
                           data = raw_data1 %>% 
                             tidyr::drop_na(lnGDP_Fd3, L1D3NFD_GDP) %>% 
                             filter(mainsmp1 == 1))
result5_8 <- fixest::feols(s_Imc_d1 ~ D1HHD_GDP + D1NFD_GDP |
                             CountryCode,
                           data = raw_data1 %>% 
                             tidyr::drop_na(lnGDP_Fd3, L1D3NFD_GDP) %>% 
                             filter(mainsmp1 == 1))
result5_9 <- fixest::feols(s_Xc_d1 ~ D1HHD_GDP + D1NFD_GDP |
                             CountryCode,
                           data = raw_data1 %>% 
                             tidyr::drop_na(lnGDP_Fd3, L1D3NFD_GDP) %>% 
                             filter(mainsmp1 == 1))

rows <- dplyr::tribble(
  ~term, ~result5_1, ~result5_2, ~result5_3, ~result5_4, ~result5_5, ~result5_6, ~result5_7, ~result5_8, ~result5_9,
  "Country fixed effects", "✓", "✓", "✓", "✓", "✓", "✓", "✓", "✓", "✓"
)
attr(rows, "position") <- 5

table5 <- modelsummary::modelsummary(
  list(
    result5_1, result5_2, result5_3, result5_4, result5_5, result5_6, result5_7, result5_8, result5_9
  ),
  stars = TRUE,
  gof_map = c("r2.within", "nobs"),
  add_rows = rows,
  vcov = ~ year + CountryCode
)
