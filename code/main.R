library(tidyverse)
library(dplyr)
library(haven)
library(psych)
library(fastDummies)
library(vars)
library(lpirfs)
library(xtable)
library(stargazer)

raw_data1 <- haven::read_dta("data/masterdata_annual.dta")

data1 <- raw_data1 %>% 
  tidyr::drop_na(c(lnGDP_Fd3, L1D3HHD_GDP, L1D3NFD_GDP)) %>% 
  dplyr::select(
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
    L1toL3hys
  )


data2 <- read.csv("data/VAR3_levels_fixedeffects_IVprctiles.csv")

p <- 5
lag_max <-  10
index_y <- c(7, 8, 9)
index_x <- 10:ncol(data2)

amat <- diag(3)
amat[2,1] <- NA
amat[3,1] <- NA
amat[3,2] <- NA

var_model <- vars::VAR(data2[index_y], p = p, type = "const", exogen = data2[index_x], lag.max = lag_max)
svar_model <- vars::SVAR(var_model, estmethod = "scoring", Amat = amat, Bmat = NULL, hessian = TRUE, method="BFGS")
irf1 <- vars::irf(svar_model, boot = TRUE)
fig1 <- plot(irf1)
