# =============================================================================
# VAr.R  --  The Impact of Mining Costs on the Cryptocurrency Market
# -----------------------------------------------------------------------------
# Main empirical analysis for the thesis. Estimates the bidirectional, short-
# and long-run relationship between the Bitcoin price and the electricity
# consumed mining it, using a (structural) vector autoregression.
#
# Pipeline:
#   1. Load and merge daily price and energy-consumption series
#   2. Test stationarity (ADF) and cointegration (Johansen)
#   3. Difference the series and pick the VAR lag order
#   4. Granger-causality + ITCV robustness, HAC-robust coefficient tests
#   5. Impulse-response functions (reduced-form and structural SVAR)
#
# NOTE: this is an exploratory script meant to be run interactively in chunks,
#       not strictly top-to-bottom. A few blocks reference objects that are
#       created further down (see flags marked "ORDER:") and a couple of
#       leftover lines reference undefined objects (marked "BUG:"). They are
#       left unchanged to preserve the original analysis; comments flag them.
# =============================================================================

# ---- Packages ---------------------------------------------------------------
library(dynlm)            # dynamic linear models with L() lag operator
library(stargazer)        # publication-ready regression/summary tables
library("readxl")         # read .xlsx input files
library(TimeSeries.OBeu)  # time-series helpers
library("writexl")        # write .xlsx (not strictly required here)
library(SciViews)         # provides ln(); descr() etc.
library(tseries)          # adf.test() (augmented Dickey-Fuller)
library(lmtest)           # coeftest(), Breusch-Godfrey, etc.
library(sjmisc)           # descr() descriptive statistics
library(zoo)              # ordered/irregular time series
library(xts)              # extensible time series
library(outliers)         # outlier tests (Grubbs etc.)
library(AICcmodavg)       # AICc model selection tables
library(vars)             # VAR(), SVAR(), irf(), serial.test(), arch.test()
library(forecast)         # auto.arima()
library(ggplot2)          # plotting
library(patchwork)        # combine ggplots with + and /
library(sandwich)         # vcovHAC() heteroskedasticity/autocorr-robust vcov
library(hrbrthemes)       # theme_ipsum() plot theme
library(urca)             # ca.jo() Johansen cointegration test
library(konfound)         # ITCV (impact threshold of a confounding variable)
library(xtable)           # export tables to LaTeX

# Working directory holding the input data files.
# PORTABILITY: hardcoded absolute path -- change this to your local folder
# (or better, use an RStudio Project / here::here() so the script is portable).
setwd("D:/Downloads")

# ---- 1. Load data -----------------------------------------------------------
# Bitcoin daily price (source: CoinMarketCap export).
Price <- read.csv("BTC-Daily.csv")
Price <- Price[,c(2,7)]                              # keep columns 2 (date) and 7 (close)
colnames(Price)[1] ="Date"
dates <- Price$Date
Price$Date <- as.Date(dates, "%Y-%m-%d %H:%M:%S")    # parse the timestamp to a Date

# Bitcoin energy consumption (source: Cambridge CBECI footprint export).
Energy <- read_xlsx("BTC_Footprints_v1.xlsx")
Energy <- Energy[,c(1,4)]                            # keep date + the BTCENEGUE column
dates <- Energy$Date
Energy$Date <- as.Date(dates, tz = "Hongkong")
Energy$BTCENEGUE <- Energy$BTCENEGUE/1e6             # rescale raw units -> GW used in the thesis

# Alternative input (pre-computed returns) -- kept for reference, not used.
#Price <- read_xlsx("PrRet.xlsx")
#colnames(Price)[1] ="Date"
#dates <- Price$Date
#Price$Date <- as.Date(dates, tz = "Hongkong")
#Energy <- read_xlsx("EnRet.xlsx")
#dates <- Energy$Date
#Energy$Date <- as.Date(dates, tz = "Hongkong")


# Inner join on Date -> only dates present in BOTH series are kept.
df <- merge(Price, Energy)
colnames(df)[2]='Price'
colnames(df)[3]='Energy'

# Quick raw-series check -- kept for reference.
#pp <- ggplot(df, aes(x=Date, y = close)) + geom_line()
#ep <- ggplot(df, aes(x=Date, y = BTCENEGUE)) + geom_line()
#pp + ep

# Convert each column to a ts object for the time-series routines.
p <- ts(df$Price)
e <- ts(df$Energy)

# ---- 2. Cointegration (Johansen) --------------------------------------------
# Johansen trace test on the two levels series, with a linear trend in the
# cointegration relation. Used to decide VECM vs. (S)VAR.
# NOTE: K = 18 lags here; the thesis text mentions 14 -- worth reconciling so
#       the reported table matches the code exactly.
tsdf <- df[-1]                                       # drop the Date column -> price + energy only
vecm.model <- ca.jo(
  tsdf, ecdet = "trend",
  type  = "trace", K = 18, spec = "transitory", dumvar = NULL)
summary(vecm.model)                                  # trace stats vs critical values (see thesis Table)

# ORDER: this block plots `p`/`e` but feeds them to ggplot as if data frames;
#        kept as in the original (the publication chart is the block below).
pp <- ggplot(p, aes(x=df$Date, y = p)) + geom_line()
ep <- ggplot(e, aes(x=df$Date, y = e)) + geom_line()
pp + ep

# Stationarity of the LEVELS. High p-value for price => non-stationary.
adf.test(p)
adf.test(e)

# ---- 3. Difference to stationarity ------------------------------------------
# First differences (daily changes). The thesis works with these I(0) series.
pd <- diff(p, differences = 1)
ed <- diff(e, differences = 1)

# ORDER: this dual-axis chart uses `dfc`, which is only created later (line ~126).
#        Run it AFTER dfc exists. It reproduces the "Vars.png" figure.
pp <- ggplot(dfc, aes(x=Date)) + geom_line(aes(y = price), size = 0.8, color = priceColor)
ep <- ggplot(dfc, aes(x=Date)) + geom_line(aes(y = energy), size = 0.8, color = EnergyColor)
pp /
ep

#graph

# Value used to transform the data
coeff <- 2                                           # scale factor so price shares the energy axis

# A few constants
EnergyColor <- "firebrick"
  priceColor <- "dodgerblue"

  # Publication chart: energy (left axis) vs price*coeff (right axis = price/coeff).
  ggplot(df, aes(x=Date)) +

    geom_line( aes(y=Energy), size=0.8, color=EnergyColor) +
    geom_line( aes(y=Price * coeff), size=0.8, color=priceColor) +

    scale_y_continuous(

      # Features of the first axis
      name = "Energy consumption (GW)",

      # Add a second axis and specify its features
      sec.axis = sec_axis(~./coeff, name="Price ($)")
    ) +

    theme_ipsum() +

    theme(
      axis.title.y = element_text(color = EnergyColor, size=13),
      axis.title.y.right = element_text(color = priceColor, size=13)
    )

# Stationarity of the FIRST DIFFERENCES -> both stationary at 1% in the thesis.
adf.test(pd)
adf.test(ed)

# Descriptive statistics (feeds the "Descriptive statistics" table).
descr(df)
summary(df)

# FIXED: removed two broken lines `adf.test(pdc)` / `adf.test(edc)` -- the objects
#        `pdc` / `edc` were never defined (leftovers) and would error. The
#        stationarity of the first differences is already tested above via
#        adf.test(pd) / adf.test(ed).

# ---- 4a. Granger causality (necessary condition for two-way link) -----------
# Tested both directions. Low p-values => bidirectional Granger causality.
grangertest(pd ~ ed, order = 1)
grangertest(ed ~ pd, order = 1)

# Reference: best univariate ARIMA for each differenced series.
bestp <- auto.arima(pd)
summary(bestp)
beste <- auto.arima(ed)
summary(beste)

# ---- 4b. Build the modelling data frame -------------------------------------
dfc <- data.frame(price = pd, energy = ed)           # the two stationary series
dfc$Date <- df$Date[-1]                               # FIXED: drop the first date lost to diff(); adapts to any sample length

# ---- 4c. VAR lag selection and estimation -----------------------------------
# FIXED: estimate on an explicit 2-variable frame so the Date column is NOT
#        pulled into the VAR as a spurious third endogenous variable.
dfvar <- dfc[, c("price", "energy")]
VARselect(dfvar, lag.max = 120, type = "const")      # information criteria by lag
varm <- VAR(dfvar, p = 18, type = "const", ic = c("AIC", "HQ", "SC", "FPE"))  # SC picks p = 18
tmp <- summary(varm)
xtable(tmp$varresult$price)                          # export price-equation coefficients to LaTeX
xtable(tmp$varresult$B)

# Residual diagnostics: serial correlation + ARCH effects.
serial.test(varm, type = 'PT.adjusted')
arch.test(varm)

# ---- 4d. ITCV robustness on the simple linear link --------------------------
linm <- lm(energy ~ price, dfc)                      # simple OLS used only for the ITCV test

# Full 18-lag ARDL for the price equation (manual dynlm specification).
plm <- dynlm(price~L(price)+L(price,2)+L(price,3)+L(price,4)+L(price,5)+L(price,6)+L(price,7)+L(price,8)+L(price,9)+L(price,10)+L(price,11)+L(price,12)+L(price,13)+L(price,14)+L(price,15)+L(price,16)+L(price,17)+L(price,18)+L(energy)+L(energy,2)+L(energy,3)+L(energy,4)+L(energy,5)+L(energy,6)+L(energy,7)+L(energy,8)+L(energy,9)+L(energy,10)+L(energy,11)+L(energy,12)+L(energy,13)+L(energy,14)+L(energy,15)+L(energy,16)+L(energy,17)+L(energy,18), dfc)
summary(plm)

# ITCV: how strong a confounder would have to be to overturn the price effect.
konfound_output <- konfound(linm, price, to_return = c("raw_output", "thresh_plot", "corr_plot"))
konfound_output$thresh_plot                          # reproduces "Frank.png"

# HAC-robust coefficient test (the model has both heteroskedasticity & autocorr).
coeftest(plm, vcov = vcovHAC(plm))

# Same 18-lag ARDL for the energy equation.
ple <- dynlm(energy~L(price)+L(price,2)+L(price,3)+L(price,4)+L(price,5)+L(price,6)+L(price,7)+L(price,8)+L(price,9)+L(price,10)+L(price,11)+L(price,12)+L(price,13)+L(price,14)+L(price,15)+L(price,16)+L(price,17)+L(price,18)+L(energy)+L(energy,2)+L(energy,3)+L(energy,4)+L(energy,5)+L(energy,6)+L(energy,7)+L(energy,8)+L(energy,9)+L(energy,10)+L(energy,11)+L(energy,12)+L(energy,13)+L(energy,14)+L(energy,15)+L(energy,16)+L(energy,17)+L(energy,18), dfc)
summary(ple)
coeftest(ple, vcov = vcovHAC(ple))

# Parsimonious 7-lag versions of both equations (re-uses names plm / ple).
plm <- dynlm(price~L(price)+L(price,2)+L(price,3)+L(price,4)+L(price,5)+L(price,6)+L(price,7)+L(energy)+L(energy,2)+L(energy,3)+L(energy,4)+L(energy,5)+L(energy,6)+L(energy,7), dfc)
summary(plm)
coeftest(plm, vcov = vcovHAC(plm))
ple <- dynlm(energy~L(price)+L(price,2)+L(price,3)+L(price,4)+L(price,5)+L(price,6)+L(price,7)+L(energy)+L(energy,2)+L(energy,3)+L(energy,4)+L(energy,5)+L(energy,6)+L(energy,7), dfc)
summary(ple)
coeftest(ple, vcov = vcovHAC(ple))

# ---- 5a. Reduced-form impulse-response functions ----------------------------
# n.ahead = 30 days, 90% CI. FIXED: cumulative now passed as a logical
# (FALSE/TRUE) instead of the strings "false"/"true".
# Non-cumulative IRFs (per-day response):
irfp <- irf(varm, impulse ="energy", n.ahead = 30, response = "price", ci = 0.9, cumulative = FALSE)
plot(irfp)
irfe <- irf(varm, impulse ="price", n.ahead = 30, response = "energy", ci = 0.9, cumulative = FALSE)
plot(irfe)
irfp <- irf(varm, impulse ="price", n.ahead = 30, response = "price", ci = 0.9, cumulative = FALSE)
plot(irfp)
irfe <- irf(varm, impulse ="energy", n.ahead = 30, response = "energy", ci = 0.9, cumulative = FALSE)
plot(irfe)

# Cumulative (accumulated) IRFs -- these are the four figures in the thesis
# (enpr.png, pren.png, prpr.png, enen.png).
irfp <- irf(varm, impulse ="energy", n.ahead = 30, response = "price", ci = 0.9, cumulative = TRUE)
plot(irfp)
irfe <- irf(varm, impulse ="price", n.ahead = 30, response = "energy", ci = 0.9, cumulative = TRUE)
plot(irfe)
irfp <- irf(varm, impulse ="price", n.ahead = 30, response = "price", ci = 0.9, cumulative = TRUE)
plot(irfp)
irfe <- irf(varm, impulse ="energy", n.ahead = 30, response = "energy", ci = 0.9, cumulative = TRUE)
plot(irfe)


# ---- 5b. Structural VAR (SVAR) ----------------------------------------------
#SVAR
# A-model: estimate contemporaneous matrix A. diag = NA (free), and amat[1,2]=NA
# leaves the (1,2) element free while amat[2,1] stays 0 -> the identifying
# restriction that an energy shock has NO contemporaneous effect on price.
amat <- diag(2)
diag(amat) <- NA
amat[1,2] <- NA


final <- SVAR(varm, Amat = amat)                     # structural model
# Accumulated structural IRFs (30- and 60-day horizons).
sirfp <- irf(final, impulse ="energy", n.ahead = 30, response = "price", ci = 0.9, cumulative = TRUE)
plot(sirfp)
sirfe <- irf(final, impulse ="price", n.ahead = 30, response = "energy", ci = 0.9, cumulative = TRUE)
plot(sirfe)
sirfp <- irf(final, impulse ="price", n.ahead = 60, response = "price", ci = 0.9, cumulative = TRUE)
plot(sirfp)
sirfe <- irf(final, impulse ="energy", n.ahead = 60, response = "energy", ci = 0.9, cumulative = TRUE)
plot(sirfe)
summary(final)
