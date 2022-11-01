library(tidyverse)
library(readxl)
library(zoo)
library(NonCompart)
library(gtsummary)
library(sasLM)
library(nlme)

my_theme <-
        list(
                # round large p-values to two places
                "pkgwide-fn:pvalue_fun" = function(x) style_pvalue(x, digits = 2),
                "pkgwide-fn:prependpvalue_fun" = function(x) style_pvalue(x, digits = 2, prepend_p = TRUE),
                # report median (IQR) and n (percent) as default stats in `tbl_summary()`
                "tbl_summary-str:continuous_stat" = "{mean} ± {sd}",
                "tbl_summary-str:categorical_stat" = "{n} ({p}%)"
        )
set_gtsummary_theme(my_theme)
theme_gtsummary_compact()



## Dapa prep
data <- read_excel('Data/PK/Dapagliflozin.xlsx', sheet = 5, skip = 1)
"%ni%" <- Negate("%in%")

dapa <- data %>%
    mutate(Period = ifelse(is.na(Period), na.locf(Period), Period)) %>%
    rename("Time" = 2) %>%
    gather("ID", "Conc", -c("Period", "Time")) %>%
    mutate(Conc = as.numeric(Conc)) %>%
    filter(ID %ni% c("B140", "B070")) %>%
    as.data.frame()

# DB 수령 후 Real-time 적용 예정

## NCA calculation
dapa_NCA <- tblNCA(dapa, key = c("Period", "ID"), colTime = "Time", colConc = "Conc", dose = 100, adm = "Extravascular", R2ADJ = -1, concUnit = 'ng/mL') %>%
    select(Period, ID, CMAX, TMAX, LAMZHL, AUCLST, AUCIFO, CLFO, VZFO) 

## Create summary table 
dapa_NCA %>%
    select(-ID) %>%
    tbl_summary(by = Period, 
                type = TMAX ~ "continuous", 
                statistic = c("TMAX") ~ "{median} ({min}- {max})")

## Comparative PK
dapa_BE_raw <- dapa_NCA  %>%
    mutate(LCMAX = log10(CMAX), LAUCLST = log10(AUCLST))

f <- LCMAX ~ Period  # LCMAX, LAUCLST
# f <- LAUCLST ~ Period 

BE <- lme(f, random = ~1|ID, data = dapa_BE_raw)    
ci <- intervals(BE, 0.9)
exp(ci$fixed["Period2기", ])    ## 90% CI result 

GLM(f, dapa_BE_raw)$ANOVA     ## Anova result