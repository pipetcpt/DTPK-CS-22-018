library(tidyverse)
library(readxl)
library(zoo)
library(NonCompart)
library(gtsummary)
library(sasLM)
library(nlme)

setwd("C:/Users/cmc/Desktop/아주PK")

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
datad <- read_excel('Dapagliflozin.xlsx', sheet = 5, skip = 1)
"%ni%" <- Negate("%in%")

dapa <- datad %>%
  mutate(Period = ifelse(is.na(Period), na.locf(Period), Period)) %>%
  rename("Time" = 2) %>%
  gather("ID", "Conc", -c("Period", "Time")) %>%
  mutate(Conc = as.numeric(Conc)) %>%
  filter(ID %ni% c("B140", "B070")) %>%
  as.data.frame()

# DB 수령 후 Real-time 적용 예정

## NCA calculation
dapa_NCA <- tblNCA(dapa, key = c("Period", "ID"), colTime = "Time", colConc = "Conc", dose = 10, adm = "Extravascular", R2ADJ = -1, concUnit = 'ng/mL') %>%
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

BEd <- lme(f, random = ~1|ID, data = dapa_BE_raw)    
cid <- intervals(BEd, 0.9)
exp(cid$fixed["Period2기", ])    ## 90% CI result 

GLM(f, dapa_BE_raw)$ANOVA     ## Anova result










## Metformin prep(Cohort B)
datamb <- read_excel('Metformin.xlsx', sheet = 5, skip = 1)
"%ni%" <- Negate("%in%")

metforminb <- datamb %>%
  mutate(Period = ifelse(is.na(Period), na.locf(Period), Period)) %>%
  rename("Time" = 2) %>%
  gather("ID", "Conc", -c("Period", "Time")) %>%
  mutate(Conc = as.numeric(Conc)) %>%
  filter(ID %ni% c("B140", "B070")) %>%
  as.data.frame()

# DB 수령 후 Real-time 적용 예정

## NCA calculation
metformin_NCAb <- tblNCA(metforminb, key = c("Period", "ID"), colTime = "Time", colConc = "Conc", dose = 1000, adm = "Extravascular", R2ADJ = -1, concUnit = 'ng/mL') %>%
  select(Period, ID, CMAX, TMAX, LAMZHL, AUCLST, AUCIFO, CLFO, VZFO) 

## Create summary table 
metformin_NCAb %>%
  select(-ID) %>%
  tbl_summary(by = Period, 
              type = TMAX ~ "continuous", 
              statistic = c("TMAX") ~ "{median} ({min}- {max})")

## Comparative PK
metformin_BE_rawb <- metformin_NCAb  %>%
  mutate(LCMAX = log10(CMAX), LAUCLST = log10(AUCLST))

f <- LCMAX ~ Period  # LCMAX, LAUCLST
# f <- LAUCLST ~ Period 

BEmb <- lme(f, random = ~1|ID, data = metformin_BE_rawb)    
cimb <- intervals(BEmb, 0.9)
exp(cimb$fixed["Period2기", ])    ## 90% CI result 

GLM(f, metformin_BE_rawb)$ANOVA     ## Anova result










## Metformin prep(Cohort C)
datamc <- read_excel('Metformin.xlsx', sheet = 6, skip = 1)
"%ni%" <- Negate("%in%")

metforminc <- datamc %>%
  mutate(Period = ifelse(is.na(Period), na.locf(Period), Period)) %>%
  rename("Time" = 2) %>%
  gather("ID", "Conc", -c("Period", "Time")) %>%
  mutate(Conc = as.numeric(Conc)) %>%
  filter(ID %ni% c("C040", "C190", "C230")) %>%
  as.data.frame()

# DB 수령 후 Real-time 적용 예정

## NCA calculation
metformin_NCAc <- tblNCA(metforminc, key = c("Period", "ID"), colTime = "Time", colConc = "Conc", dose = 1000, adm = "Extravascular", R2ADJ = -1, concUnit = 'ng/mL') %>%
  select(Period, ID, CMAX, TMAX, LAMZHL, AUCLST, AUCIFO, CLFO, VZFO) 

## Create summary table 
metformin_NCAc %>%
  select(-ID) %>%
  tbl_summary(by = Period, 
              type = TMAX ~ "continuous", 
              statistic = c("TMAX") ~ "{median} ({min}- {max})")

## Comparative PK
metformin_BE_rawc <- metformin_NCAc  %>%
  mutate(LCMAX = log10(CMAX), LAUCLST = log10(AUCLST))

f <- LCMAX ~ Period  # LCMAX, LAUCLST
# f <- LAUCLST ~ Period 

BEmc <- lme(f, random = ~1|ID, data = metformin_BE_rawc)    
cimc <- intervals(BEmc, 0.9)
exp(cimc$fixed["Period2기", ])    ## 90% CI result 

GLM(f, metformin_BE_rawc)$ANOVA     ## Anova result





## Metformin NCA (Arm B, C bind)
metformin_NCAbc<-rbind(metformin_NCAb,metformin_NCAc)

## Create summary table 
metformin_NCAbc %>%
  select(-ID) %>%
  tbl_summary(by = Period, 
              type = TMAX ~ "continuous", 
              statistic = c("TMAX") ~ "{median} ({min}- {max})")









## Valsartan prep(Cohort A)
datava <- read_excel('Valsartan.xlsx', sheet = 5, skip = 1)
"%ni%" <- Negate("%in%")

valsartana <- datava %>%
  mutate(Period = ifelse(is.na(Period), na.locf(Period), Period)) %>%
  rename("Time" = 2) %>%
  gather("ID", "Conc", -c("Period", "Time")) %>%
  mutate(Conc = as.numeric(Conc)) %>%
  filter(ID %ni% c("A030", "A120", "A210")) %>%
  as.data.frame()

# DB 수령 후 Real-time 적용 예정

## NCA calculation
valsartan_NCAa <- tblNCA(valsartana, key = c("Period", "ID"), colTime = "Time", colConc = "Conc", dose = 160, adm = "Extravascular", R2ADJ = -1, concUnit = 'ng/mL') %>%
  select(Period, ID, CMAX, TMAX, LAMZHL, AUCLST, AUCIFO, CLFO, VZFO) 

## Create summary table 
valsartan_NCAa %>%
  select(-ID) %>%
  tbl_summary(by = Period, 
              type = TMAX ~ "continuous", 
              statistic = c("TMAX") ~ "{median} ({min}- {max})")

## Comparative PK
valsartan_BE_rawa <- valsartan_NCAa  %>%
  mutate(LCMAX = log10(CMAX), LAUCLST = log10(AUCLST))

f <- LCMAX ~ Period  # LCMAX, LAUCLST
# f <- LAUCLST ~ Period 

BEva <- lme(f, random = ~1|ID, data = valsartan_BE_rawa)    
civa <- intervals(BEva, 0.9)
exp(civa$fixed["Period2기", ])    ## 90% CI result 

GLM(f, valsartan_BE_rawa)$ANOVA     ## Anova result






## Valsartan prep(Cohort C)
datavc <- read_excel('Valsartan.xlsx', sheet = 6, skip = 1)
"%ni%" <- Negate("%in%")

valsartanc <- datavc %>%
  mutate(Period = ifelse(is.na(Period), na.locf(Period), Period)) %>%
  rename("Time" = 2) %>%
  gather("ID", "Conc", -c("Period", "Time")) %>%
  mutate(Conc = as.numeric(Conc)) %>%
  filter(ID %ni% c("C040", "C190", "C230")) %>%
  as.data.frame()

# DB 수령 후 Real-time 적용 예정

## NCA calculation
valsartan_NCAc <- tblNCA(valsartanc, key = c("Period", "ID"), colTime = "Time", colConc = "Conc", dose = 160, adm = "Extravascular", R2ADJ = -1, concUnit = 'ng/mL') %>%
  select(Period, ID, CMAX, TMAX, LAMZHL, AUCLST, AUCIFO, CLFO, VZFO) 

## Create summary table 
valsartan_NCAc %>%
  select(-ID) %>%
  tbl_summary(by = Period, 
              type = TMAX ~ "continuous", 
              statistic = c("TMAX") ~ "{median} ({min}- {max})")

## Comparative PK
valsartan_BE_rawc <- valsartan_NCAc  %>%
  mutate(LCMAX = log10(CMAX), LAUCLST = log10(AUCLST))

f <- LCMAX ~ Period  # LCMAX, LAUCLST
# f <- LAUCLST ~ Period 

BEvc <- lme(f, random = ~1|ID, data = valsartan_BE_rawc)    
civc <- intervals(BEvc, 0.9)
exp(civc$fixed["Period2기", ])    ## 90% CI result 

GLM(f, valsartan_BE_rawc)$ANOVA     ## Anova result



## Valsartan NCA (Arm A, C bind)
valsartan_NCAac<-rbind(valsartan_NCAa,valsartan_NCAc)

## Create summary table 
valsartan_NCAac %>%
  select(-ID) %>%
  tbl_summary(by = Period, 
              type = TMAX ~ "continuous", 
              statistic = c("TMAX") ~ "{median} ({min}- {max})")