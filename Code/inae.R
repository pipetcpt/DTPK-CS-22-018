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
datad <- read_excel('Data/PK/Dapagliflozin.xlsx', sheet = 5, skip = 1)
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
dapa_NCA <- tblNCA(dapa, key = c("Period", "ID"), colTime = "Time", colConc = "Conc", dose = 10, adm = "Extravascular", R2ADJ = -1, concUnit = 'mg/mL') %>%
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

fc <- LCMAX ~ Period  # LCMAX


BEdc <- lme(fc, random = ~1|ID, data = dapa_BE_raw)    
cidc <- intervals(BEdc, 0.9)
exp(cidc$fixed["Period2기", ])    ## 90% CI result 

GLM(fc, dapa_BE_raw)$ANOVA     ## Anova result




## Comparative PK
dapa_BE_raw <- dapa_NCA  %>%
  mutate(LCMAX = log10(CMAX), LAUCLST = log10(AUCLST))

fa <- LAUCLST ~ Period  # LAUCLST

BEda <- lme(fa, random = ~1|ID, data = dapa_BE_raw)    
cida <- intervals(BEda, 0.9)
exp(cida$fixed["Period2기", ])    ## 90% CI result 

GLM(fa, dapa_BE_raw)$ANOVA     ## Anova result








## Metformin prep(Cohort B)
datamb <- read_excel('Data/PK/Metformin.xlsx', sheet = 5, skip = 1)
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

fmc <- LCMAX ~ Period  # LCMAX

BEmbc <- lme(fmc, random = ~1|ID, data = metformin_BE_rawb)    
cimbc <- intervals(BEmbc, 0.9)
exp(cimbc$fixed["Period2기", ])    ## 90% CI result 

GLM(fmc, metformin_BE_rawb)$ANOVA     ## Anova result





## Comparative PK
metformin_BE_rawb <- metformin_NCAb  %>%
  mutate(LCMAX = log10(CMAX), LAUCLST = log10(AUCLST))

fma <- LAUCLST ~ Period  # LAUCLST

BEmba <- lme(fma, random = ~1|ID, data = metformin_BE_rawb)    
cimba <- intervals(BEmba, 0.9)
exp(cimba$fixed["Period2기", ])    ## 90% CI result 

GLM(fma, metformin_BE_rawb)$ANOVA     ## Anova result








## Metformin prep(Cohort C)
datamc <- read_excel('Data/PK/Metformin.xlsx', sheet = 6, skip = 1)
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

## Comparative PK(Cmax)
metformin_BE_rawc <- metformin_NCAc  %>%
  mutate(LCMAX = log10(CMAX), LAUCLST = log10(AUCLST))

fmc <- LCMAX ~ Period  # LCMAX 

BEmc <- lme(fmc, random = ~1|ID, data = metformin_BE_rawc)    
cimc <- intervals(BEmc, 0.9)
exp(cimc$fixed["Period2기", ])    ## 90% CI result 

GLM(fmc, metformin_BE_rawc)$ANOVA     ## Anova result

## Comparative PK(AUCLST)
metformin_BE_rawc <- metformin_NCAc  %>%
  mutate(LCMAX = log10(CMAX), LAUCLST = log10(AUCLST))

fma <- LAUCLST ~ Period  # LAUCLST

BEma <- lme(fma, random = ~1|ID, data = metformin_BE_rawc)    
cima <- intervals(BEma, 0.9)
exp(cima$fixed["Period2기", ])    ## 90% CI result 

GLM(fma, metformin_BE_rawc)$ANOVA     ## Anova result


## Metformin NCA (Arm B, C bind)
metformin_NCAbc<-rbind(metformin_NCAb,metformin_NCAc)

## Create summary table 
metformin_NCAbc %>%
  select(-ID) %>%
  tbl_summary(by = Period, 
              type = TMAX ~ "continuous", 
              statistic = c("TMAX") ~ "{median} ({min}- {max})")









## 다시 수정 필요함/ Cmax, AUCLST




## Valsartan prep(Cohort A)
datava <- read_excel('Data/PK/Valsartan.xlsx', sheet = 5, skip = 1)
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
datavc <- read_excel('Data/PK/Valsartan.xlsx', sheet = 6, skip = 1)
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