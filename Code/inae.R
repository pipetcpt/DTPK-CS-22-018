library(ggplot2)
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


setwd("C:/Users/cmc/Documents/GitHub/DTPK-CS-22-018")

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


head(dapa)

dapa %>%
  mutate(Conc = ifelse(Time == 0, 0, Conc)) %>%
  group_by(Period, Time) %>%
  summarise(mean = mean(Conc, na.rm = T), sd = sd(Conc, na.rm = T)) %>%
  ungroup() %>%
  ggplot() +
  geom_line(aes(x = Time, y = mean, col = Period)) +
  geom_point(aes(x = Time, y = mean, col = Period)) + 
  geom_errorbar(aes(x = Time, ymax = mean + sd, ymin = mean, col = Period)) +
  theme_bw() +
  labs(y = "Plasma concentration of dapagliflozin(ng/mL)")




# Real-time 적용
db <- list.files('Data/DB', pattern = "xlsx", full.names = T)

db_rn <- read_excel(db, sheet = "RN")
db_pk <- read_excel(db, sheet = "PB")

db_rn


view(dapa)

db_pk_tidy <- db_pk %>%
  mutate(PBDETM = ifelse(SEQ == 1, 0, PBDETM), PBDETM = as.numeric(PBDETM)) %>%
  left_join(db_rn %>% select(SUBJID, RNNO), by = "SUBJID") %>%
  rename(ID = RNNO) %>%
  mutate(Period = ifelse(VISIT %in% c(3, 4), "1기", "2기")) %>%
  mutate(Time = parse_number(PBNT)) 

new_dapa <- dapa %>%
  left_join(db_pk_tidy %>% select(ID, Time, PBDETM, Period), by = c("ID", "Period", "Time")) %>%
  mutate(RTime = Time + PBDETM/60) %>%
  arrange(ID)




## NCA calculation
dapa_NCA <- tblNCA(new_dapa, key = c("Period", "ID"), colTime = "RTime", colConc = "Conc", dose = 10, adm = "Extravascular", R2ADJ = -1, concUnit = 'ng/mL') %>%
  select(Period, ID, CMAX, TMAX, LAMZHL, AUCLST, AUCIFO, CLFO, VZFO) 

## Create summary table 
dapa_NCA %>%
  select(-ID) %>%
  tbl_summary(by = Period, 
              type = TMAX ~ "continuous", 
              statistic = c("TMAX") ~ "{median} ({min}- {max})")

## Comparative PK(CMAX)
dapa_BE_raw <- dapa_NCA  %>%
  mutate(LCMAX = log(CMAX), LAUCLST = log(AUCLST))

fc <- LCMAX ~ Period  # LCMAX


BEdc <- lme(fc, random = ~1|ID, data = dapa_BE_raw)    
cidc <- intervals(BEdc, 0.9)
exp(cidc$fixed["Period2기", ])  %>% round(4)  ## 90% CI result 

GLM(fc, dapa_BE_raw)$ANOVA     ## Anova result




## Comparative PK(AUCLST)
dapa_BE_raw <- dapa_NCA  %>%
  mutate(LCMAX = log(CMAX), LAUCLST = log(AUCLST))

fa <- LAUCLST ~ Period  # LAUCLST

BEda <- lme(fa, random = ~1|ID, data = dapa_BE_raw)    
cida <- intervals(BEda, 0.9)
exp(cida$fixed["Period2기", ])  %>% round(4)   ## 90% CI result 

GLM(fa, dapa_BE_raw)$ANOVA     ## Anova result








## Metformin prep(Cohort B)
datamb <- read_excel('Data/PK/metformin.xlsx', sheet = 5, skip = 1)
"%ni%" <- Negate("%in%")

metb <- datamb %>%
  mutate(Period = ifelse(is.na(Period), na.locf(Period), Period)) %>%
  rename("Time" = 2) %>%
  gather("ID", "Conc", -c("Period", "Time")) %>%
  mutate(Conc = as.numeric(Conc)) %>%
  filter(ID %ni% c("B140", "B070")) %>%
  as.data.frame()





# DB 수령 후 Real-time 적용 예정

view(metb)

db_pk_tidy <- db_pk %>%
  mutate(PBDETM = ifelse(SEQ == 1, 0, PBDETM), PBDETM = as.numeric(PBDETM)) %>%
  left_join(db_rn %>% select(SUBJID, RNNO), by = "SUBJID") %>%
  rename(ID = RNNO) %>%
  mutate(Period = ifelse(VISIT %in% c(3, 4), "1기", "2기")) %>%
  mutate(Time = parse_number(PBNT)) 

new_metb <- metb %>%
  left_join(db_pk_tidy %>% select(ID, Time, PBDETM, Period), by = c("ID", "Period", "Time")) %>%
  mutate(RTime = Time + PBDETM/60) %>%
  arrange(ID)




## NCA calculation
metb_NCA <- tblNCA(new_metb, key = c("Period", "ID"), colTime = "RTime", colConc = "Conc", dose = 1000, adm = "Extravascular", R2ADJ = -1, concUnit = 'ng/mL') %>%
  select(Period, ID, CMAX, TMAX, LAMZHL, AUCLST, AUCIFO, CLFO, VZFO) 

## Create summary table 
metb_NCA %>%
  select(-ID) %>%
  tbl_summary(by = Period, 
              type = TMAX ~ "continuous", 
              statistic = c("TMAX") ~ "{median} ({min}- {max})")

## Comparative PK(CMAX)
metb_BE_raw <- metb_NCA  %>%
  mutate(LCMAX = log(CMAX), LAUCLST = log(AUCLST))

fmbc <- LCMAX ~ Period  # LCMAX


BEmbc <- lme(fmbc, random = ~1|ID, data = metb_BE_raw)    
cimbc <- intervals(BEmbc, 0.9)
exp(cimbc$fixed["Period2기", ])  %>% round(4)   ## 90% CI result 

GLM(fmbc, metb_BE_raw)$ANOVA     ## Anova result




## Comparative PK(AUCLST)
metb_BE_raw <- metb_NCA  %>%
  mutate(LCMAX = log(CMAX), LAUCLST = log(AUCLST))

fmba <- LAUCLST ~ Period  # LAUCLST

BEmba <- lme(fmba, random = ~1|ID, data = metb_BE_raw)    
cimba <- intervals(BEmba, 0.9)
exp(cimba$fixed["Period2기", ])  %>% round(4)   ## 90% CI result 

GLM(fmba, metb_BE_raw)$ANOVA     ## Anova result














## Metformin prep(Cohort C)
datamc <- read_excel('Data/PK/Metformin.xlsx', sheet = 6, skip = 1)
"%ni%" <- Negate("%in%")

metc <- datamc %>%
  mutate(Period = ifelse(is.na(Period), na.locf(Period), Period)) %>%
  rename("Time" = 2) %>%
  gather("ID", "Conc", -c("Period", "Time")) %>%
  mutate(Conc = as.numeric(Conc)) %>%
  filter(ID %ni% c("C040", "C190", "C230")) %>%
  as.data.frame()

# DB 수령 후 Real-time 적용 예정

view(metc)

db_pk_tidy <- db_pk %>%
  mutate(PBDETM = ifelse(SEQ == 1, 0, PBDETM), PBDETM = as.numeric(PBDETM)) %>%
  left_join(db_rn %>% select(SUBJID, RNNO), by = "SUBJID") %>%
  rename(ID = RNNO) %>%
  mutate(Period = ifelse(VISIT %in% c(3, 4), "1기", "2기")) %>%
  mutate(Time = parse_number(PBNT)) 

new_metc <- metc %>%
  left_join(db_pk_tidy %>% select(ID, Time, PBDETM, Period), by = c("ID", "Period", "Time")) %>%
  mutate(RTime = Time + PBDETM/60) %>%
  arrange(ID)




## NCA calculation
metc_NCA <- tblNCA(new_metc, key = c("Period", "ID"), colTime = "RTime", colConc = "Conc", dose = 1000, adm = "Extravascular", R2ADJ = -1, concUnit = 'ng/mL') %>%
  select(Period, ID, CMAX, TMAX, LAMZHL, AUCLST, AUCIFO, CLFO, VZFO) 

## Create summary table 
metc_NCA %>%
  select(-ID) %>%
  tbl_summary(by = Period, 
              type = TMAX ~ "continuous", 
              statistic = c("TMAX") ~ "{median} ({min}- {max})")

## Comparative PK(CMAX)
metc_BE_raw <- metc_NCA  %>%
  mutate(LCMAX = log(CMAX), LAUCLST = log(AUCLST))

fmcc <- LCMAX ~ Period  # LCMAX


BEmcc <- lme(fmcc, random = ~1|ID, data = metc_BE_raw)    
cimcc <- intervals(BEmcc, 0.9)
exp(cimcc$fixed["Period2기", ])   %>% round(4)  ## 90% CI result 

GLM(fmcc, metc_BE_raw)$ANOVA     ## Anova result




## Comparative PK(AUCLST)
metc_BE_raw <- metc_NCA  %>%
  mutate(LCMAX = log(CMAX), LAUCLST = log(AUCLST))

fmca <- LAUCLST ~ Period  # LAUCLST

BEmca <- lme(fmca, random = ~1|ID, data = metc_BE_raw)    
cimca <- intervals(BEmca, 0.9)
exp(cimca$fixed["Period2기", ])  %>% round(4)   ## 90% CI result 

GLM(fmca, metc_BE_raw)$ANOVA     ## Anova result

















## Valsartan prep(Cohort A)
datava <- read_excel('Data/PK/Valsartan.xlsx', sheet = 5, skip = 1)
"%ni%" <- Negate("%in%")

valsaa <- datava %>%
  mutate(Period = ifelse(is.na(Period), na.locf(Period), Period)) %>%
  rename("Time" = 2) %>%
  gather("ID", "Conc", -c("Period", "Time")) %>%
  mutate(Conc = as.numeric(Conc)) %>%
  filter(ID %ni% c("A030", "A120", "A210")) %>%
  as.data.frame()

# DB 수령 후 Real-time 적용 예정

view(valsaa)

db_pk_tidy <- db_pk %>%
  mutate(PBDETM = ifelse(SEQ == 1, 0, PBDETM), PBDETM = as.numeric(PBDETM)) %>%
  left_join(db_rn %>% select(SUBJID, RNNO), by = "SUBJID") %>%
  rename(ID = RNNO) %>%
  mutate(Period = ifelse(VISIT %in% c(3, 4), "1기", "2기")) %>%
  mutate(Time = parse_number(PBNT)) 

new_valsaa <- valsaa %>%
  left_join(db_pk_tidy %>% select(ID, Time, PBDETM, Period), by = c("ID", "Period", "Time")) %>%
  mutate(RTime = Time + PBDETM/60) %>%
  arrange(ID)




## NCA calculation
valsaa_NCA <- tblNCA(new_valsaa, key = c("Period", "ID"), colTime = "RTime", colConc = "Conc", dose = 160, adm = "Extravascular", R2ADJ = -1, concUnit = 'ng/mL') %>%
  select(Period, ID, CMAX, TMAX, LAMZHL, AUCLST, AUCIFO, CLFO, VZFO) 

## Create summary table 
valsaa_NCA %>%
  select(-ID) %>%
  tbl_summary(by = Period, 
              type = TMAX ~ "continuous", 
              statistic = c("TMAX") ~ "{median} ({min}- {max})")

## Comparative PK(CMAX)
valsaa_BE_raw <- valsaa_NCA  %>%
  mutate(LCMAX = log(CMAX), LAUCLST = log(AUCLST))

fvac <- LCMAX ~ Period  # LCMAX


BEvac <- lme(fvac, random = ~1|ID, data = valsaa_BE_raw)    
civac <- intervals(BEvac, 0.9)
exp(civac$fixed["Period2기", ])   %>% round(4)  ## 90% CI result 

GLM(fvac, valsaa_BE_raw)$ANOVA     ## Anova result


## Comparative PK(AUCLST)
valsaa_BE_raw <- valsaa_NCA  %>%
  mutate(LCMAX = log(CMAX), LAUCLST = log(AUCLST))

fvaa <- LAUCLST ~ Period  # LAUCLST

BEvaa <- lme(fvaa, random = ~1|ID, data = valsaa_BE_raw)    
civaa <- intervals(BEvaa, 0.9)
exp(civaa$fixed["Period2기", ])  %>% round(4)   ## 90% CI result 

GLM(fvaa, valsaa_BE_raw)$ANOVA     ## Anova result











## Valsartan prep(Cohort C)
datavc <- read_excel('Data/PK/Valsartan.xlsx', sheet = 6, skip = 1)
"%ni%" <- Negate("%in%")

valsac <- datavc %>%
  mutate(Period = ifelse(is.na(Period), na.locf(Period), Period)) %>%
  rename("Time" = 2) %>%
  gather("ID", "Conc", -c("Period", "Time")) %>%
  mutate(Conc = as.numeric(Conc)) %>%
  filter(ID %ni% c("C040", "C190", "C230")) %>%
  as.data.frame()

# DB 수령 후 Real-time 적용 예정

view(valsac)

db_pk_tidy <- db_pk %>%
  mutate(PBDETM = ifelse(SEQ == 1, 0, PBDETM), PBDETM = as.numeric(PBDETM)) %>%
  left_join(db_rn %>% select(SUBJID, RNNO), by = "SUBJID") %>%
  rename(ID = RNNO) %>%
  mutate(Period = ifelse(VISIT %in% c(3, 4), "1기", "2기")) %>%
  mutate(Time = parse_number(PBNT)) 

new_valsac <- valsac %>%
  left_join(db_pk_tidy %>% select(ID, Time, PBDETM, Period), by = c("ID", "Period", "Time")) %>%
  mutate(RTime = Time + PBDETM/60) %>%
  arrange(ID)




## NCA calculation
valsac_NCA <- tblNCA(new_valsac, key = c("Period", "ID"), colTime = "RTime", colConc = "Conc", dose = 160, adm = "Extravascular", R2ADJ = -1, concUnit = 'ng/mL') %>%
  select(Period, ID, CMAX, TMAX, LAMZHL, AUCLST, AUCIFO, CLFO, VZFO) 

## Create summary table 
valsac_NCA %>%
  select(-ID) %>%
  tbl_summary(by = Period, 
              type = TMAX ~ "continuous", 
              statistic = c("TMAX") ~ "{median} ({min}- {max})")

## Comparative PK(CMAX)
valsac_BE_raw <- valsac_NCA  %>%
  mutate(LCMAX = log(CMAX), LAUCLST = log(AUCLST))

fvcc <- LCMAX ~ Period  # LCMAX


BEvcc <- lme(fvcc, random = ~1|ID, data = valsac_BE_raw)    
civcc <- intervals(BEvcc, 0.9)
exp(civcc$fixed["Period2기", ])  %>% round(4)   ## 90% CI result 

GLM(fvcc, valsac_BE_raw)$ANOVA     ## Anova result


## Comparative PK(AUCLST)
valsac_BE_raw <- valsac_NCA  %>%
  mutate(LCMAX = log(CMAX), LAUCLST = log(AUCLST))

fvca <- LAUCLST ~ Period  # LAUCLST

BEvca <- lme(fvca, random = ~1|ID, data = valsac_BE_raw)    
civca <- intervals(BEvca, 0.9)
exp(civca$fixed["Period2기", ])  %>% round(4)   ## 90% CI result 

GLM(fvca, valsac_BE_raw)$ANOVA     ## Anova result
