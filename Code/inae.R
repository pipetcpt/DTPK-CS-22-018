library('psych')
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

#plot 생성

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


## geometric.mean
dapageo <- dapa_NCA %>%
  group_by(Period) %>%
  summarise_at(vars(CMAX, TMAX, LAMZHL, AUCLST, AUCIFO, VZFO), geometric.mean, na.rm = TRUE) %>%
  mutate_at(vars(-Period), round, 2)%>%
  as.data.frame()

write.csv(dapageo, 'dapageo.csv', row.names = F, fileEncoding = 'cp949')



## CV 계산!!

pktable <- dapa_NCA %>% 
 # select(1:2, "AUC~last~(hr*ng/mL)"=AUCLST, "C~max~(ng/mL)"=CMAX, 
#         "T~max~(hr)"=TMAX, "t~1/2~(hr)"=LAMZHL, "CL/F(L/hr)"=CLFO, "V~d~/F(L)"=VZFO) %>% 
  gather(param, value,CMAX:VZFO) %>% 
  na.omit() %>% #결측치제거
  group_by(Period, param) %>% 
  summarise_at(vars(value), lst(mean, sd, median, min, max)) %>% 
  ungroup() %>% 
  mutate_at(-(1:2), round, 2) %>% 
  # mutate("Mean±SD" = ifelse(param=="Tmax(hr)", 
  #                       sprintf("%0.2f(%0.2f-%0.2f)", median, min, max),
  #                       sprintf("%0.2f±%0.2f", mean, sd))) %>% 
  select(1:2, mean, sd, Median=median, Min=min, Max=max) %>% 
  mutate(CV=round(sd/mean*100,2)) %>%
  mutate(param = factor(param, levels = c("CMAX","TMAX","LAMZHL","AUCLST","AUCIFO","CLFO","VZFO"))) %>%
  arrange(Period,param)
#  gather(stat, value, 3:7) %>% 
#  unite(trtstat, trt, stat, sep="_") %>%
#  mutate(trtstat=factor(trtstat, level=c("R_mean", "R_sd", "R_Median","R_Min","R_Max",
#                                         "T_mean", "T_sd", "T_Median","T_Min","T_Max" ))) %>% 
#  spread(trtstat, value) %>% 
#  split(.$substance)

pktable[,2]<-pktable[,2]%>%factor(levels = c('CMAX','TMAX','LAMZHL','AUCLST','AUCIFO','CLFO','VZFO'))

CV(mean, sd)
#도저히 모르겠당...CV 구하기!!!!!!!!
# sd(dapa_NCA$CMAX) / mean(dapa_NCA%CMAX) * 100
### 일단 1기와 2기를 나눠보자..
dapa_group1 <- dapa_NCA %>%
  filter(Period %in% "1기")

dapa_group2 <- dapa_NCA %>%
  filter(Period %in% "2기")

#CMAX_CV구하기
dapacv1 <- sd(dapa_group1$CMAX) / mean(dapa_group1$CMAX) * 100
dapacv2 <- sd(dapa_group2$CMAX) / mean(dapa_group2$CMAX) * 100

# 아니면
dapacmax1 <- CV(mean=mean(dapa_group1$CMAX), sd=sd(dapa_group1$CMAX))
dapacmax2 <- CV(mean=mean(dapa_group2$CMAX), sd=sd(dapa_group2$CMAX))

table4.1 <- dapa_NCA %>% 
  group_by(Period) %>% 
  summarise(N=n(), Mean=mean(conc), SD= sd(conc), Median=median(conc), Min=min(conc), Max=max(conc)) %>% 
  ungroup() %>% 
  mutate("CV(%)" = SD/Mean*100) %>% 
  mutate_at(-(1:2), round, 2) %>% 
  mutate_at(-(1:2), ~sprintf("%0.2f", .x)) %>% 
  unite("Mean ± SD", Mean, SD, sep = " ± ") %>% 
  select(1:6, 10, everything())





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

#plot 생성
head(metb)

metb %>%
  mutate(Conc = ifelse(Time == 0, 0, Conc)) %>%
  group_by(Period, Time) %>%
  summarise(mean = mean(Conc, na.rm = T), sd = sd(Conc, na.rm = T)) %>%
  ungroup() %>%
  ggplot() +
  geom_line(aes(x = Time, y = mean, col = Period)) +
  geom_point(aes(x = Time, y = mean, col = Period)) + 
  geom_errorbar(aes(x = Time, ymax = mean + sd, ymin = mean, col = Period)) +
  theme_bw() +
  labs(y = "Plasma concentration of Metformin (Cohort B)(ng/mL)")


# Real-time 적용

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


#plot 생성


head(metc)

metc %>%
  mutate(Conc = ifelse(Time == 0, 0, Conc)) %>%
  group_by(Period, Time) %>%
  summarise(mean = mean(Conc, na.rm = T), sd = sd(Conc, na.rm = T)) %>%
  ungroup() %>%
  ggplot() +
  geom_line(aes(x = Time, y = mean, col = Period)) +
  geom_point(aes(x = Time, y = mean, col = Period)) + 
  geom_errorbar(aes(x = Time, ymax = mean + sd, ymin = mean, col = Period)) +
  theme_bw() +
  labs(y = "Plasma concentration of Metformin (Cohort c)(ng/mL)")


# Real-time 적용

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



#plot 생성

head(valsaa)

valsaa %>%
  mutate(Conc = ifelse(Time == 0, 0, Conc)) %>%
  group_by(Period, Time) %>%
  summarise(mean = mean(Conc, na.rm = T), sd = sd(Conc, na.rm = T)) %>%
  ungroup() %>%
  ggplot() +
  geom_line(aes(x = Time, y = mean, col = Period)) +
  geom_point(aes(x = Time, y = mean, col = Period)) + 
  geom_errorbar(aes(x = Time, ymax = mean + sd, ymin = mean, col = Period)) +
  theme_bw() +
  labs(y = "Plasma concentration of Valsartan (Cohort a)(ng/mL)")


#  Real-time 적용 

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

#plot 생성

head(valsac)

valsac %>%
  mutate(Conc = ifelse(Time == 0, 0, Conc)) %>%
  group_by(Period, Time) %>%
  summarise(mean = mean(Conc, na.rm = T), sd = sd(Conc, na.rm = T)) %>%
  ungroup() %>%
  ggplot() +
  geom_line(aes(x = Time, y = mean, col = Period)) +
  geom_point(aes(x = Time, y = mean, col = Period)) + 
  geom_errorbar(aes(x = Time, ymax = mean + sd, ymin = mean, col = Period)) +
  theme_bw() +
  labs(y = "Plasma concentration of Valsartan (Cohort C)(ng/mL)")

# Real-time 

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
