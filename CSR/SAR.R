#library(tidyverse)
#library(lubridate)
#library(kableExtra)
#library(readxl)
#if(!require("qwraps2")){
#  install.packages("qwraps2")
#  library(qwraps2)
#}

#### 0. install packages ####
ipak <-function(pkg){
  new.pkg<-pkg[!(pkg %in% installed.packages()[,"Package"])]
  if(length(new.pkg))
    install.packages(new.pkg,dependencies=TRUE)
  sapply(pkg,require,character.only=TRUE)
}
pkg<-c("tidyverse","lubridate","kableExtra","readxl","qwraps2","DiagrammeR", "dplyr","data.table", "formattable")
ipak(pkg)

setwd("C:/Users/Owner/Documents/GitHub/DTPK-CS-22-018")

# group check! 
IP <- read_excel("Data/DTC21IP070_20221028160911_EXCEL/DataSetExcel_20221028160911/DTC21IP070_DataCenter_DataSet_List_20221028160911.xlsx", sheet = "IP")
group_list <- IP %>% 
  group_by(SUBJID) %>% 
  summarise(group_num = paste(unique(IPARM), collapse = ','))

group_list <- group_list %>% 
  mutate(IP = case_when(group_num == "1,2" ~ "A",
                        group_num == "3,5" ~ "B",
                        group_num == "4,5" ~ "C", 
                        group_num == "1" ~ "A",
                        group_num == "3" ~ "B",
                        group_num == "4" ~ "C"))

table(group_list$group_num, group_list$IP)
# A-1 : 3, A 1-2 : 21 -> 24, 21
# B-3: 1 , B 3-5 : 15 -> 16, 15
# C-4 : 3 , C4-5 : 21 -> 24, 21 


#### Table 1 ####
# scr : Screening summary
# scrf : Reason for screening failure
# dropout : Reason for drop-out

Anal_set<-read_excel("Data/[21HD10503] CSR_Apendix_Template_2022.10.25.xlsx", sheet = "임상대상자_리스트_yj")

# Anal_set + group_num 붙이기 
Anal_set <- left_join(Anal_set, group_list[c(1,3)], by = c("SID" = "SUBJID"))


DS <- read_excel("Data/DTC21IP070_20221028160911_EXCEL/DataSetExcel_20221028160911/DTC21IP070_DataCenter_DataSet_List_20221028160911.xlsx", sheet="DS") 
DS_set<-left_join(DS, Anal_set, by = c("SUBJID"="SID"))%>%
  select(SUBJID, DSYN, DSDTC, DSREAS, DS, SS, PS,RID, IP)

DS_setA<-DS_set %>% filter(IP == "A") #24
DS_setB<-DS_set %>% filter(IP =="B")  #16
DS_setC<-DS_set %>% filter(IP =="C")  #24

DS_dropout_all <- DS_set %>%
  group_by(DS, DSREAS) %>%
  summarize(N = n()) %>%
  ungroup() %>%
  filter(!is.na(DSREAS)) %>%
  full_join(expand.grid(DSREAS = 1:7, DS = c(0, 1)), key = c("DSREAS", "DS")) %>%
  arrange(DSREAS) %>%
  mutate(N = ifelse(is.na(N), 0, N))

TotalScr <- unique(length(DS_set$SUBJID)) #99
TotalRan <- sum(DS_set$DS) #64
TotalAE <- sum(DS_set$SS) #56
TotalPK <- sum(DS_set$PS) # 64

Scr <- c(TotalScr, TotalRan, TotalScr - TotalRan)

ScrF <- DS_dropout_all %>%
  filter(DSREAS %in% c(1,3,7), DS == 0) %>%
  pull(N) # 수기 작성

TotalRanA <- sum(DS_setA$DS)
TotalAEA <- sum(DS_setA$SS)
TotalPKA <- sum(DS_setA$PS)
AlloA <- c(TotalRanA, TotalAEA, TotalRanA - TotalAEA, TotalPKA, TotalRanA - TotalPKA)

TotalRanB <- sum(DS_setB$DS)
TotalAEB <- sum(DS_setB$SS)
TotalPKB <- sum(DS_setB$PS)
AlloB <- c(TotalRanB, TotalAEB, TotalRanB - TotalAEB, TotalPKB, TotalRanB - TotalPKB)

TotalRanC <- sum(DS_setC$DS)
TotalAEC <- sum(DS_setC$SS)
TotalPKC <- sum(DS_setC$PS)
AlloC <- c(TotalRanC, TotalAEC, TotalRanC - TotalAEC, TotalPKC, TotalRanC - TotalPKC)

dropout <- DS_dropout_all %>%
  filter(DSREAS %in% c(2:7), DS == 1) %>%
  pull(N)

table1_names <- read.csv("CSR/Table/Table1_names.csv", stringsAsFactors = F, header = T)
table1 <- data.frame(table1_names, N = c(Scr, ScrF, AlloA, AlloB, AlloC,dropout))

table1
write.csv(table1,"CSR/Table/Table1.csv")

#### Table 2 ####
# code reference : https://cran.r-project.org/web/packages/qwraps2/vignettes/summary-statistics.html

DM <- read_excel("Data/DTC21IP070_20221028160911_EXCEL/DataSetExcel_20221028160911/DTC21IP070_DataCenter_DataSet_List_20221028160911.xlsx",sheet="DM") 
LS <- read_excel("Data/DTC21IP070_20221028160911_EXCEL/DataSetExcel_20221028160911/DTC21IP070_DataCenter_DataSet_List_20221028160911.xlsx",sheet="LS")
VS <- read_excel("Data/DTC21IP070_20221028160911_EXCEL/DataSetExcel_20221028160911/DTC21IP070_DataCenter_DataSet_List_20221028160911.xlsx",sheet="VS")
EG <- read_excel("Data/DTC21IP070_20221028160911_EXCEL/DataSetExcel_20221028160911/DTC21IP070_DataCenter_DataSet_List_20221028160911.xlsx",sheet="EG")

Total <- left_join(DM, LS, by = c("SUBJID","VISIT"))%>%
  left_join(., VS, by = c("SUBJID","VISIT"))%>%
  left_join(., EG, by = c("SUBJID","VISIT"))%>%
  left_join(.,Anal_set,  by = c("SUBJID"="SID"))%>%
  select(SUBJID,VISIT,AGE,SEX,HT,WT,IBW, EXSMKYN,EXAHOLYN,EXCAFFYN, SYSBP,DIABP,PULSE,TEMP, EGVR,EGPR,EGQRSD,EGQT,EGQTC, DS,SS,PS,IP,SUBJID)%>%
  filter(DS==1) 

summary_dm <-
  list("Demographics" = 
         list("Age" = ~paste(round(mean(AGE),2),"\u00B1",round(sd(AGE),2),"(",min(AGE),"-",max(AGE),")"),
              "Height" = ~paste(round(mean(HT),2),"\u00B1",round(sd(HT),2),"(",min(HT),"-",max(HT),")"),
              "Weight" = ~paste(round(mean(WT),2),"\u00B1",round(sd(WT),2),"(",min(WT),"-",max(WT),")"),
              "IBW" = ~paste(round(mean(IBW),2),"\u00B1",round(sd(IBW),2),"(",min(IBW),"-",max(IBW),")")
         ),
       "Vital sign" = 
         list("Systolic blood pressure" = ~paste(round(mean(SYSBP),2),"\u00B1",round(sd(SYSBP),2),"(",min(SYSBP),"-",max(SYSBP),")"),
              "Diastolic blood pressure" = ~paste(round(mean(DIABP),2),"\u00B1",round(sd(DIABP),2),"(",min(DIABP),"-",max(DIABP),")"),
              "Pulse rate" = ~paste(round(mean(PULSE),2),"\u00B1",round(sd(PULSE),2),"(",min(PULSE),"-",max(PULSE),")"),
              "Body temperature" = ~paste(round(mean(TEMP),2),"\u00B1",round(sd(TEMP),2),"(",min(TEMP),"-",max(TEMP),")")
         ),
       "ECG" = 
         list("Ventricular rate" = ~paste(round(mean(EGVR),2),"\u00B1",round(sd(EGVR),2),"(",min(EGVR),"-",max(EGVR),")"),
              "PR interval" = ~paste(round(mean(EGPR),2),"\u00B1",round(sd(EGPR),2),"(",min(EGPR),"-",max(EGPR),")"),
              "QRSD" = ~paste(round(mean(EGQRSD),2),"\u00B1",round(sd(EGQRSD),2),"(",min(EGQRSD),"-",max(EGQRSD),")"),
              "QT" = ~paste(round(mean(EGQT),2),"\u00B1",round(sd(EGQT),2),"(",min(EGQT),"-",max(EGQT),")"),
              "QTc" = ~paste(round(mean(EGQTC),2),"\u00B1",round(sd(EGQTC),2),"(",min(EGQTC),"-",max(EGQTC),")")
         )
  )
whole <-summary_table(Total, summary_dm) 
by_group <- summary_table(group_by(Total,IP),summary_dm)
table2 <-cbind(by_group, whole) 

write.csv(as.data.frame(table2),"CSR/Table/Table2.csv", fileEncoding = "euc-kr")

#### Table 3 ####
MH <- read_excel("Data/DTC21IP070_20221028160911_EXCEL/DataSetExcel_20221028160911/DTC21IP070_DataCenter_DataSet_List_20221028160911.xlsx", sheet="MH") 
MH <- MH %>% 
  group_by(SUBJID) %>%
  summarise(SOC=paste(SOC,collapse=","))

PE <- read_excel("Data/DTC21IP070_20221028160911_EXCEL/DataSetExcel_20221028160911/DTC21IP070_DataCenter_DataSet_List_20221028160911.xlsx", sheet="PE")
PE <- PE %>% 
  filter(PENOR == 2) %>% 
  select(SUBJID, PENOR, PETEST)
# penor -> 0건

Total <- Total %>% 
  mutate(mh=ifelse(SUBJID %in% MH$SUBJID,1,0)) %>%
  mutate(PE=0)%>% #PENOR == 2인 경우가 0건이므로, 임의로 0으로 모두 붙여줌 
  left_join(.,MH, by = c("SUBJID"="SUBJID")) %>%
  select(SUBJID, EXSMKYN, EXAHOLYN,EXCAFFYN, IP, mh,SOC,PE) 

Total$SOC[is.na(Total$SOC)] <-0

summary_mh <- 
  list("A" = 
         list("Number of subjects no medical history" = ~n_perc(mh==0 ,digits = getOption("qwraps2_frmt_digits", 1)),
              "Number of subjects having clinically NOT significant medical history" = ~n_perc(mh==1,digits = getOption("qwraps2_frmt_digits", 1))
         ),
       list("Blood and lymphatic system disorders" = ~n_perc(grepl("Blood and lymphatic system disorders",SOC),digits = getOption("qwraps2_frmt_digits", 1)),
            "Cardiac disorders" = ~n_perc(grepl("Cardiac disorders",SOC)),
            "Congenital, familial and genetic disorders" = ~n_perc(grepl("Congenital, familial and genetic disorders",SOC),digits = getOption("qwraps2_frmt_digits", 1)),
            "Ear and labyrinth disorders" = ~n_perc(grepl("Ear and labyrinth disorders",SOC),digits = getOption("qwraps2_frmt_digits", 1)),
            "Endocrine disorders" = ~n_perc(grepl("Endocrine disorders",SOC),digits = getOption("qwraps2_frmt_digits", 1)),
            "Eye disorders" = ~n_perc(grepl("Eye disorders",SOC),digits = getOption("qwraps2_frmt_digits", 1)),
            "Gastrointestinal disorders" = ~n_perc(grepl("Gastrointestinal disorders",SOC),digits = getOption("qwraps2_frmt_digits", 1)),
            "General disorders and administration site conditions" = ~n_perc(grepl("General disorders and administration site conditions",SOC),digits = getOption("qwraps2_frmt_digits", 1)),
            "Hepatobiliary disorders" = ~n_perc(grepl("Hepatobiliary disorders",SOC),digits = getOption("qwraps2_frmt_digits", 1)),
            "Infections and infestation" = ~n_perc(grepl("Infections and infestations",SOC),digits = getOption("qwraps2_frmt_digits", 1)),
            "Immune system disorders" = ~n_perc(grepl("Immune system disorders",SOC), digits = getOption("qwraps2_frmt_digits", 1)),
            "Injury, poisoning and procedural complications" = ~n_perc(grepl("Injury, poisoning and procedural complications",SOC),digits = getOption("qwraps2_frmt_digits", 1)),
            "Investigations" = ~n_perc(grepl("Investigations",SOC),digits = getOption("qwraps2_frmt_digits", 1)),
            "Metabolism and nutrition disorders" = ~n_perc(grepl("Metabolism and nutrition disorders",SOC),digits = getOption("qwraps2_frmt_digits", 1)),
            "Musculoskeletal and connective tissue disorders" = ~n_perc(grepl("Musculoskeletal and connective tissue disorders",SOC),digits = getOption("qwraps2_frmt_digits", 1)),
            "Neoplasms benign, malignant and unspecified (incl cysts and polyps)" = ~n_perc(grepl("Neoplasms benign, malignant and unspecified",SOC),digits = getOption("qwraps2_frmt_digits", 1)),
            "Nervous system disorders" = ~n_perc(grepl("Nervous system disorders",SOC),digits = getOption("qwraps2_frmt_digits", 1)),
            "Pregnancy, puerperium and perinatal conditions" = ~n_perc(grepl("Pregnancy, puerperium and perinatal conditions",SOC),digits = getOption("qwraps2_frmt_digits", 1)),
            "Product issues" = ~n_perc(grepl("Product issues",SOC),digits = getOption("qwraps2_frmt_digits", 1)),
            "Psychiatric disorders" = ~n_perc(grepl("Psychiatric disorders",SOC),digits = getOption("qwraps2_frmt_digits", 1)),
            "Renal and urinary disorders" = ~n_perc(grepl("Renal and urinary disorders",SOC),digits = getOption("qwraps2_frmt_digits", 1)),
            "Reproductive system and breast disorders" = ~n_perc(grepl("Reproductive system and breast disorders",SOC),digits = getOption("qwraps2_frmt_digits", 1)),
            "Respiratory, thoracic and mediastinal disorders" = ~n_perc(grepl("Respiratory, thoracic and mediastinal disorders",SOC),digits = getOption("qwraps2_frmt_digits", 1)),
            "Skin and subcutaneous tissue disorders" = ~n_perc(grepl("Skin and subcutaneous tissue disorders",SOC),digits = getOption("qwraps2_frmt_digits", 1)),
            "Social circumstances" = ~n_perc(grepl("Social circumstances",SOC),digits = getOption("qwraps2_frmt_digits", 1)),
            "Surgical and medical procedures" = ~n_perc(grepl("Surgical and medical procedures",SOC),digits = getOption("qwraps2_frmt_digits", 1)),
            "Vascular disorders" = ~n_perc(grepl("Vascular disorders",SOC),digits = getOption("qwraps2_frmt_digits", 1))
       ),
       list( "General" = ~n_perc(PE=="General",digits = getOption("qwraps2_frmt_digits", 1)),
             "Nutrition" = ~n_perc(PE=="Nutrition",digits = getOption("qwraps2_frmt_digits", 1)),
             "Integumentary system (skin/mucosa)" = ~n_perc(PE=="Integumentary system (skin/mucosa)",digits = getOption("qwraps2_frmt_digits", 1)),
             "Ophthalmologic system (eye, excluding decrease visual acuity)" = ~n_perc(PE=="Ophthalmologic system (eye, excluding decrease visual acuity)",digits = getOption("qwraps2_frmt_digits", 1)),
             "Ear, nose, & throat" = ~n_perc(PE=="Ear, nose, & throat",digits = getOption("qwraps2_frmt_digits", 1)),
             "Thyroid" = ~n_perc(PE=="Thyroid",digits = getOption("qwraps2_frmt_digits", 1)),
             "Respiratory system" = ~n_perc(PE=="Respiratory system",digits = getOption("qwraps2_frmt_digits", 1)),
             "Cardiovascular system" = ~n_perc(PE=="Cardiovascular system",digits = getOption("qwraps2_frmt_digits", 1)),
             "Abdomen" = ~n_perc(PE=="Abdomen",digits = getOption("qwraps2_frmt_digits", 1)),
             "Kidney / Genitourinary system" = ~n_perc(PE=="Kidney / Genitourinary system",digits = getOption("qwraps2_frmt_digits", 1)),
             "Neuropsychiatry" = ~n_perc(PE=="Neuropsychiatry",digits = getOption("qwraps2_frmt_digits", 1)),
             "Vertebra / Limb / Any malignancies" = ~n_perc(PE=="Vertebra / Limb / Any malignancies",digits = getOption("qwraps2_frmt_digits", 1)),
             "Peripheral blood supply" = ~n_perc(PE=="Peripheral blood supply",digits = getOption("qwraps2_frmt_digits", 1)),
             "Lymphatics" = ~n_perc(PE=="Lymphatics",digits = getOption("qwraps2_frmt_digits", 1)),
             "Others" = ~n_perc(PE=="Others",digits = getOption("qwraps2_frmt_digits", 1))
       ),
       list("Number of smokers (NOT exceeding the amount indicated in exclusion criteria)" = ~n_perc(EXSMKYN ==1,digits = getOption("qwraps2_frmt_digits", 1)),
            "Number of non-smokiers" = ~n_perc(EXSMKYN == 2,digits = getOption("qwraps2_frmt_digits", 1))
       ),
       list("Number of subjects consuming alcohol (NOT exceeding the amount indicated in exclusion criteria)" = ~n_perc(EXAHOLYN ==1,digits = getOption("qwraps2_frmt_digits", 1)),
            "Number of subjects NOT consuming alcohol" = ~n_perc(EXAHOLYN==2,digits = getOption("qwraps2_frmt_digits", 1))
       ),
       list("Number of subjects consuming caffeine (NOT exceeding the amount indicated in exclusion criteria)" = ~n_perc(EXCAFFYN ==1,digits = getOption("qwraps2_frmt_digits", 1)),
            "Number of subjects NOT consuming caffeine" = ~n_perc(EXCAFFYN==2,digits = getOption("qwraps2_frmt_digits", 1))
       )
  )
whole2 <-summary_table(Total,summary_mh)
by_group2 <- summary_table(group_by(Total,IP),summary_mh)
table3 <-cbind(by_group2, whole2)
        
write.csv(as.data.frame(table3),"CSR/Table/Table3.csv")

#### AE #####
AE <- read_excel("Data/DTC21IP070_20221028160911_EXCEL/DataSetExcel_20221028160911/DTC21IP070_DataCenter_DataSet_List_20221028160911.xlsx", sheet="AE")
# AE <- AE %>% 
#   select(SUBJID, AETYPE, AETERM, AESTDTC, AESTTC, AEENDTC, AEENTC, AESER, AESEV, AEREL, AEOUT, PT, SOC)

AE <- left_join(AE, Anal_set[c(1,2,6)], by = c("SUBJID" = "SID"))
AE$AESTDTC <- ymd(AE$AESTDTC)

# add period date
group <- c("A","B","C")
p1_date <- c(20220727, 20220727, 20220816)
p2_date <- c(20220803, 20220803, 20220823)
df<- data.frame(group, p1_date, p2_date)

df$p1_date <- ymd(df$p1_date)
df$p2_date <- ymd(df$p2_date)
str(df)

AE <- left_join(AE, df, by = c("IP" = "group"))         
AE <- AE %>% 
  mutate(Period_1 = ifelse(AESTDTC >= p1_date & AESTDTC < p2_date, 1,0),
         Period_2 = ifelse(AESTDTC >= p2_date, 1,0),
         Period = ifelse(AESTDTC >= p1_date & AESTDTC < p2_date, "Period_1", "Period_2"))

# 인애쌤 전달 AE 파일
write.csv(AE, "AE.csv", row.names = FALSE, fileEncoding = "euc-kr")
#########################################################################################

# Table 8, Table 11, Table 14
# AE_A <- AE %>%
#   filter(IP == "A")
# AE_B <- AE %>%
#   filter(IP == "B")
# AE_C <- AE %>%
#   filter(IP == "C")
# 
# 
# # 중복제거
# AE_A_filter <- AE_A[!duplicated(AE_A[,c("SUBJID","Period")]),]
# AE_B_filter <- AE_B[!duplicated(AE_B[,c("SUBJID","Period")]),]
# AE_C_filter <- AE_C[!duplicated(AE_C[,c("SUBJID","Period")]),]
# 
# summary_AE1 <-
#   list("A" =
#          list("Number of adverse event" = ~n_perc(!is.na(IP)),
#               "Number of adverse drug reaction" = ~n_perc(AEREL%in% c(2,3,4,5)),
#               "Number of serious adverse event" = ~n_perc(AESER == 1),
#               "Number of serious adverse drug reaction" = ~n_perc(AEREL %in% c(2,3,4,5) & AESER == 1)
#               )
#        )
# 
# t8_whole <- summary_table(AE_A_filter, summary_AE1)
# by_period <- summary_table(group_by(AE_A_filter, Period), summary_AE1)
# table_8 <- cbind(by_period, t8_whole)
# 
# t11_whole <- summary_table(AE_B_filter, summary_AE)
# by_period <- summary_table(group_by(AE_B_filter, Period), summary_AE)
# table_11 <- cbind(by_period, t11_whole)
# 
# t14_whole <- summary_table(AE_C, summary_AE)
# by_period <- summary_table(group_by(AE_C, Period), summary_AE)
# table_14 <- cbind(by_period, t14_whole)
# 
# write.csv(table_8, "CSR/Table/table_8.csv", row.names = TRUE)
# write.csv(table_11, "CSR/Table/table_11.csv", row.names = TRUE)
# write.csv(table_14, "CSR/Table/table_14.csv", row.names = TRUE)


# p-value 



