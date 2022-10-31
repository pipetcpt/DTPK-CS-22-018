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
pkg<-c("tidyverse","lubridate","kableExtra","readxl","qwraps2","DiagrammeR")
ipak(pkg)

#### Table 1 ####
Anal_set<-read_excel("data/DKF-5122_appendix_Analysis_variable.xlsx",sheet=2)
DS <- read_excel("data/exceldata/DWPDN11_P01R_DataCenter_DataSet_List_20220817084019.xlsx",sheet="DS") 
DS_set<-left_join(DS, Anal_set, by = c("SUBJID"="SID"))%>%
  select(SUBJID, DSYN, DSDTC, DSREAS, DS, SS, PS)

DS_dropout_all <- DS_set %>%
  group_by(DS, DSREAS) %>%
  summarize(N = n()) %>%
  ungroup() %>%
  filter(!is.na(DSREAS)) %>%
  full_join(expand.grid(DSREAS = 1:7, DS = c(0, 1)), key = c("DSREAS", "DS")) %>%
  arrange(DSREAS) %>%
  mutate(N = ifelse(is.na(N), 0, N))

TotalScr <- unique(length(DS_set$SUBJID)) #18
TotalRan <- sum(DS_set$DS) #8
TotalAE <- sum(DS_set$SS)  #8
TotalPK <- sum(DS_set$PS)  #6

Scr <- c(TotalScr, TotalRan, TotalScr - TotalRan) #18, 8, 10

ScrF <- DS_dropout_all %>%
  filter(DSREAS %in% c(1,2,7), DS == 0) %>% #10,0,0
  pull(N)

Allo <- c(TotalRan, TotalAE, TotalRan - TotalAE, TotalPK, TotalAE - TotalPK) #8 8 0 6 2

dropout <- DS_dropout_all %>%
  filter(DSREAS %in% c(2:6), DS == 1) %>%
  pull(N)

table1_names <- read.csv("CSR/Table/Table1.csv", stringsAsFactors = F, header = T)
table1 <- data.frame(table1_names, N = c(Scr, ScrF, Allo, dropout))

table1
#table1 <-data.frame(table1)
#write.csv(table1, "table1.csv")


#### Table 2 ####
# code reference : https://cran.r-project.org/web/packages/qwraps2/vignettes/summary-statistics.html

DM <- read_excel("data/exceldata/DWPDN11_P01R_DataCenter_DataSet_List_20220817084019.xlsx",sheet="DM") 
LS <- read_excel("data/exceldata/DWPDN11_P01R_DataCenter_DataSet_List_20220817084019.xlsx",sheet="LS")
VS <- read_excel("data/exceldata/DWPDN11_P01R_DataCenter_DataSet_List_20220817084019.xlsx",sheet="VS")
EG <- read_excel("data/exceldata/DWPDN11_P01R_DataCenter_DataSet_List_20220817084019.xlsx",sheet="EG")
EYE <- read_excel("data/exceldata/DWPDN11_P01R_DataCenter_DataSet_List_20220817084019.xlsx",sheet="EYE")

# EYE값 전처리 / 시력검사 및 안압검사)
EYE$OD_v1 <- round(EYE$VAOD01 / EYE$VAOD02,2)
EYE$OS_v1 <- round(EYE$VAOS01 / EYE$VAOS02,2)

EYE_iop <- DS_set %>% 
  filter(DS==1) %>% 
  select(SUBJID)

EYE_2 <- EYE %>% 
  filter(VISIT == 1 & SUBJID %in% EYE_iop$SUBJID) %>% 
  group_by(SUBJID)

EYE_iop <- EYE_iop %>% 
  mutate(IOPOD_v1 = ifelse(is.na(EYE_2$IOPOD[EYE_2$IOPNUM == "중앙값"]), EYE_2$IOPOD[EYE_2$IOPNUM=="평균"], EYE_2$IOPOD[EYE_2$IOPNUM == "중앙값"])) %>% 
  mutate (IOPOS_v1 =  ifelse(is.na(EYE_2$IOPOS[EYE_2$IOPNUM == "중앙값"]), EYE_2$IOPOS[EYE_2$IOPNUM=="평균"], EYE_2$IOPOS[EYE_2$IOPNUM == "중앙값"]))


EYE <- left_join(EYE, EYE_iop, by= "SUBJID") %>% 
  filter(VISIT == 1 & IOPNUM== "평균")  #필요한 데이터가 모두 같으므로 하나만 뽑아오기 

DM_set <- left_join(DM, LS, by = c("SUBJID","VISIT"))%>%
  left_join(.,VS, by = c("SUBJID","VISIT")) %>% 
  left_join(.,EG, by = c("SUBJID","VISIT"))%>% 
  left_join(.,EYE, by = c("SUBJID","VISIT")) %>%
  left_join(.,Anal_set,  by = c("SUBJID" = "SID")) %>%
  select(SUBJID,VISIT,AGE,SEX,HT,WT,LSSMKYN, LSAHOLYN, LSCAFFYN,SYSBP,DIABP,PULSE,TEMP,EGHR, EGPR, EGQRS, EGQT, EGQTC,OD_v1, OS_v1,
         IOPOD_v1,IOPOS_v1,TBUTOD, TBUTOS, SCHOD, SCHOS, DS,SS,PS,SUBJID)%>%
  filter(DS==1) 

DM_set$PULSE <-as.numeric(DM_set$PULSE)

summary_dm <-
  list("Demographics" = 
         list("Age" = ~paste(round(mean(AGE),2),"\u00B1",round(sd(AGE),2),"(",min(AGE),"-",max(AGE),")"),
              "Height" = ~paste(round(mean(HT),2),"\u00B1",round(sd(HT),2),"(",min(HT),"-",max(HT),")"),
              "Weight" = ~paste(round(mean(WT),2),"\u00B1",round(sd(WT),2),"(",min(WT),"-",max(WT),")")
         ),
       "Vital sign" = 
         list("Systolic blood pressure" = ~paste(round(mean(SYSBP),2),"\u00B1",round(sd(SYSBP),2),"(",min(SYSBP),"-",max(SYSBP),")"),
              "Diastolic blood pressure" = ~paste(round(mean(DIABP),2),"\u00B1",round(sd(DIABP),2),"(",min(DIABP),"-",max(DIABP),")"),
              "Pulse rate" = ~paste(round(mean(PULSE),2),"\u00B1",round(sd(PULSE),2),"(",min(PULSE),"-",max(PULSE),")"),
              "Body temperature" = ~paste(round(mean(TEMP),2),"\u00B1",round(sd(TEMP),2),"(",min(TEMP),"-",max(TEMP),")")
         ),
       "ECG" = 
         list("Ventricular rate" = ~paste(round(mean(EGHR),2),"\u00B1",round(sd(EGHR),2),"(",min(EGHR),"-",max(EGHR),")"),
              "PR interval" = ~paste(round(mean(EGPR),2),"\u00B1",round(sd(EGPR),2),"(",min(EGPR),"-",max(EGPR),")"),
              "QRSD" = ~paste(round(mean(EGQRS),2),"\u00B1",round(sd(EGQRS),2),"(",min(EGQRS),"-",max(EGQRS),")"),
              "QT" = ~paste(round(mean(EGQT),2),"\u00B1",round(sd(EGQT),2),"(",min(EGQT),"-",max(EGQT),")"),
              "QTc" = ~paste(round(mean(EGQTC),2),"\u00B1",round(sd(EGQTC),2),"(",min(EGQTC),"-",max(EGQTC),")")
         ),
       "Ophthalmologic test"=
         list("Eye sight OD (decimal)" = ~paste(round(mean(OD_v1),2),"\u00B1",round(sd(OD_v1),2),"(",min(OD_v1),"-",max(OD_v1),")"),
              "Eye sight OS (decimal)" = ~paste(round(mean(OS_v1),2),"\u00B1",round(sd(OS_v1),2),"(",min(OS_v1),"-",max(OS_v1),")"),
              "IOP OD" = ~paste(round(mean(IOPOD_v1),2),"\u00B1",round(sd(IOPOD_v1),2),"(",min(IOPOD_v1),"-",max(IOPOD_v1),")"),
              "IOP OS" = ~paste(round(mean(IOPOS_v1),2),"\u00B1",round(sd(IOPOS_v1),2),"(",min(IOPOS_v1),"-",max(IOPOS_v1),")"),
              "Tear break-up time test OD" = ~paste(round(mean(TBUTOD),2),"\u00B1",round(sd(TBUTOD),2),"(",min(TBUTOD),"-",max(TBUTOD),")"),
              "Tear break-up time test OS" = ~paste(round(mean(TBUTOS),2),"\u00B1",round(sd(TBUTOS),2),"(",min(TBUTOS),"-",max(TBUTOS),")"),
              "Schirmer’s test OD" = ~paste(round(mean(SCHOD),2),"\u00B1",round(sd(SCHOD),2),"(",min(SCHOD),"-",max(SCHOD),")"),
              "Schirmer’s test OS" = ~paste(round(mean(SCHOS),2),"\u00B1",round(sd(SCHOS),2),"(",min(SCHOS),"-",max(SCHOS),")") 
               )
)

       
       
table2 <-summary_table(DM_set, summary_dm)
#table2 <-data.frame(table2)
table2
#write.csv(table1, "table2.csv")


#### Table 3 ####
MH <- read_excel("data/exceldata/DWPDN11_P01R_DataCenter_DataSet_List_20220817084019.xlsx",sheet="MH") 
MH <- MH %>% group_by(SUBJID) %>%summarise(SOC=paste(SOC,collapse=","))

MH_set <- DM_set %>% 
  mutate(mh=ifelse(SUBJID %in% MH$SUBJID,1,0)) %>%
  mutate(PE=0)%>%
  left_join(.,MH, by = c("SUBJID"="SUBJID")) %>%
  select(SUBJID, LSSMKYN, LSAHOLYN, LSCAFFYN, mh,SOC,PE) 

MH_set$SOC[is.na(MH_set$SOC)] <-0

summary_mh <- 
  list("A" = 
         list("Number of subjects no medical history" = ~n_perc(mh==0),
              "Number of subjects having clinically NOT significant medical history" = ~n_perc(mh==1)
         ),
       list("Blood and lymphatic system disorders" = ~n_perc(grepl("Blood and lymphatic system disorders",SOC)),
            "Cardiac disorders" = ~n_perc(grepl("Cardiac disorders",SOC)),
            "Congenital, familial and genetic disorders" = ~n_perc(grepl("Congenital, familial and genetic disorders",SOC)),
            "Ear and labyrinth disorders" = ~n_perc(grepl("Ear and labyrinth disorders",SOC)),
            "Endocrine disorders" = ~n_perc(grepl("Endocrine disorders",SOC)),
            "Eye disorders" = ~n_perc(grepl("Eye disorders",SOC)),
            "Gastrointestinal disorders" = ~n_perc(grepl("Gastrointestinal disorders",SOC)),
            "General disorders and administration site conditions" = ~n_perc(grepl("General disorders and administration site conditions",SOC)),
            "Hepatobiliary disorders" = ~n_perc(grepl("Hepatobiliary disorders",SOC)),
            "Infections and infestation" = ~n_perc(grepl("Infections and infestations",SOC)),
            "Immune system disorders" = ~n_perc(grepl("Immune system disorders",SOC)),
            "Injury, poisoning and procedural complications" = ~n_perc(grepl("Injury, poisoning and procedural complications",SOC)),
            "Investigations" = ~n_perc(grepl("Investigations",SOC)),
            "Metabolism and nutrition disorders" = ~n_perc(grepl("Metabolism and nutrition disorders",SOC)),
            "Musculoskeletal and connective tissue disorders" = ~n_perc(grepl("Musculoskeletal and connective tissue disorders",SOC)),
            "Neoplasms benign, malignant and unspecified (incl cysts and polyps)" = ~n_perc(grepl("Neoplasms benign, malignant and unspecified",SOC)),
            "Nervous system disorders" = ~n_perc(grepl("Nervous system disorders",SOC)),
            "Pregnancy, puerperium and perinatal conditions" = ~n_perc(grepl("Pregnancy, puerperium and perinatal conditions",SOC)),
            "Product issues" = ~n_perc(grepl("Product issues",SOC)),
            "Psychiatric disorders" = ~n_perc(grepl("Psychiatric disorders",SOC)),
            "Renal and urinary disorders" = ~n_perc(grepl("Renal and urinary disorders",SOC)),
            "Reproductive system and breast disorders" = ~n_perc(grepl("Reproductive system and breast disorders",SOC)),
            "Respiratory, thoracic and mediastinal disorders" = ~n_perc(grepl("Respiratory, thoracic and mediastinal disorders",SOC)),
            "Skin and subcutaneous tissue disorders" = ~n_perc(grepl("Skin and subcutaneous tissue disorders",SOC)),
            "Social circumstances" = ~n_perc(grepl("Social circumstances",SOC)),
            "Surgical and medical procedures" = ~n_perc(grepl("Surgical and medical procedures",SOC)),
            "Vascular disorders" = ~n_perc(grepl("Vascular disorders",SOC))
       ),
       list( "General" = ~n_perc(PE=="General"),
             "Nutrition" = ~n_perc(PE=="Nutrition"),
             "Integumentary system (skin/mucosa)" = ~n_perc(PE=="Integumentary system (skin/mucosa)"),
             "Ophthalmologic system (eye, excluding decrease visual acuity)" = ~n_perc(PE=="Ophthalmologic system (eye, excluding decrease visual acuity)"),
             "Ear, nose, & throat" = ~n_perc(PE=="Ear, nose, & throat"),
             "Thyroid" = ~n_perc(PE=="Thyroid"),
             "Respiratory system" = ~n_perc(PE=="Respiratory system"),
             "Cardiovascular system" = ~n_perc(PE=="Cardiovascular system"),
             "Abdomen" = ~n_perc(PE=="Abdomen"),
             "Kidney / Genitourinary system" = ~n_perc(PE=="Kidney / Genitourinary system"),
             "Neuropsychiatry" = ~n_perc(PE=="Neuropsychiatry"),
             "Vertebra / Limb / Any malignancies" = ~n_perc(PE=="Vertebra / Limb / Any malignancies"),
             "Peripheral blood supply" = ~n_perc(PE=="Peripheral blood supply"),
             "Lymphatics" = ~n_perc(PE=="Lymphatics"),
             "Others" = ~n_perc(PE=="Others")
       ),
       list("Number of smokers (NOT exceeding the amount indicated in exclusion criteria)" = ~n_perc(LSSMKYN ==1),
            "Number of non-smokiers" = ~n_perc(LSSMKYN == 2)
       ),
       list("Number of subjects consuming alcohol (NOT exceeding the amount indicated in exclusion criteria)" = ~n_perc( LSAHOLYN ==1),
            "Number of subjects NOT consuming alcohol" = ~n_perc( LSAHOLYN==2)
       ),
       list("Number of subjects consuming caffeine (NOT exceeding the amount indicated in exclusion criteria)" = ~n_perc(LSCAFFYN ==1),
            "Number of subjects NOT consuming caffeine" = ~n_perc(LSCAFFYN==2)
       )
  )

table3 <-summary_table(MH_set,summary_mh)
#table3    
table3 <-data.frame(table3)
table3
#write.csv(table1, "table3.csv")

#### Table 9 ####
#aeR<-c(11,24)
#aeT<-c(9,26)
#adrR<-c(12,23)
#adrT<-c(12,23)
#ae<-cbind(aeR,aeT)
#adr<-cbind(adrR,adrT)
#chisq.test(adr,correct = FALSE)