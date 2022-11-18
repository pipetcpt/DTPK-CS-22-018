#memory.limit(size=10000000024)
#rm(list=ls())


library(officer)
library(flextable) 
library(tidyverse) 
library(sas7bdat) 
library(maditr) 
library(knitr) 
library(R.utils)
library(rvg) 
library(reshape2)
library(chron) 
library(magrittr)
library(icesTAF) 
library(moonBook)
library(haven) 
library(readxl)

getwd()
# Save CSV Files

list.files('sasdata')
#####
# CSV file exist
# subDir  <- "CSV"
# sas.files <- list.files('sasdata')[grepl(".sas7bdat", list.files('sasdata'))]     
# for(files_ in sas.files){
#  name_ <- toupper(strsplit(files_, ".", fixed = TRUE)[[1]][1])
#  outfile <- write.csv(as.data.frame(read_sas(paste0('sasdata/',files_)), NULL), 
#                       file = paste0(file.path(subDir), "/", name_, ".csv"))
#}
#####

######### FUNCTIONS #########
source('appendix_code/function/count.mean.sd.r')
source('appendix_code/function/createtable_v2.r')
source('appendix_code/function/flextable.r')
source('appendix_code/function/others.r')

# ======== Starting From Here ===========

list.files(getwd())

# Load CSV Files
data.path <- paste0(getwd(), "/Data/DTC21IP070_20221028160627_SAS/SAS_20221028160627")
data.files <- list.files(data.path)
data.files

#SUBJID load
subid <- read_excel("Data/[21HD10503] CSR_Apendix_Template_2022.10.25.xlsx", sheet = '임상대상자_리스트_yj')
subid$SID <- as.character(subid$SID)
subid$RID <- as.character(subid$RID)
str(subid)
subid <- subid %>% select(SID, RID)

VST <- as.data.frame(read_sas(paste0(data.path, "/", data.files[grepl("\\<visit\\>", data.files)]), NULL))
VST$VISIT_CODE <- as.integer(VST$VISIT_CODE)
VST$VISIT_LABEL <- as.character(VST$VISIT_LABEL)


#16.2.1 중도탈락자
DS0 <- as.data.frame(read_sas(paste0(data.path, "/", data.files[grepl("\\<ds\\>", data.files)]), NULL))
DM0 <- as.data.frame(read_sas(paste0(data.path, "/", data.files[grepl("\\<dm\\>", data.files)]), NULL))
RN <- as.data.frame(read_sas(paste0(data.path, "/", data.files[grepl("\\<rn\\>", data.files)]), NULL))

DS <- DS0 %>%
  left_join(RN[,c(1,5)], by = "SUBJID") %>%
  filter(DSYN == 2) %>%
  mutate(SID = SUBJID,
         RID = RNNO,
         Final_visit = DSDTC,
         Reason_for_discontinuation = case_when(DSREAS == 1 ~ "스크리닝 탈락",
                                                DSREAS == 2 ~ "대상자가 임상시험용 의약품의 안전성이나 약동학적 특성을 평가하는데 영향을 줄 것으로 예상되는 의약품을 투약한 경우",
                                                DSREAS == 3 ~ "대상자가 임상시험용 의약품의 투약 중단을 요구하거나, 시험 참여 동의를 철회하는 경우",
                                                DSREAS == 4 ~ "중대한 이상반응/약물이상반응이 발생하여 시험자가 시험을 계속 할 수 없다고 판단하는 경우",
                                                DSREAS == 5 ~ "임상시험 중 선정/제외 기준 등 중대한 계획서 위반 사항이 새롭게 발견되는 경우",
                                                DSREAS == 6 ~ "시험자가 시험을 중지하여야 한다고 판단한 경우",
                                                DSREAS == 7 ~ "기타")) %>%
  select(SID, RID, Reason_for_discontinuation, Final_visit) %>%
  arrange(SID)

DM <- DM0 %>%
  dplyr::rename(SID = SUBJID,
                Sex = SEX,
                Age = AGE) %>%
  mutate(Sex = case_when(Sex == 1 ~ '남',
                         Sex == 2 ~ '여')) %>%
  select(SID, Sex, Age)

data <- merge(DS, DM, by="SID")
data <- data[,c("SID", 
                "RID",
                "Sex", 
                "Age", 
                "Final_visit", 
                "Reason_for_discontinuation")]

MyFTable_16.2.1 <- flex.table.fun(data)


# 16.2.4 시험대상자 특성표
SV0 <- as.data.frame(read_sas(paste0(data.path, "/", data.files[grepl("\\<sv\\>", data.files)]), NULL))
LS0 <- as.data.frame(read_sas(paste0(data.path, "/", data.files[grepl("\\<ls\\>", data.files)]), NULL))

SV <- SV0 %>%
  filter (VISIT == 1) %>% 
  mutate(SID = SUBJID) %>%
  select(SID, SVDTC) %>%
  arrange(SID)

DM <- DM0 %>% 
  mutate(SID = SUBJID,
         SEX = case_when(SEX == 1 ~ '남',
                         SEX == 2 ~ '여')) %>%
  select(SID, SEX, AGE, HT,WT) %>%
  arrange(SID)

LS <- LS0 %>%
  filter(VISIT == 1) %>% 
  mutate(SID = SUBJID,
         EXSMKYN = case_when(EXSMKYN == 1 ~ '아니오',
                             EXSMKYN == 2 ~ '10개피/일 이하',
                             EXSMKYN == 3 ~ '10개피/일 초과'),
         EXAHOLYN = case_when(EXAHOLYN == 1 ~ '아니오',
                              EXAHOLYN == 2 ~ '210g/주 이하',
                              EXAHOLYN == 3 ~ '210g/주 초과'),
         EXCAFFYN = case_when(EXCAFFYN == 1 ~ '아니오',
                              EXCAFFYN == 2 ~ '5컵/일 이하',
                              EXCAFFYN ==3 ~ '5컵/일 초과')) %>%
  select(SID, EXSMKYN, EXAHOLYN, EXCAFFYN) %>%
  arrange(SID)

DMSVLS <- merge(subid,DM,by=c("SID")) %>% 
  merge(., LS, by = c("SID"))  %>% 
  merge(., SV, by = c("SID"))


DMSVLS <- DMSVLS[c("SID", 
               "RID",
               "SVDTC",
               "AGE",
               "SEX",
               "EXSMKYN",
               "EXAHOLYN",
               "EXCAFFYN",
               "HT",
               "WT")]

colnames(DMSVLS)=c("SID", 
                 "RID",
                 "Screening Date",
                 "Age",
                 "Sex",
                 "Smoking", 
                 "Alcohol",
                 "Caffeine",
                 "Height", 
                 "Weight" )

# DMSV[is.na(DMSV)] <- " "

str(DMSVLS)
MyFTable_16.2.4 <- flex.table.fun(DMSVLS)


# 16.2.5.1 ###시험대상자별 투약시간
IP0 <- as.data.frame(read_sas(paste0(data.path, "/", data.files[grepl("\\<ip\\>", data.files)]), NULL))
IP0 <- visit.match.fun(infile = IP0, visitData = VST)

# ARM 나누기 
group_list <- IP0 %>% 
  group_by(SUBJID) %>% 
  summarise(group_num = paste(unique(IPARM), collapse = ','))

group_list <- group_list %>% 
  mutate(IP = case_when(group_num == "1,2" ~ "A",
                        group_num == "3,5" ~ "B",
                        group_num == "4,5" ~ "C", 
                        group_num == "1" ~ "A",
                        group_num == "3" ~ "B",
                        group_num == "4" ~ "C"))

IP0 <- left_join(IP0, group_list[c(1,3)] , by = ("SID" = "SUBJID"))

IP.A <- IP0 %>% 
  left_join(RN[,c(1,5)], by = "SUBJID") %>% 
  filter(IP == "A") %>% 
  mutate(scheduled = VISIT) %>% 
  select(SID=SUBJID, RID=RNNO,scheduled, done = IPTC) %>% 
  spread(key=scheduled, value =done)

IP.B <- IP0 %>% 
  left_join(RN[,c(1,5)], by = "SUBJID") %>% 
  filter(IP == "B") %>% 
  mutate(scheduled = VISIT) %>% 
  select(SID=SUBJID, RID=RNNO,scheduled, done = IPTC) %>% 
  spread(key=scheduled, value =done)

IP.C <- IP0 %>% 
  left_join(RN[,c(1,5)], by = "SUBJID") %>% 
  filter(IP == "C") %>% 
  mutate(scheduled = VISIT) %>% 
  select(SID=SUBJID, RID=RNNO,scheduled, done = IPTC) %>% 
  spread(key=scheduled, value =done)

MyFTable_16.2.5.1.a <- flex.table.fun(IP.A)
MyFTable_16.2.5.1.b <- flex.table.fun(IP.B)
MyFTable_16.2.5.1.c <- flex.table.fun(IP.C)

################################################################################################################
###################################  Park maria T part.  #######################################################
################################################################################################################
#####
# dosingtime <- as.data.frame(read_sas(paste0(data.path, "/", data.files[grepl("\\<ip\\>", data.files)]), NULL))%>% 
#   select(1,4,7) %>% 
#   left_join(RN, by="SUBJID") %>%
#   select(SID=SUBJID, RID=RNNO, "투여일"=IPDTC, "투여시간"=IPTC)
#  MyFTable_16.2.5.1 <- flex.table.fun(dosingtime)
 
# 
# #16.2.5.2 시험대상자별 혈장 내 농도
# if(!require(readxl)) install.packages("readxl");library(readxl) 
# 
# 
# N <- read_excel("data/pk/niclosamide.xlsx", sheet=5) %>% select(-1)
# HN <- read_excel("data/pk/3-hydroxyN.xlsx", sheet=5) %>% select(-1)
# MyFTable_16.2.5.2 <- flex.table.fun(N)
# 
# 
# #16.2.6.1 시험대상자별 비구획 분석 결과
# nca <- tblNCA(as.data.frame(atdt4nca), 
#               key=c("substance", "subj", "trt"),
#               colTime="actualtime",
#               colConc="conc",
#               dose=subj_dose,
#               R2ADJ=0.1) 
# 
# 
# nca_appendix <- nca %>% rename(RNNO=subj) %>% left_join(RN, by="RNNO") %>% 
#   select(1,40,2,3, AUCLST, CMAX, TMAX, LAMZHL, CLFO, VZFO) %>% 
#   arrange(substance)
# colnames(nca_appendix) <- c("Substance", "SID", "RID", "Treatment",
#                             "AUClast (hr*ng/mL)", "Cmax (ng/mL)", 
#                             "Tmax (hr)", "t1/2 (hr)", "CL/F (L/hr)", "Vd/F (L)")  

################################################################################################################
################################################################################################################
################################################################################################################

#16.2.6.5  ###시험대상자별 채혈수행시각
pb0 <- as.data.frame(read_sas(paste0(data.path, "/", data.files[grepl("\\<pb\\>", data.files)]), NULL))
pb0 <- visit.match.fun(infile = pb0, visitData = VST)

pb0 <- left_join(pb0, group_list[c(1,3)], by = "SUBJID")

pb.a <- pb0 %>% 
  left_join(RN[,c(1,5)], by="SUBJID") %>%
  filter(IP == "A") %>% 
  mutate(scheduled2=paste(parse_number(PBNT),'hr',sep=""), scheduled=paste(VISIT,scheduled2,sep="_")) %>%
  select(SID=SUBJID, RID=RNNO, scheduled, done=PBAT) %>%
  spread(key=scheduled, value=done) 

pb.b <- pb0 %>% 
  left_join(RN[,c(1,5)], by="SUBJID") %>%
  filter(IP == "B") %>% 
  mutate(scheduled2=paste(parse_number(PBNT),'hr',sep=""), scheduled=paste(VISIT,scheduled2,sep="_")) %>%
  select(SID=SUBJID, RID=RNNO, scheduled, done=PBAT) %>%
  spread(key=scheduled, value=done) 

pb.c <- pb0 %>% 
  left_join(RN[,c(1,5)], by="SUBJID") %>%
  filter(IP == "C") %>% 
  mutate(scheduled2=paste(parse_number(PBNT),'hr',sep=""), scheduled=paste(VISIT,scheduled2,sep="_")) %>%
  select(SID=SUBJID, RID=RNNO, scheduled, done=PBAT) %>%
  spread(key=scheduled, value=done) 

MyFTable_16.2.6.5.a1 <- flex.table.fun(pb.a %>% select(1,2,6,3,4,5,10,7,11:17, 8,9,18:19))
MyFTable_16.2.6.5.a2 <- flex.table.fun(pb.a %>% select(1,2,23,20:22,27,24,28:34,25,26,35:36))

MyFTable_16.2.6.5.b1 <- flex.table.fun(pb.b %>% select(1,2,6,3,4,5,10,7,11:17, 8,9,18:19))
MyFTable_16.2.6.5.b2 <- flex.table.fun(pb.b %>% select(1,2,23,20:22,27,24,28:34,25,26,35:36))

MyFTable_16.2.6.5.c1 <- flex.table.fun(pb.c %>% select(1,2,6,3,4,5,10,7,11:17, 8,9,18:19))
MyFTable_16.2.6.5.c2 <- flex.table.fun(pb.c %>% select(1,2,23,20:22,27,24,28:34,25,26,35:36))


# 16.2.7 AE-  시험대상자별 이상반응 목록
AE0 <- as.data.frame(read_sas(paste0(data.path, "/", data.files[grepl("\\<ae\\>", data.files)]), NULL))

AE <- AE0 %>%
  mutate(SID = SUBJID,
         AESEV = case_when(AESEV == 1 ~ "경증(Mild)",
                           AESEV == 2 ~ "중등증(Moderate)",
                           AESEV == 3 ~ "중증(Severe)"),
         AESER = case_when(AESER == 1 ~ "예",
                           AESER == 2 ~ "아니오"),
         AEOUT = case_when(AEOUT == 1 ~ "회복됨/해결됨 (Recovered/Resolved)",
                           AEOUT == 2 ~ "회복되었으나 후유증이 남음/해결되었으나 후유증이 남음(recovered/resolved with sequela)",
                           AEOUT == 3 ~ "회복중임/해결중임(recovering/resolving)",
                           AEOUT == 4 ~ "회복되지 않음/해결되지 않음(not recovered/not resolved)",
                           AEOUT == 5 ~ "알 수 없음(Unknown)",
                           AEOUT == 6 ~ "사망(Death)"),
         AEREL = case_when(AEREL == 1 ~ "관련성이 없음 (Not related)",
                           AEREL == 2 ~ "관련성이 적음 (Unlikely related)", 
                           AEREL == 3 ~ "관련성이 의심됨 (Possibly related)", 
                           AEREL == 4 ~ "관련성이 많음 (Probably related)", 
                           AEREL == 5 ~ "관련성이 명백 (Definitely related)", 
                           AEREL == 6 ~ "평가 불가능 (Unassessable/unclassifiable)"),
         AECON = case_when(AECON == 1 ~ "투여중지(Drug Interrupted)",
                           AECON == 2 ~ "용량변화 없음(Dose Maintained)",
                           AECON == 3 ~ "알 수 없음(Unknown)",
                           AECON == 4 ~ "해당사항 없음(Not Applicable)"),
         AEACN = case_when(AEACN == 1 ~ "없음(Not Done)",
                           AEACN == 2 ~ "약물투여(Medication)",
                           AEACN == 3 ~ "기타(Others)")) %>%
  select(SID,
         PT,
         AESTDTC,
         AESTTC,
         AEENDTC,
         AEENTC,
         AESEV,
         AESER,
         AEOUT,
         AEREL,
         AECON,
         AEACN) %>%
  arrange(SID) 

subAE <- merge(subid,AE,by=c("SID")) 

colnames(subAE) <- c("SID",
                     "RID",
                     "이상반응명",
                     "시작일",
                     "시작시간",
                     "종료일",
                     "종료시간",
                     "중증도",
                     "중대성",
                     "결과",
                     "의약품과의 인과관계",
                     "임상시험용의약품에 대한 조치",
                     "이상반응에 대한 조치")


MyFTable_16.2.7 <- flex.table.fun(subAE)


#16.2.8 AbNormal Indentity Data - 피험자별 임상검사 비정상치 자료
LB0 <- as.data.frame(read_sas(paste0(data.path, "/", data.files[grepl("\\<lb\\>", data.files)]), NULL))
subLB <- LB0 %>% 
      mutate(SID=SUBJID) %>% 
  merge(subid,.,by=c("SID")) 


# 16.4.1 antecedent drug (선행/병행약물)
CM0 <- as.data.frame(read_sas(paste0(data.path, "/", data.files[grepl("\\<cm\\>", data.files)]), NULL))

CM<- CM0 %>%
  mutate(SID = SUBJID,
         CMDOSU = case_when(CMDOSU == 1 ~ "Tablets",
                            CMDOSU == 2 ~ "Capsules",
                            CMDOSU == 3 ~ "Ampules",
                            CMDOSU == 4 ~ "Vials",
                            CMDOSU == 5 ~ "Gram",
                            CMDOSU == 6 ~ "Milligram",
                            CMDOSU == 7 ~ "Microgram",
                            CMDOSU == 8 ~ "Drops",
                            CMDOSU == 9 ~ "Liters",
                            CMDOSU == 10 ~ "Milliliters",
                            CMDOSU == 11 ~ "Puffs",
                            CMDOSU == 12 ~ "Suppository",
                            CMDOSU == 13 ~ "Ointment",
                            CMDOSU == 14 ~ "Finger Tip Unit",
                            CMDOSU == 15 ~ "Other(Specify)",
                            CMDOSU == 16 ~ "UK(Unknown)"),
         
         CMFRQ = case_when(CMFRQ == 1 ~ "OD(Once daily)",
                           CMFRQ == 2 ~ "BID(Twice daily)",
                           CMFRQ == 3 ~ "TID(Three times daily)",
                           CMFRQ == 4 ~ "QOD(Every other day)",
                           CMFRQ == 5 ~ "PRN(per requested need)",
                           CMFRQ == 6 ~ "Stat(Single immediatedose)",
                           CMFRQ == 7 ~ "Cont(Continuous)",
                           CMFRQ == 8 ~ "Other(Specify)",
                           CMFRQ == 9 ~ "UK(Unknown)"),
         
         CMROUTE = case_when(CMROUTE == 1 ~ "PO(Per Oral)",
                             CMROUTE == 2 ~ "SL(Sublingual)",
                             CMROUTE == 3 ~ "IM(Intramuscular)",
                             CMROUTE == 4 ~ "IV(Intravenous)",
                             CMROUTE == 5 ~ "SC(Subcutaneous)",
                             CMROUTE == 6 ~ "Inhalation",
                             CMROUTE == 7 ~ "Transdermal",
                             CMROUTE == 8 ~ "Topical",
                             CMROUTE == 9 ~ "Ophthalmic",
                             CMROUTE == 10 ~ "NA(Not Applicable)",
                             CMROUTE == 11 ~ "Other(Specify)",
                             CMROUTE == 12 ~ "UK(Unknown)"),
         
         CMINDC = case_when(CMINDC == 1 ~ "병력",
                            CMINDC == 2 ~ "이상반응",
                            CMINDC == 3 ~ "예방",
                            CMINDC == 4 ~ "기타(Specify)")) %>%
  select(SID,
         CMTRT,
         CMDOSTOT,
         CMDOSU,
         CMFRQ,
         CMROUTE,
         CMINDC,
         CMSTDTC,
         CMENDTC,
         ATCCD,
         ATCLV1,
         ATCLV2,
         ATCLV3,
         ATCLV4,
         ATCLV5) %>%
  arrange(SID) 

subCM <- merge(subid,CM,by=c("SID")) 

colnames(subCM) <- c("SID",
                     "RID",
                     "약물명",
                     "투여량",
                     "단위",
                     "빈도",
                     "경로",
                     "사유",
                     "시작일",
                     "종료일",
                     "ATCCD",
                     "ATCLV1",
                     "ATCLV2",
                     "ATCLV3",
                     "ATCLV4",
                     "ATCLV5")


MyFTable_16.4.1 <- flex.table.fun(subCM)


#16.4.2. ~ Chemistry lab data - 실험실적 검사
LB0 <- as.data.frame(read_sas(paste0(data.path, "/", data.files[grepl("\\<lb\\>", data.files)]), NULL))

LB <- LB0 %>% 
  mutate(SID = SUBJID,
         LBTEST = case_when(endsWith(LBTEST, "-GT") ~ "GGT",
                          TRUE ~ as.character(LBTEST))) %>%
dplyr::select(SID, VISIT, LBTEST, LBORRES) %>% dplyr::rename(TEST = LBTEST, 
                                                                RES = LBORRES) %>% 
  merge(subid,.,by=c("SID")) 


# DB Code  ;  Appendix Code 
name_map1 <- data.frame(rbind(
  # hematology test
  c("WBC", "WBC"),
  c("RBC", "RBC"),
  c("hemoglobin", "Hemoglobin"),
  c("hematocrit", "Hematocrit"),
  c("platelets", "Platelets"),
  c("Seg. neutrophil", "Seg.neutrophils"),
  c("Lymphocyte", "Lymphocyte"),
  c("Monocyte", "Monocyte"),
  c("Eosinophil", "Eosinophil"),
  c("Basophil", "Basophil"),
  c("MCV", "MCV"),
  c("MCH", "MCH"),
  c("MCHC", "MCHC"),
  
 
  # blood chemistry test
  c("Glucose", "Glucose"),
  c("BUN", "BUN"),
  c("creatinine", "creatinine"),
  c("MDRD-eGFR", "MDRD-eGFR"),
  c("total protein", "total protein"),
  c("albumin", "albumin"),
  c("AST", "AST"),
  c("ALT", "ALT"),
  c("alkaline phosphatase", "alkaline phosphatase"),
  c("total bilirubin", "total bilirubin"),
  c("GGT", "γ-GT"),
  c("uric acid", "uric acid"),
  c("Calcium", "Calcium"),
  c("phosphorus", "phosphorus"),
  c("Na", "Na"),
  c("K", "K"),
  c("Cl", "Cl"),
  c("LDH", "LDH"),
  c("CPK", "CPK"),
  c("Total cholesterol", "Total cholesterol")
  )
  )
  
colnames(name_map1) <- c("AnalyteName", "Fullname")

MyFTable_16.4.2.1To16.4.2.33 <- list()
filter_var <- "TEST"
value_var <- "RES"

for(i in 1:length(name_map1$AnalyteName)){
  MyFTable_16.4.2.1To16.4.2.33[[i]] <- flex.table.fun(
    data = create.table(csvfile    = LB,
                        code_id    = name_map1$AnalyteName[i], 
                        fullname   = name_map1$Fullname[i],
                        filter_var = filter_var,
                        value_var  = value_var,
                        period_    = NULL,
                        type_      = "numeric"),
    
    calc_ = TRUE)
  print(paste("Creating Table for", name_map1$Fullname[i], "on subject code"))
}

length(name_map1$AnalyteName) #33

#check
MyFTable_16.4.2.1To16.4.2.33[[10]]
name_map1$AnalyteName[24]


 #16.4.2.34 ~ 16.4.2.46 - 뇨검사
# Main 
name_map2 <- data.frame(rbind(
  # urine test
  c("glucose(UA)", "Glucose(UA)"),
  c("color", "Color"),
  c("leukocyte", "Leukocyte"),
  c("bilirubin", "Bilirubin"),
  c("ketone", "Ketone"),
  c("Specific gravity", "Specific gravity"),
  c("occult blood", "Occult blood"),
  c("pH", "PH"),
  c("protein", "Protein"),
  c("urobilinogen", "Urobilinogen"),
  c("nitrite", "Nitrite"),
  
  #  urine test_Microscopy
  c("WBC_microscopy", "Microscopy WBC"),
  c("RBC_microscopy", "Microscopy RBC"))

)

colnames(name_map2) <- c("AnalyteName", "Fullname")

MyFTable_16.4.2.34To16.4.2.46 <- list()
filter_var <- "TEST"
value_var <- "RES"

for(i in 1:length(name_map2$AnalyteName)){
  MyFTable_16.4.2.34To16.4.2.46[[i]] <- flex.table.fun(
    data = create.table(csvfile    = LB, 
                        code_id    = name_map2$AnalyteName[i], 
                        fullname   = name_map2$Fullname[i],
                        filter_var = filter_var,
                        value_var  = value_var,
                        period_    = NULL,
                        type_      = FALSE),
    calc_ = FALSE)
  
  print(paste("Table created for", name_map2$Fullname[i]))
} 

length(name_map2$AnalyteName) #13

#check
MyFTable_16.4.2.34To16.4.2.46[[12]]
name_map2$AnalyteName[12]

MyFTable_16.4.2.34To16.4.2.46[[13]]
name_map2$AnalyteName[13]
#####
# for(i in 1:length(name_map$AnalyteName)){
#   data = create.table(csvfile    = LB,
#                       code_id    = name_map$AnalyteName[i],
#                       fullname   = name_map$Fullname[i],
#                       filter_var = filter_var,
#                       value_var  = value_var,
#                       period_    = NULL,
#                       type_      = FALSE)
# 
#   #write.csv(data, paste0("Appendix_code/Table/16.4.2.3 요검사_", name_map$Fullname[i], ".csv"),row.names=F,fileEncoding = "cp949")
# 
# }
MyFTable_16.4.2.34To16.4.2.46[[1]]
#####

#16.4.2.44 blood coagulation test - 혈액응고검사 

BC.LAB <- LB %>% 
  filter(TEST == "PT(INR)" |
           TEST == "aPTT") %>%
  select(SID,RID, TEST, RES) %>%
  arrange(SID)
AnalyteName.pivot2 <- dcast(BC.LAB, SID ~ TEST, value.var = "RES")
subBCT <- merge(subid,AnalyteName.pivot2,by=c("SID")) 
MyFTable_16.4.2.44 <- flex.table.fun(subBCT)

#16.4.2.45 serum test - 혈청검사;

SR.LAB <- LB %>% 
  filter(TEST == "HBsAg" |
          TEST == "anti-HCV Ab" |
          TEST == "HIV Ag/Ab" |
          TEST == "Syphilis reagin test") %>%
  select(SID,RID, TEST, RES) %>%
  arrange(SID)

SR.LAB[SR.LAB==1] <- "Negative"
SR.LAB[SR.LAB==2] <- "Positive"
AnalyteName.pivot <- dcast(SR.LAB, SID ~ TEST, value.var = "RES")
subST <- merge(subid,AnalyteName.pivot,by=c("SID")) 
MyFTable_16.4.2.45 <- flex.table.fun(subST)

######
#16.4.2.48 urine drug test
# UD.LAB <- LB %>% 
#   filter(TEST == "Amphetamine" |
#            TEST == "Cocaine" |
#            TEST == "Barbiturate" |
#            TEST == "Benzodiaepine") %>%
#   select(SID,RID, TEST, RES) %>%
#   arrange(SID)
# AnalyteName.pivot3 <- dcast(UD.LAB, SID ~ TEST, value.var = "RES")
# subUDT <- merge(subid,AnalyteName.pivot3,by=c("SID")) 
# MyFTable_16.4.2.46 <- flex.table.fun(subUDT)
# 
# 
# write.csv(subUDT,"appendix/Table/16.4.2.4 소변약물 검사.csv",row.names=F,fileEncoding = "cp949")
######

#16.4.3 Medical History - 병력 
MH0 <- as.data.frame(read_sas(paste0(data.path, "/", data.files[grepl("\\<mh\\>", data.files)]), NULL))

MH <- MH0 %>%
  mutate(SID = SUBJID) %>%
  select(SID, MHSTDTC, PT, SOC)
subMH <- merge(subid,MH,by=c("SID")) 


MyFTable_16.4.3 <- flex.table.fun(subMH)


# 16.4.4. 시험대상자의 선정 
IE0 <- as.data.frame(read_sas(paste0(data.path, "/", data.files[grepl("\\<ie\\>", data.files)]), NULL))

IE <- IE0 %>%
  mutate(SID = SUBJID) %>%
  select(SID, IEYN) 
#IE.pivot <- dcast(IE, SUBJID ~ IETEST_STD, value.var = "IEORRES_STD")
subIE <- merge(subid,IE,by=c("SID"))

subIE[subIE==1] <- "Yes"
subIE[subIE==2] <- "No"


MyFTable_16.4.4 <- flex.table.fun(subIE) 


# 16.4.5 Physical examination - 신체검사 
PE0 <- as.data.frame(read_sas(paste0(data.path, "/", data.files[grepl("\\<pe\\>", data.files)]), NULL))
PE <- visit.match.fun(infile = PE0, visitData = VST)

PE <- left_join(PE, group_list[c(1,3)], by ="SUBJID")

PE <- PE %>% 
  mutate(SID = SUBJID,
         PENT = case_when(grepl("0d|7d", VISIT) & PENT == 1 ~ "투약전(0h)",
                          grepl("1d|8d", VISIT) & PENT == 1 ~ "투약 후 24hr",
                          grepl("0d|7d", VISIT) & PENT == 2 ~ "투약 후 4hr",
                          grepl("1d|8d", VISIT) & PENT == 2 ~ "투약 후 36hr",
                          PENT == 3 ~ "투약 후 8hr",
                          PENT == 4 ~ "투약 후 12hr",
                          TRUE ~ PENT))
                         
PE <- PE %>% 
  mutate(time = ifelse(!is.na(PE$PENT),paste(PE$VISIT, PE$PENT, sep = '_'), PE$VISIT))

data <- dcast(PE, SID ~ factor(time, levels = str_sort(unique(PE$time), numeric = T)), value.var = "PENOR")
data[data==1] <- "Normal Indent"
data[data==2] <- "AbNormal Indent"
data[is.na(data)] <- "NA"

subPE <- merge(subid,data,by=c("SID")) %>% 
  select("SID", "RID", starts_with("Screening"), starts_with("1기"), starts_with("2기"), everything())

MyFTable_16.4.5.1 <- flex.table.fun(subPE %>% select(1,2,3,4,8,5,6,7,9,10))
MyFTable_16.4.5.2 <- flex.table.fun(subPE %>% select(1,2,12,16,13,14,15,17,18,19,21))


# 16.4.6.1 ~ 16.4.6.5 vital signs - 활력징후
# Main
#VS0 <- as.data.frame(read_sas(paste0(data.path, "/", data.files[grepl("\\<vs\\>", data.files)]), NULL))
VS0 <- read_excel("Data/DTC21IP070_20221028160911_EXCEL/DataSetExcel_20221028160911/DTC21IP070_DataCenter_DataSet_List_20221028160911.xlsx", sheet = 'VS')
VS0 <- visit.match.fun(infile = VS0, visitData = VST)
VS0 <- VS0 %>%  filter(!is.na(VSNT)) #데이터 확인 후 결측치 값들 모두 제거해도 이상 없어서 제거함.

VS <- VS0 %>%
  mutate(SID = SUBJID, 
         newvisit= ifelse(!(is.na(VSNT)), paste(VISIT, VSNT, sep='_'),VISIT)) %>% 
  merge(subid,.,by=c("SID"))

#SYSBP
VS.SYSBP <- flex.table.fun(
  dcast(VS, SID ~ factor(newvisit,
                         levels = str_sort(
                           unique(VS$newvisit), numeric = T
                         )),
        value.var = "SYSBP") %>%
    merge(subid,.,by=c("SID")) %>% 
    select("SID", "RID", starts_with("Screening"), starts_with("1기"), everything()),  calc_ = TRUE)


# BPSYS
VS.BPSYS <- flex.table.fun(dcast(VS, SID ~ factor(newvisit,levels = str_sort(unique(VS$newvisit), numeric = T)),value.var = "DIABP") %>%
                             merge(subid,.,by=c("SID")) %>% 
                             select("SID", "RID", starts_with("Screening"), starts_with("1기"), everything()), calc_ = TRUE)


# PULSE
VS.PL <- flex.table.fun(
  dcast(VS, SID ~ factor(newvisit,
                               levels = str_sort(
                                 unique(VS$newvisit), numeric = T)),
        value.var = "PULSE") %>%
    merge(subid,.,by=c("SID")) %>% 
    select("SID", "RID", starts_with("Screening"), starts_with("1기"), everything()),calc_ = TRUE)


# TEMP
VS.TEMP <- flex.table.fun(
  dcast(VS, SID ~ factor(newvisit,
                               levels = str_sort(
                                 unique(VS$newvisit), numeric = T)),
        value.var = "TEMP") %>%
    merge(subid,.,by=c("SID")) %>% 
    select("SID", "RID", starts_with("Screening"), starts_with("1기"), everything()),calc_ = TRUE)

# 16.4.7 - LEAD 심전도 
# Main
EG0 <- as.data.frame(read_sas(paste0(data.path, "/", data.files[grepl("\\<eg\\>", data.files)]), NULL))
EG <- visit.match.fun(infile = EG0, visitData = VST) %>% filter(is.na(EGND))
EG <- EG %>%  rename("SID" = "SUBJID")


EGVR <- flex.table.fun(dcast(EG, SID~factor(VISIT, levels=str_sort(unique(EG$VISIT), numeric=T)), value.var = "EGVR")%>%
                          merge(subid,.,by=c("SID")) %>% 
                         select("SID", "RID", "Screening",starts_with("1기"), starts_with("2기"), everything()), calc_= TRUE)

EGPR <- flex.table.fun(dcast(EG, SID~factor(VISIT, levels=str_sort(unique(EG$VISIT), numeric=T)), value.var = "EGPR")%>%
                          merge(subid,.,by=c("SID")) %>% 
                         select("SID", "RID", "Screening",starts_with("1기"), starts_with("2기"),everything()),calc_= TRUE)

EGQRSD <- flex.table.fun(dcast(EG, SID~factor(VISIT, levels=str_sort(unique(EG$VISIT), numeric=T)), value.var = "EGQRSD")%>%
                          merge(subid,.,by=c("SID")) %>% 
                           select("SID", "RID", "Screening",starts_with("1기"), starts_with("2기"), everything()),calc_= TRUE)

EGQT <- flex.table.fun(dcast(EG, SID~factor(VISIT, levels=str_sort(unique(EG$VISIT), numeric=T)), value.var = "EGQT")%>%
                          merge(subid,.,by=c("SID")) %>% 
                         select("SID", "RID", "Screening",starts_with("1기"), starts_with("2기"), everything()),calc_= TRUE)

EGQTC <- flex.table.fun(dcast(EG, SID~factor(VISIT, levels=str_sort(unique(EG$VISIT), numeric=T)), value.var = "EGQTC")%>%
                          merge(subid,.,by=c("SID")) %>% 
                          select("SID", "RID", "Screening",starts_with("1기"), starts_with("2기"), everything()),calc_= TRUE)

ECGNORM <- EG %>%
  select(SID, VISIT, EGNOR) %>%
  distinct(SID, VISIT, EGNOR) %>%
  mutate(EGNOR = case_when(EGNOR == 1 ~ "정상",
                           EGNOR == 2 ~ "비정상/NCS",
                           EGNOR == 3 ~ "비정상/CS",
                           TRUE ~ as.character(EGNOR)))

ECGNORM1 <- dcast(ECGNORM, SID~VISIT, value.var = "EGNOR") %>% merge(subid,.,by=c("SID")) %>% 
  select("SID", "RID", "Screening", starts_with("1기"), starts_with("2기"), everything())

ECGNORM1[is.na(ECGNORM1)] <- "NA"

ECGNORM <- flex.table.fun(ECGNORM1, calc_ = FALSE)



# 16.4.8 혈당
# 16.4.8.1 BST
BS0 <- read_excel("Data/DTC21IP070_20221028160911_EXCEL/DataSetExcel_20221028160911/DTC21IP070_DataCenter_DataSet_List_20221028160911.xlsx", sheet = 'BS')
#BS0 <- as.data.frame(read_sas(paste0(data.path, "/", data.files[grepl("\\<bs\\>", data.files)]), NULL))

BS0 <- visit.match.fun(infile = BS0, visitData = VST)
BS0 <- BS0 %>%  filter(!is.na(BSNT))

BS <- BS0 %>% 
  mutate(SID = SUBJID,
         newvisit = ifelse(!(is.na(SEQ)), paste(VISIT, BSNT, sep='_'),VISIT)) %>% 
  merge(subid,.,by=c("SID"))

BST <- flex.table.fun(
  dcast(BS, SID ~ factor(newvisit, levels = str_sort(unique(BS$newvisit), numeric = T)), value.var = "BSVAL") %>%
    merge(subid,., by=c("SID")) %>% 
    select(SID,RID, starts_with("1기"), starts_with("2기"), everything()),  calc_ = TRUE)


#arm 별로
#####
# # arm A 
# BS.a <- BS %>%
#   filter(IP == "A") %>% 
#   mutate(SID = SUBJID, 
#          newvisit= ifelse(!(is.na(SEQ)), paste(VISIT, BSNT, sep='_'),VISIT)) %>% 
#   merge(subid,.,by=c("SID"))
# 
# 
# BST.a <- flex.table.fun(
#   dcast(BS.a, SID ~ factor(newvisit,
#                          levels = str_sort(
#                            unique(BS.a$newvisit), numeric = T)),
#         value.var = "BSVAL") %>%
#     merge(subid,., by=c("SID")) %>% 
#     select(SID,RID, starts_with("1기"), everything()),  calc_ = TRUE)
# 
# 
# # arm B 
# BS.b <- BS %>%
#   filter(IP == "B") %>% 
#   mutate(SID = SUBJID, 
#          newvisit= ifelse(!(is.na(SEQ)), paste(VISIT, BSNT, sep='_'),VISIT)) %>% 
#   merge(subid,.,by=c("SID"))
# 
# 
# BST.b <- flex.table.fun(
#   dcast(BS.b, SID ~ factor(newvisit,levels = str_sort(unique(BS.b$newvisit), numeric = T)),
#         value.var = "BSVAL") %>%
#     merge(subid,., by=c("SID")) %>% 
#     select("SID","RID", starts_with("1기"), everything()),  calc_ = TRUE)
# 
# # arm C 
# BS.c <- BS %>%
#   filter(IP == "C") %>% 
#   mutate(SID = SUBJID, 
#          newvisit= ifelse(!(is.na(SEQ)), paste(VISIT, BSNT, sep='_'),VISIT)) %>% 
#   merge(subid,.,by=c("SID"))
# 
# 
# BST.c <- flex.table.fun(
#   dcast(BS.c, SID ~ factor(newvisit,levels = str_sort(unique(BS.c$newvisit), numeric = T)),
#         value.var = "BSVAL") %>%
#     merge(subid,., by=c("SID")) %>% 
#     select("SID","RID", starts_with("1기"), everything()),  calc_ = TRUE)



# # 16.4.8.2 BST follow up
# BF0 <- as.data.frame(read_sas(paste0(data.path, "/", data.files[grepl("\\<bf\\>", data.files)]), NULL))
# BF0 <- visit.match.fun(infile = BF0, visitData = VST)
# BF <- BF0[,c(1,2,4,5)] %>% arrange(VISIT)
# colnames(BF) <- c("SID","visit","측정시점","결과(mg/dL)")
# subBF <- flex.table.fun(merge(subid,BF,by=c("SID")))



#안과검사
# 안과 검사 
#####
# 16.4.8.1. ~ 16.4.8.2. 시력검사 (검사 1번)
# EYE0 <- as.data.frame(read_sas(paste0(data.path, "/", data.files[grepl("\\<eye\\>", data.files)]), NULL))
# 
# 
# EYE0$OD_v1 <- round(EYE0$VAOD01 / EYE0$VAOD02,2) #시력검사 결과 산출
# EYE0$OS_v1 <- round(EYE0$VAOS01 / EYE0$VAOS02,2) #시력검사 결과 산출 
# EYE0 <-visit.match.fun(infile= EYE0, visitData = VST)
# EYE0$time <-  paste(EYE0$VISIT,EYE0$IOPNUM, sep="_")
# 
# EYE1 <- EYE0 %>% 
#   group_by(SUBJID, VISIT) %>% 
#   select(SID=SUBJID, VISIT, OD_v1, OS_v1) %>% 
#   summarise(OD_v1 = unique(OD_v1), OS_v1 = unique(OS_v1)) %>% 
#   ungroup() 
# 
# #OD
# MyFTable_16.4.8.1 <- flex.table.fun(
#   spread(data = EYE1[,c(-4)],key = VISIT, value = OD_v1) %>% 
#     merge(subid,., by= c("SID")) %>% 
#   select("SID","RID","Screening  Visit","투여일(5d)","UV1"))
# 
# #OS    
# MyFTable_16.4.8.2 <- flex.table.fun(
#   spread(data = EYE1[,c(-3)],key = VISIT, value = OS_v1) %>% 
#     merge(subid,., by= c("SID")) %>% 
#   select("SID","RID","Screening  Visit","투여일(5d)","UV1"))
# 
# 
# # 16.4.8.3.~16.4.8.4   안압검사 (검사 여러번)
# #OD
# MyFTable_16.4.8.3 <- flex.table.fun(
#   spread(data = EYE0[,c(1,9,19)],key = time, value = IOPOD) %>% 
#     merge(subid,., by.x="SID", by.y="SUBJID") %>%
#     select("SID","RID", "Screening  Visit_1회", "Screening  Visit_2회","Screening  Visit_평균",
#            "Screening  Visit_3회","Screening  Visit_중앙값",
#            "투여일(5d)_1회", "투여일(5d)_2회", "투여일(5d)_평균",
#            "투여일(5d)_3회","투여일(5d)_중앙값",
#            "UV1_1회", "UV1_2회", "UV1_평균",
#            "UV1_3회", "UV1_중앙값", ))
#   
# #OS    
# MyFTable_16.4.8.4 <- flex.table.fun(
#   spread(data = EYE0[,c(1,10,19)],key = time, value = IOPOS) %>% 
#     merge(subid,., by.x="SID", by.y="SUBJID") %>% 
#     select("SID","RID", "Screening  Visit_1회", "Screening  Visit_2회","Screening  Visit_평균",
#            "Screening  Visit_3회","Screening  Visit_중앙값",
#            "투여일(5d)_1회", "투여일(5d)_2회", "투여일(5d)_평균",
#            "투여일(5d)_3회","투여일(5d)_중앙값", 
#            "UV1_1회", "UV1_2회", "UV1_평균",
#            "UV1_3회","UV1_중앙값"))
# 
# 
# # 16.4.8.5 ~ 16.8.6 눈물막 파괴 사건 (검사 1번)
# EYE2 <- EYE0 %>%
#   group_by(SUBJID, VISIT) %>% 
#   select(SID=SUBJID, VISIT, TBUTOD, TBUTOS, SCHOD, SCHOS, SLOD, SLOS) %>% 
#   summarise(TBUTOD = unique(TBUTOD), TBUTOS = unique(TBUTOS),SCHOD = unique(SCHOD),SCHOS = unique(SCHOS),SLOD = unique(SLOD),SLOS = unique(SLOS)) %>% 
#   ungroup() 
# 
# EYE2 <- EYE2 %>% 
#   mutate(SLOD = case_when(SLOD == 1 ~ "grade 0",
#                           SLOD == 2 ~ "grade Ⅰ",
#                           SLOD == 3 ~ "grade Ⅱ",
#                           SLOD == 4 ~ "grade Ⅲ",
#                           SLOD == 5 ~ "grade Ⅳ",
#                           SLOD == 6 ~ "grade Ⅴ"),
#          SLOS = case_when(SLOS == 1 ~ "grade 0",
#                           SLOS == 2 ~ "grade Ⅰ",
#                           SLOS == 3 ~ "grade Ⅱ",
#                           SLOS == 4 ~ "grade Ⅲ",
#                           SLOS == 5 ~ "grade Ⅳ",
#                           SLOS == 6 ~ "grade Ⅴ"))
# 
# #OD
# MyFTable_16.4.8.5 <- flex.table.fun(
#   spread(data = EYE2[,c(1,2,3)],key = VISIT, value = TBUTOD) %>% 
#     merge(subid,., by = "SID") %>%
#     select("SID","RID", "Screening  Visit","투여일(5d)","UV1"))
# 
# 
# #OS    
# MyFTable_16.4.8.6 <- flex.table.fun(
#   spread(data = EYE2[,c(1,2,4)],key = VISIT, value = TBUTOS) %>% 
#     merge(subid,., by = "SID") %>%
#     select("SID","RID", "Screening  Visit","투여일(5d)","UV1"))
# 
# 
# 
# # 16.4.8.7 ~ 16.8.8 쉬르머 검사 (검사 1번)
# #OD
# MyFTable_16.4.8.7 <- flex.table.fun(
#   spread(data = EYE2[,c(1,2,5)],key = VISIT, value = SCHOD) %>% 
#     merge(subid,., by = "SID") %>%
#     select("SID","RID", "Screening  Visit","투여일(5d)","UV1"))
# 
# 
# #OS    
# MyFTable_16.4.8.8 <- flex.table.fun(
#   spread(data = EYE2[,c(1,2,6)],key = VISIT, value = SCHOS) %>% 
#     merge(subid,., by = "SID") %>%
#     select("SID","RID", "Screening  Visit","투여일(5d)","UV1"))
# 
# 
# # 16.4.8.9 ~ 16.8.10 세극등 검사 (검사 1번)
# #OD
# MyFTable_16.4.8.9 <- flex.table.fun(
#   spread(data = EYE2[,c(1,2,7)],key = VISIT, value = SLOD) %>% 
#     merge(subid,., by = "SID") %>%
#     select("SID","RID", "Screening  Visit","투여일(5d)","UV1"))
# 
# 
# #OS    
# MyFTable_16.4.8.10 <- flex.table.fun(
#   spread(data = EYE2[,c(1,2,8)],key = VISIT, value = SLOS) %>% 
#     merge(subid,., by = "SID") %>%
#     select("SID","RID", "Screening  Visit","투여일(5d)","UV1"))



# # 16.4.8 혈당
# # 16.4.8.1 BST
# BS0 <- as.data.frame(read_sas(paste0(data.path, "/", data.files[grepl("\\<bs\\>", data.files)]), NULL))
# BS0 <- visit.match.fun(infile = BS0, visitData = VST)
# BS <- BS0 %>%
#   mutate(SID = SUBJID, 
#          newvisit= ifelse(!(is.na(SEQ)), paste(VISIT, BSNT, sep='_'),VISIT)) %>% 
#   merge(subid,.,by=c("SID"))
# 
# 
# BST <- flex.table.fun(
#   dcast(BS, SID ~ factor(newvisit,
#                          levels = str_sort(
#                            unique(BS$newvisit), numeric = T
#                          )),
#         value.var = "BSREASS") %>%
#     merge(subid,.,by=c("SID")) %>% 
#     dplyr::select(SID,RID, starts_with("1기"), everything()),  calc_ = TRUE
#   
# )
# 
# 
# # 16.4.8.2 BST follow up
# BF0 <- as.data.frame(read_sas(paste0(data.path, "/", data.files[grepl("\\<bf\\>", data.files)]), NULL))
# BF0 <- visit.match.fun(infile = BF0, visitData = VST)
# BF <- BF0[,c(1,2,4,5)] %>% arrange(VISIT)
# colnames(BF) <- c("SID","visit","측정시점","결과(mg/dL)")
# subBF <- flex.table.fun(merge(subid,BF,by=c("SID")))






##############################
###       Doc & Table      ###
##############################

# old <- set_flextable_defaults(font.family = "Times New Roman")
# do.call(set_flextable_defaults,old)

doc <- read_docx(path = "C:/Users/Owner/Documents/GitHub/Template/Template.docx")

doc <- body_add_par(doc, "Appendix", style = "heading 1") #1
doc <- body_add_par(doc, "시험정보", style = "heading 2") #1.1
doc <- body_add_par(doc, "임상시험 계획서", style = "heading 3") #1.1.1
doc <- body_add_par(doc, "Refer to ISF", style = "Normal Indent")
doc <- body_add_break(doc)
doc <- body_add_par(doc, "증례기록서 양식", style = "heading 3") #1.1.2
doc <- body_add_par(doc, "Refer to ISF", style = "Normal Indent")
doc <- body_add_break(doc)
doc <- body_add_par(doc, "임상시험 심사위원회의 명단 및 심사기록, 서면동의서와 동의를 위한 설명서 양식", style = "heading 3") #1.1.3
doc <- body_add_par(doc, "Refer to ISF", style = "Normal Indent")
doc <- body_add_break(doc)
doc <- body_add_par(doc, "시험책임자 및 담당자, 공동연구자의 명단 및 시험수행에 적합한 약력", style = "heading 3") #1.1.4
doc <- body_add_par(doc, "Refer to ISF", style = "Normal Indent")
doc <- body_add_break(doc)
doc <- body_add_par(doc, "시험책임자 및 공동연구자 또는 시험 의뢰자의 서명", style = "heading 3") #1.1.5
doc <- body_add_par(doc, "Refer to section 1.2", style = "Normal Indent")
doc <- body_add_break(doc)
doc <- body_add_par(doc, "2가지 이상의 batch에서 제조된 임상시험용의약품을 투여 받았을 경우 그 제조번호 목록 및 각 batch별 투여된 피험자의 명단", style = "heading 3") #1.1.6
doc <- body_add_par(doc, "Not Applicable", style = "Normal Indent")
doc <- body_add_break(doc)
doc <- body_add_par(doc, "무작위 배정 방법 및 배정표", style = "heading 3") #1.1.7
doc <- body_add_par(doc, "Refer to ISF", style = "Normal Indent")
doc <- body_add_break(doc)
doc <- body_add_par(doc, "점검확인서",  style = "heading 3") #1.1.8
doc <- body_add_par(doc, "Refer to Audit Certificates", style = "Normal Indent")
doc <- body_add_break(doc)
doc <- body_add_par(doc, "통계적 방법에 관한 문서", style = "heading 3") #1.1.9
doc <- body_add_par(doc, "Refer to statistical analysis plan", style = "Normal Indent")
doc <- body_add_break(doc)
doc <- body_add_par(doc, "각 실시기관의 실험실간 표준화 방법, 기타 자료의 질적 보증을 위해 사용한 방법", style = "heading 3") #1.1.10
doc <- body_add_par(doc, "Refer to ISF", style = "Normal Indent")
doc <- body_add_break(doc)
doc <- body_add_par(doc, "임상시험결과를 출판하였을 경우 출판물", style = "heading 3") #1.1.11
doc <- body_add_par(doc, "Not Applicable", style = "Normal Indent")
doc <- body_add_break(doc)
doc <- body_add_par(doc, "임상시험결과의 평가에 절대적 영향을 미친 참고문헌", style = "heading 3") #1.1.12
doc <- body_add_par(doc, "Not Applicable", style = "Normal Indent")
doc <- body_add_break(doc)
doc <- body_add_par(doc, "시험대상자 자료 목록", style = "heading 2") #1.2

doc <- body_add_par(doc, "중도탈락자", style = "heading 3")	#1.2.1
doc <- body_add_flextable(doc, MyFTable_16.2.1)
doc <- body_add_break(doc)


doc <- body_add_par(doc, "시험계획서 이탈", style = "heading 3") #1.2.2
doc <- body_add_par(doc, "Refer to section 10.2", style = "Normal Indent")
doc <- body_add_break(doc)

doc <- body_add_par(doc, "약동학 평가에서 제외된 시험대상자", style = "heading 3") #1.2.3
doc <- body_add_par(doc, "Not Applicable", style = "Normal Indent")
doc <- body_add_break(doc)

doc <- body_add_par(doc, "시험대상자 특성표", style = "heading 3") #1.2.4
doc <- body_add_flextable(doc, MyFTable_16.2.4)
doc <- body_add_break(doc)

doc <- body_end_section_continuous(doc) # 가로 시작
doc <- body_add_par(doc, "순응도 및 혈중농도 자료", style = "heading 3") #1.2.5
doc <- body_add_par(doc, "16.2.5.1 시험대상자별 개별 투약 시간", style="Normal Indent")

doc <- body_add_par(doc, "■ Arm A", style = "Normal Indent" )
doc <- body_add_flextable(doc, MyFTable_16.2.5.1.a)
doc <- body_add_break(doc)

doc <- body_add_par(doc, "■ Arm B", style = "Normal Indent" )
doc <- body_add_flextable(doc, MyFTable_16.2.5.1.b)
doc <- body_add_break(doc)

doc <- body_add_par(doc, "■ Arm C", style = "Normal Indent" )
doc <- body_add_flextable(doc, MyFTable_16.2.5.1.c)
doc <- body_add_break(doc)

doc <- body_end_section_landscape(doc) #가로 끝
######################################
doc <- body_add_par(doc, "1.2.5.2 시험대상자별 혈장 내 농도 (Period 1 Tamsulosin)", style="Normal Indent")
doc <- body_add_break(doc)
doc <- body_add_par(doc, "1.2.5.3 시험대상자별 혈장 내 농도 (Period 2 Mirabegron)", style="Normal Indent")
doc <- body_add_break(doc)
doc <- body_add_par(doc, "1.2.5.4 시험대상자별 혈장 내 농도 (Period 3 Tamsulosin)", style="Normal Indent")
doc <- body_add_break(doc)
doc <- body_add_par(doc, "1.2.5.5 시험대상자별 혈장 내 농도 (Period 3 Mirabegron)", style="Normal Indent")
doc <- body_add_break(doc)
######################################
doc <- body_add_par(doc, "시험대상자별 약동학 결과", style = "heading 3") #16.2.6
doc <- body_add_par(doc, "1.2.6.1 시험대상자별 비구획 분석 결과 (Tamsulosin)", style="Normal Indent")
doc <- body_add_break(doc)
doc <- body_add_par(doc, "1.2.6.2 시험대상자별 비구획 분석 결과 (Mirabegron)", style="Normal Indent")
doc <- body_add_break(doc)
doc <- body_add_par(doc, "1.2.6.3 시험대상자별 혈중농도 프로파일 (Tamsulosin)", style="Normal Indent")
doc <- body_add_break(doc)
doc <- body_add_par(doc, "1.2.6.4 시험대상자별 혈중농도 프로파일 (Mirabegron)", style="Normal Indent")
doc <- body_add_break(doc)
####################################
doc <- body_end_section_continuous(doc) # 가로 시작
doc <- body_add_par(doc, "1.2.6.5 시험대상자별 채혈수행시각", style="Normal Indent") 

doc <- body_add_par(doc, "■ Arm A", style = "Normal Indent" )
doc <- body_add_flextable(doc, MyFTable_16.2.6.5.a1)
doc <- body_add_break(doc)
doc <- body_add_flextable(doc, MyFTable_16.2.6.5.a2)
doc <- body_add_break(doc)

doc <- body_add_par(doc, "■ Arm B", style = "Normal Indent" )
doc <- body_add_flextable(doc, MyFTable_16.2.6.5.b1)
doc <- body_add_break(doc)
doc <- body_add_flextable(doc, MyFTable_16.2.6.5.b2)
doc <- body_add_break(doc)

doc <- body_add_par(doc, "■ Arm C", style = "Normal Indent" )
doc <- body_add_flextable(doc, MyFTable_16.2.6.5.c1)
doc <- body_add_break(doc)
doc <- body_add_flextable(doc, MyFTable_16.2.6.5.c2)
doc <- body_add_break(doc)

#doc <- body_end_section_continuous(doc) # 가로 시작
doc <- body_end_section_landscape(doc) #가로 끝
doc <- body_add_par(doc, "시험대상자별 이상반응 목록", style = "heading 3") #1.2.7
doc <- body_add_flextable(doc, MyFTable_16.2.7)
doc <- body_add_break(doc)
#doc <- body_end_section_landscape(doc) #가로 끝
doc <- body_add_par(doc, "피험자별 임상검사 비정상치 자료", style = "heading 3") #1.2.8
doc <- body_add_break(doc)


doc <- body_add_par(doc, "증례기록서", style = "heading 2") #1.3
doc <- body_add_par(doc, "사망, 또는 다른 중대한 이상반응, 다른 유의성 있는 이상반응에 대한 증례기록서", style = "heading 3") #1.3.1
doc <- body_add_par(doc, "Not Applicable", style = "Normal Indent")
doc <- body_add_par(doc, "이외의 증례기록서", style = "heading 3") #1.3.2
doc <- body_add_par(doc, "Not Applicable", style = "Normal Indent")
doc <- body_add_break(doc)

doc <- body_end_section_continuous(doc) # 가로 시작
doc <- body_add_par(doc, "시험대상자별 자료 목록", style = "heading 2") # 1.4
doc <- body_add_par(doc, "선행약물 또는 병행 약물", style = "heading 3") # 1.4.1
doc <- body_add_flextable(doc, MyFTable_16.4.1)
doc <- body_add_break(doc)
#doc <- body_end_section_landscape(doc) #가로 끝


doc <- body_add_par(doc, "실험실적 검사", style = "heading 3") 	#1.4.2
doc <- body_add_par(doc, "16.4.2.1 혈액학 검사", style="heading 3")

doc <- body_add_par(doc,paste0("■",name_map1$AnalyteName[1]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[1]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■",name_map1$AnalyteName[2]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[2]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■",name_map1$AnalyteName[3]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[3]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■",name_map1$AnalyteName[4]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[4]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■",name_map1$AnalyteName[5]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[5]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■",name_map1$AnalyteName[6]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[6]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■",name_map1$AnalyteName[7]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[7]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■",name_map1$AnalyteName[8]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[8]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■",name_map1$AnalyteName[9]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[9]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■",name_map1$AnalyteName[10]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[10]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■",name_map1$AnalyteName[11]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[11]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■",name_map1$AnalyteName[12]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[12]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■",name_map1$AnalyteName[13]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[13]])
doc <- body_add_break(doc)
#####
# doc <- body_add_par(doc, "16.4.2.2 혈액학 검사", style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[2]])
# doc <- body_add_break(doc)
# doc <- body_add_par(doc, "16.4.2.3 혈액학 검사", style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[3]])
# doc <- body_add_break(doc)
# doc <- body_add_par(doc, "16.4.2.4 혈액학 검사", style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[4]])
# doc <- body_add_break(doc)
# 
# doc <- body_add_par(doc, "16.4.2.5 혈액학 검사", style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[5]])
# doc <- body_add_break(doc)
# doc <- body_add_par(doc, "16.4.2.6 혈액학 검사", style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[6]])
# doc <- body_add_break(doc)
# doc <- body_add_par(doc, "16.4.2.7 혈액학 검사", style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[7]])
# doc <- body_add_break(doc)
# doc <- body_add_par(doc, "16.4.2.8 혈액학 검사", style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[8]])
# doc <- body_add_break(doc)
# 
# doc <- body_add_par(doc, "16.4.2.9 혈액학 검사", style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[9]])
# doc <- body_add_break(doc)
# doc <- body_add_par(doc, "16.4.2.10 혈액학 검사", style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[10]])
# doc <- body_add_break(doc)
# doc <- body_add_par(doc, "16.4.2.11 혈액학 검사", style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[11]])
# doc <- body_add_break(doc)
# doc <- body_add_par(doc, "16.4.2.12 혈액학 검사", style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[12]])
# doc <- body_add_break(doc)
# doc <- body_add_par(doc, "16.4.2.13 혈액학 검사", style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[13]])
# doc <- body_add_break(doc)
#####

doc <- body_add_par(doc, "1.4.2.2 혈액화학 검사", style="heading 3")

doc <- body_add_par(doc,paste0("■",name_map1$AnalyteName[14]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[14]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■",name_map1$AnalyteName[15]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[15]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■",name_map1$AnalyteName[16]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[16]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■",name_map1$AnalyteName[17]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[17]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■",name_map1$AnalyteName[18]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[18]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■",name_map1$AnalyteName[19]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[19]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■",name_map1$AnalyteName[20]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[20]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■",name_map1$AnalyteName[21]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[21]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■",name_map1$AnalyteName[22]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[22]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■",name_map1$AnalyteName[23]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[23]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■",name_map1$AnalyteName[24]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[24]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■",name_map1$AnalyteName[25]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[25]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■",name_map1$AnalyteName[26]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[26]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■",name_map1$AnalyteName[27]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[27]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■",name_map1$AnalyteName[28]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[28]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■",name_map1$AnalyteName[29]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[29]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■",name_map1$AnalyteName[30]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[30]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■",name_map1$AnalyteName[31]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[31]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■",name_map1$AnalyteName[32]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[32]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■",name_map1$AnalyteName[33]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[33]])
doc <- body_add_break(doc)

#####
# doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[14]])
# doc <- body_add_break(doc)
# doc <- body_add_par(doc, "16.4.2.15 혈액화학 검사", style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[15]])
# doc <- body_add_break(doc)
# doc <- body_add_par(doc, "16.4.2.16 혈액화학 검사", style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[16]])
# doc <- body_add_break(doc)
# doc <- body_add_par(doc, "16.4.2.17 혈액화학 검사", style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[17]])
# doc <- body_add_break(doc)
# doc <- body_add_par(doc, "16.4.2.18 혈액화학 검사", style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[18]])
# doc <- body_add_break(doc)
# 
# doc <- body_add_par(doc, "16.4.2.19 혈액화학 검사", style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[19]])
# doc <- body_add_break(doc)
# doc <- body_add_par(doc, "16.4.2.20 혈액화학 검사", style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[20]])
# doc <- body_add_break(doc)
# doc <- body_add_par(doc, "16.4.2.21 혈액화학 검사", style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[21]])
# doc <- body_add_break(doc)
# doc <- body_add_par(doc, "16.4.2.22 혈액화학 검사", style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[22]])
# doc <- body_add_break(doc)

#doc <- body_end_section_continuous(doc) # 가로 시작
# doc <- body_add_par(doc, "16.4.2.23 혈액화학 검사", style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[23]])
# doc <- body_add_break(doc)
#doc <- body_end_section_landscape(doc) #가로 끝

# doc <- body_add_par(doc, "16.4.2.24 혈액화학 검사", style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[24]])
# doc <- body_add_break(doc)
# doc <- body_add_par(doc, "16.4.2.25 혈액화학 검사", style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[25]])
# doc <- body_add_break(doc)
# doc <- body_add_par(doc, "16.4.2.26 혈액화학 검사", style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[26]])
# doc <- body_add_break(doc)
# doc <- body_add_par(doc, "16.4.2.27 혈액화학 검사", style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[27]])
# doc <- body_add_break(doc)
# doc <- body_add_par(doc, "16.4.2.28 혈액화학 검사", style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[28]])
# doc <- body_add_break(doc)
# doc <- body_add_par(doc, "16.4.2.29 혈액화학 검사", style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[29]])
# doc <- body_add_break(doc)
# doc <- body_add_par(doc, "16.4.2.30 혈액화학 검사", style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[30]])
# doc <- body_add_break(doc)
# 
# doc <- body_add_par(doc, "16.4.2.31 혈액화학 검사", style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[31]])
# doc <- body_add_break(doc)
# doc <- body_add_par(doc, "16.4.2.32 혈액화학 검사", style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[32]])
# doc <- body_add_break(doc)
#doc <- body_add_par(doc, "16.4.2.33 혈액화학 검사", style="Normal Indent")
#doc <- body_add_flextable(doc, MyFTable_16.4.2.1To16.4.2.33[[33]])
#doc <- body_add_break(doc)
#####

doc <- body_add_par(doc, "16.4.2.3 뇨검사",style="heading 3")

doc <- body_add_par(doc,paste0("■", name_map2$AnalyteName[1]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.34To16.4.2.46[[1]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■", name_map2$AnalyteName[2]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.34To16.4.2.46[[2]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■", name_map2$AnalyteName[3]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.34To16.4.2.46[[3]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■", name_map2$AnalyteName[4]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.34To16.4.2.46[[4]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■", name_map2$AnalyteName[5]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.34To16.4.2.46[[5]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■", name_map2$AnalyteName[6]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.34To16.4.2.46[[6]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■", name_map2$AnalyteName[7]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.34To16.4.2.46[[7]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■", name_map2$AnalyteName[8]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.34To16.4.2.46[[8]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■", name_map2$AnalyteName[9]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.34To16.4.2.46[[9]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■", name_map2$AnalyteName[10]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.34To16.4.2.46[[10]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■", name_map2$AnalyteName[11]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.34To16.4.2.46[[11]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■", name_map2$AnalyteName[12]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.34To16.4.2.46[[12]])
doc <- body_add_break(doc)

doc <- body_add_par(doc,paste0("■", name_map2$AnalyteName[13]), style = "Normal Indent")
doc <- body_add_flextable(doc, MyFTable_16.4.2.34To16.4.2.46[[13]])
doc <- body_add_break(doc)
#####
# doc <- body_add_par(doc, "16.4.2.35 뇨검사",style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.34To16.4.2.46[[2]])
# doc <- body_add_break(doc)
# doc <- body_add_par(doc, "16.4.2.36 뇨검사",style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.34To16.4.2.46[[3]])
# doc <- body_add_break(doc)
# doc <- body_add_par(doc, "16.4.2.37 뇨검사",style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.34To16.4.2.46[[4]])
# doc <- body_add_break(doc)
# doc <- body_add_par(doc, "16.4.2.38 뇨검사",style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.34To16.4.2.46[[5]])
# doc <- body_add_break(doc)
# doc <- body_add_par(doc, "16.4.2.39 뇨검사",style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.34To16.4.2.46[[6]])
# doc <- body_add_break(doc)
# doc <- body_add_par(doc, "16.4.2.40 뇨검사",style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.34To16.4.2.46[[7]])
# doc <- body_add_break(doc)
# doc <- body_add_par(doc, "16.4.2.41 뇨검사",style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.34To16.4.2.46[[8]])
# doc <- body_add_break(doc)
# doc <- body_add_par(doc, "16.4.2.42 뇨검사",style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.34To16.4.2.46[[9]])
# doc <- body_add_break(doc)
# doc <- body_add_par(doc, "16.4.2.43 뇨검사",style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.34To16.4.2.46[[10]])
# doc <- body_add_break(doc)
# doc <- body_add_par(doc, "16.4.2.44 뇨검사",style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.34To16.4.2.46[[11]])
# doc <- body_add_break(doc)
# doc <- body_add_par(doc, "16.4.2.45 뇨검사",style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.34To16.4.2.46[[12]])
# doc <- body_add_break(doc)
# doc <- body_add_par(doc, "16.4.2.46 뇨검사",style="Normal Indent")
# doc <- body_add_flextable(doc, MyFTable_16.4.2.34To16.4.2.46[[13]])
# doc <- body_add_break(doc)
#####

doc <- body_add_par(doc, "16.4.2.4 혈액응고검사",  style="heading 3") 
doc <- body_add_flextable(doc, MyFTable_16.4.2.44)
doc <- body_add_break(doc)

doc <- body_add_par(doc, "16.4.2.5 혈청검사",  style="heading 3")
doc <- body_add_flextable(doc, MyFTable_16.4.2.45)
doc <- body_add_break(doc)

#doc <- body_add_par(doc, "16.4.2.49 소변약물 검사",  style="heading 3")
#doc <- body_add_flextable(doc, MyFTable_16.4.2.46)
#doc <- body_add_break(doc)

doc <- body_add_par(doc, "병력", style = "heading 3") 	#16.4.3
doc <- body_add_flextable(doc, MyFTable_16.4.3)
doc <- body_add_break(doc)


doc <- body_add_par(doc, "시험대상자의 선정", style = "heading 3") #16.4.4
doc <- body_add_flextable(doc, MyFTable_16.4.4)
doc <- body_add_break(doc)

doc <- body_end_section_continuous(doc) # 가로 시작
doc <- body_add_par(doc, "신체 검사", style = "heading 3") 	#16.4.5
doc <- body_add_flextable(doc, MyFTable_16.4.5.1)
doc <- body_add_break(doc)
doc <- body_add_flextable(doc, MyFTable_16.4.5.2)
doc <- body_add_break(doc)

doc <- body_add_par(doc, "활력 징후", style = "heading 3") # 16.4.6
doc <- body_add_par(doc, "16.4.6.1 활력징후[SBP]", style="heading 4")
doc <- body_add_flextable(doc, VS.SYSBP)
doc <- body_add_break(doc)

doc <- body_add_par(doc, "16.4.6.2 활력징후[DBP]", style="heading 4")
doc <- body_add_flextable(doc, VS.BPSYS)
doc <- body_add_break(doc)

doc <- body_add_par(doc, "16.4.6.3 활력징후[PULSE]", style="heading 4")
doc <- body_add_flextable(doc, VS.PL)
doc <- body_add_break(doc)

doc <- body_add_par(doc, "16.4.6.4 활력징후[BODY TEMPERATURE]", style="heading 4")
doc <- body_add_flextable(doc, VS.TEMP)
doc <- body_add_break(doc)
doc <- body_end_section_landscape(doc) # 가로 끝

doc <- body_add_par(doc, "12-Lead 심전도", style = "heading 3") # 16.4.7
doc <- body_add_par(doc, "16.4.7.1 12-Lead 심전도[Ventricular rate]", style = "heading 4") 
doc <- body_add_flextable(doc, EGVR)
doc <- body_add_break(doc)

doc <- body_add_par(doc, "16.4.7.2 12-Lead 심전도[PR Interval]", style = "heading 4") 
doc <- body_add_flextable(doc, EGPR)
doc <- body_add_break(doc)

doc <- body_add_par(doc, "16.4.7.3 12-Lead 심전도[QRSD Interval]", style = "heading 4") 
doc <- body_add_flextable(doc, EGQRSD)
doc <- body_add_break(doc)

doc <- body_add_par(doc, "16.4.7.4 12-Lead 심전도[QT Interval]", style = "heading 4") 
doc <- body_add_flextable(doc, EGQT)
doc <- body_add_break(doc)

doc <- body_add_par(doc, "16.4.7.5 12-Lead 심전도[QTc Interval]", style = "heading 4") 
doc <- body_add_flextable(doc, EGQTC)
doc <- body_add_break(doc)

doc <- body_add_par(doc, "16.4.7.6 12-Lead 심전도[Over all]", style = "heading 4") 
doc <- body_add_flextable(doc, ECGNORM)
doc <- body_add_break(doc)
doc <- body_end_section_landscape(doc) #끝

doc <- body_add_par(doc, "혈당검사", style = "heading 3") # 16.4.8.1
doc <- body_add_par(doc, "16.4.8.1 혈당검사", style = "heading 4") 
doc <- body_add_flextable(doc, BST)
doc <- body_add_break(doc)


appendix.name <- strsplit(getwd(), "/")[[1]][6]
print(doc, target = paste0(dirname(getwd()), "/", appendix.name, "/",appendix.name,"_Appendix_plus.docx"))
