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


dapa %>%
    ggplot(aes(x = Time, y = Conc, col = Period)) +
    geom_line()

# DB 수령 후 Real-time 적용 예정

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
dapa_NCA <- tblNCA(new_dapa, key = c("Period", "ID"), colTime = "RTime", colConc = "Conc", dose = 100, adm = "Extravascular", R2ADJ = -1, concUnit = 'ng/mL') %>%
    select(Period, ID, CMAX, TMAX, LAMZHL, AUCLST, AUCIFO, CLFO, VZFO) 


## Create summary table 
dapa_NCA %>%
    select(-ID) %>%
    tbl_summary(by = Period, 
                type = TMAX ~ "continuous", 
                statistic = c("TMAX") ~ "{median} ({min}- {max})")

## Comparative PK
dapa_BE_raw <- dapa_NCA  %>%
    mutate(LCMAX = log(CMAX), LAUCLST = log(AUCLST))

f <- LCMAX ~ Period  # LCMAX, LAUCLST
# f <- LAUCLST ~ Period 

BE <- lme(f, random = ~1|ID, data = dapa_BE_raw)    
ci <- intervals(BE, 0.9)
exp(ci$fixed["Period2기", ])   ## 90% CI result  

head(dapa_BE_raw) %>% write.csv('data.csv', row.names = F)
GLM(f, dapa_BE_raw)$ANOVA %>% round(4)    ## Anova result

