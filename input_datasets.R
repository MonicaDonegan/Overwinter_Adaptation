library(dplyr)
library(tidyverse)
library(data.table)

setwd("/Users/monicadonegan/Downloads/TRADEOFF_MANUSCRIPT/HREC_analysis")
source("qpcr_processing_adaptation.R")

## input qpcr dataset
all_data<- input_qpcr_files("data_files/qpcr_record.csv", "qpcr_files")
all_data<- merge(all_data, read.csv('../data_files/inoculation_treatment.csv'), by = c("Row", 'Vine'))

summary<- all_data %>%
  group_by(Sample.Name,Sample.Date, Variety, Treatment, Row, Vine, Side, Run.Date) %>%
  summarise(mean = mean(Ct), n=n(), sd = sd(Ct)) %>%
  determine_positives(37) %>% 
  mutate(VineID = paste(Row, Vine))

## add culturing data 
culturing<- read.csv('../data_files/culturing_results.csv') %>% 
  mutate(VineID = paste(Row, Vine)) %>% 
  filter(Response == 1)

culturing$year<- format(as.Date(culturing$Sampled, format="%m/%d/%y"),"%Y")
culturing$month<- format(as.Date(culturing$Sampled, format="%m/%d/%y"),"%m")

##group by year
sum_byvine_2021<- sum_by_vine_general(summary, 2021)
sum_byvine_2021<- check_culture(sum_byvine_2021, culturing, summary)

sum_byvine_2022<- sum_by_vine_general(summary, 2022)

sum_byvine_2023<- sum_by_vine_general(summary, 2023)
sum_byvine_2023<- check_culture(sum_byvine_2023, culturing, summary)

overwinter<- merge(sum_byvine_2021, 
                   sum_byvine_2022,
                   by = c('Variety', 'Treatment', 'VineID', 'Row', 'Vine'), all.x = T,all.y = T) %>%
  merge(sum_byvine_2023, by = c('Variety', 'Treatment', 'VineID', 'Row', 'Vine'), all.x = T)

overwinter<- overwinter %>%
  dplyr::select("Variety", "Treatment", "VineID", "Row", "Vine", "pos_any.x", "pos_any.y", "pos_any") %>%
  rename(Pos2021 = "pos_any.x", Pos2022 = "pos_any.y", Pos2023= "pos_any")


### fix overwinter 
## if sample was missing for 1 year (NA) but other 2 years were negative - assumed negative
## exception - 4-9 was pos in 2021 and 2022, no sample in 2023 because dead 
overwinter$Pos2021[is.na(overwinter$Pos2021)] <- 0 
overwinter$Pos2022[is.na(overwinter$Pos2022)] <- 0 
overwinter$Pos2023[overwinter$VineID== '4 9'] <- 1
overwinter$Pos2023[is.na(overwinter$Pos2023)] <- 0 

overwinter$Recovery_22 <- 0
overwinter$Recovery_23 <- 0 

# determine recovery, put 2 for false negative
for(i in 1:dim(overwinter)[1]) {
  if (overwinter$Pos2021[i]==1 & overwinter$Pos2022[i]==0 & is.na(overwinter$Pos2023[i]) ) {
    overwinter$Recovery_22[i] =1 
    overwinter$Recovery_23[i] = 'rec22'
  } else if (overwinter$Pos2021[i]== 0 & overwinter$Pos2022[i] == 0 & is.na(overwinter$Pos2023[i]) ) {
    overwinter$Recovery_22[i] = 'never'
    overwinter$Recovery_23[i] = 'never'
  } else if(is.na(overwinter$Pos2023[i])) {
    overwinter$Recovery_23[i] = NA
  } else if (overwinter$Pos2021[i]== 0 & overwinter$Pos2022[i] == 0 & overwinter$Pos2023[i]==0) {
    overwinter$Recovery_22[i] = 'never'
    overwinter$Recovery_23[i] = 'never'
  } else if (overwinter$Pos2021[i] == 1 & overwinter$Pos2022[i] == 0 & overwinter$Pos2023[i]== 0) {
    overwinter$Recovery_22[i] = 1
    overwinter$Recovery_23[i] = 'rec22'
  } else if (overwinter$Pos2022[i]==1 & overwinter$Pos2023[i]==0) {
    overwinter$Recovery_23[i] = 1
    if (overwinter$Pos2021[i]==0 & overwinter$Pos2022[i]==1) {
     overwinter$Pos2021[i] = 2
    } 
  } else if (overwinter$Pos2021[i]==0 & overwinter$Pos2022[i]==1) {
    overwinter$Pos2021[i] <- 2
  } else if (overwinter$Pos2022[i]==0 & overwinter$Pos2023[i]==1 & overwinter$Pos2021[i]==0) {
    overwinter$Pos2022[i] = 2
    overwinter$Pos2021[i] = 2
  } else if (overwinter$Pos2022[i]==0 & overwinter$Pos2023[i]==1) {
    overwinter$Pos2022[i] = 2
  } else {
    print('leaving a 0')
  }
} 

overwinter_fn<- overwinter 
false_negatives<- overwinter_fn %>% filter(Pos2021 == 2 | Pos2022 == 2  ) %>% group_by(Variety, Pos2022, Pos2021) %>% summarize(n=n()) %>% pivot_longer(cols = c(Pos2021, Pos2022)) %>% filter(value ==2)# for graphing false negatives

overwinter<- overwinter %>% 
  mutate(Pos2021 = ifelse(Pos2021>0, 1, 0),
         Pos2022 = ifelse(Pos2022>0, 1, 0))

summary$Treatment<- as.factor(summary$Treatment)
summary$Treatment <- relevel(summary$Treatment, 'Hopland')

## summarize by year, strain 
summarize_allvar_tr<- overwinter %>%
  pivot_longer(Pos2021:Pos2023, values_to = 'pos', names_to = 'year') %>%
  mutate(year = gsub('Pos', '', year)) %>% 
  group_by(year, Treatment) %>%  
  dplyr::summarize(infected = sum(pos ==1 , na.rm =T), n = n()) %>%
  mutate(percent_infected = infected / n * 100) 

## summarize by year, strain, variety
summarize_tr_var<- overwinter %>%
  pivot_longer(Pos2021:Pos2023, values_to = 'pos', names_to = 'year') %>%
  mutate(year = gsub('Pos', '', year)) %>% 
  group_by(year, Treatment, Variety) %>%  
  dplyr::summarize(infected = sum(pos ==1 , na.rm =T), n = n()) %>%
  mutate(percent_infected = infected / n * 100) 

## filter summary for vines that were tested in july 2022 
vine_id_22<- summary %>% filter(Sample.Date == '7/7/22')
summary_filter<- summary %>% filter(VineID %in% vine_id_22$VineID, Side==1)

## movement
movement <- summary %>%
  group_by(year,Side, Variety, month) %>%  
  summarize(infected = sum(Status == 'Xf+'), n = n()) %>%
  mutate(percent_infected = infected / n * 100) %>% 
  filter( Variety %in% c('Albarino', 'Tinta Francisca', 'Tannat')) %>%
  mutate(Variety = factor(Variety, levels = c('Tannat', 'Tinta Francisca', 'Albarino'))) %>%
  mutate(Side = case_when(
    Side == 1 ~ 'W', 
    Side == 2 ~ 'WC',
    Side == 3 ~ 'EC', 
    Side == 4 ~ 'E'
  )) %>% mutate(Side = factor(Side, levels = c('W', 'WC', 'EC', 'E')))

## grouped by strain
movement_tr <- summary %>%
  group_by(year,Side, Variety, month, Treatment) %>%  
  summarize(infected = sum(Status == 'Xf+'), n = n()) %>%
  mutate(percent_infected = infected / n * 100) %>% 
  filter( Variety %in% c('Albarino', 'Tinta Francisca', 'Tannat')) %>%
  mutate(Side = case_when(
    Side == 1 ~ 'W', 
    Side == 2 ~ 'WC',
    Side == 3 ~ 'EC', 
    Side == 4 ~ 'E'
  )) %>% mutate(Side = factor(Side, levels = c('W', 'WC', 'EC', 'E')))

## model of movement 
movement_full<- summary %>% 
  filter( Treatment != 'control', Variety %in% c('Albarino', 'Tinta Francisca', 'Tannat'), Side != 'X', month == '09') 

## add rows for inferred negatives on other sides
movement_full_sum<- movement_full %>% group_by(VineID, Treatment, Variety) %>% summarize(n=n())
mock<- data.frame(Treatment = c(rep('Bakersfield', 12), rep('Hopland', 12)), year = (rep(c(rep(2021, 4), rep(2022, 4), rep(2023,4)), 2)), Side = rep(c(1,2,3,4), 6))
movement_full_sum<- merge(movement_full_sum, mock, by = 'Treatment')

movement_full_sum<- merge(movement_full_sum, movement_full, by =c('Treatment', 'year', 'Variety', 'Side', "VineID"), all.x =T) 
movement_full_sum$Response[is.na(movement_full_sum$Response)]<- 0

## Susceptibility - year 1 positives 
infection_by_var<- overwinter %>% filter(Treatment != 'control') %>% group_by(Variety) %>% summarize( pos = sum(Pos2021), n = n(), perc = pos / n * 100) 

year1_no_control<- overwinter %>%
  filter(Treatment != 'control')

## ct values pre winter 
all_nocontrols<- summary %>%
  filter(Treatment != 'control', Side == 1, year == 2021, month == '09', Response == 1)

## log data is approximately normal 
all_nocontrols %>% ggplot(aes(log)) + geom_histogram(bins = 20)

sum_ct<- all_nocontrols %>% group_by(Variety) %>% summarize(log_ct = mean(log)) 
sum_ct<- merge(sum_ct, infection_by_var, by = "Variety")  
sum_ct$perc<- sum_ct$perc / 100

## post winter model of population size
## only use vines from 2022 that had been sampled at all three months 
all_nocontrols_2022<- summary_filter %>%
  filter(Treatment != 'control', Side == 1, year ==2022, Response==1)

## plotting overwinter
overwinter_sum22_tr<- overwinter %>%
  filter(Recovery_22 != 'never', ! is.na(Recovery_22), Treatment!= 'control') %>%
  group_by(Treatment) %>% # add var / treat
  dplyr::summarize(recovered = sum(Recovery_22 == 1), n = n()) %>%
  mutate(survival_perc = 100 - recovered / n * 100)%>%
  mutate(year = '2022')

overwinter_sum23_tr<- overwinter %>%
  filter(Recovery_23 != 'never', ! is.na(Recovery_23), Recovery_23!= 'rec22',  Treatment!= 'control') %>%
  group_by(Treatment) %>%  ## add var
  dplyr::summarize(recovered = sum(Recovery_23 == 1), n = n()) %>%
  mutate(survival_perc = 100 - recovered / n * 100) %>%
  mutate(year = '2023')

overwinter_sum_tr<- rbind(overwinter_sum22_tr, overwinter_sum23_tr) %>%
    mutate(infected = n - recovered) 

overwinter_sum22_var<- overwinter %>%
  filter(Recovery_22 != 'never', ! is.na(Recovery_22), Treatment!= 'control') %>%
  group_by(Variety) %>% 
  dplyr::summarize(recovered = sum(Recovery_22 == 1), n = n()) %>%
  mutate(survival_perc = 100 - recovered / n * 100)%>%
  mutate(year = '2022')

overwinter_sum23_var<- overwinter %>%
  filter(Recovery_23 != 'never', ! is.na(Recovery_23), Recovery_23!= 'rec22',  Treatment!= 'control') %>%
  group_by(Variety) %>% 
  dplyr::summarize(recovered = sum(Recovery_23 == 1), n = n()) %>%
  mutate(survival_perc = 100 - recovered / n * 100) %>%
  mutate(year = '2023')

overwinter_sum_var<- rbind(overwinter_sum22_var, overwinter_sum23_var) %>%
  mutate(infected = n - recovered) 

overwinter_sum_var$Variety<- factor(overwinter_sum_var$Variety, level = c("Tempranillo", "Tinta Amarella", "Sagrantino","Mencia", "Albarino", "Petit Manseng", "Greco di Tufo","Ciliegiolo", "Periquita","Teroldego","Tinta Francisca", "Falanghina",  "Tannat"))

overwinter_sum_total<- overwinter %>% 
  group_by(Variety) %>% 
  filter(Treatment != 'control') %>% 
  dplyr::summarize(pos_2021 = sum(Pos2021 == 1), pos_2023 = sum(Pos2023==1), n = n()) %>%
  mutate(survival_perc =  (pos_2023 / pos_2021 * 100))

## correlate this with susceptibility 
sum_ct <- merge(overwinter_sum_total, sum_ct, by = "Variety")

season_21<- summary %>% 
  filter(year == '2021') %>%
  group_by(VineID, Treatment, Variety) %>% 
  arrange(desc(month)) %>% 
  summarize(log_min = max(log))

season_22<- summary %>% 
  filter(year == '2022') %>%
  group_by(VineID, Treatment, Variety) %>%
  arrange(desc(month)) %>% 
  summarize(log_min = max(log))

overwinter_withCT<- merge(overwinter, season_21, by=c('Treatment', 'Variety', 'VineID'), all.x = T)
overwinter_withCT<- merge(overwinter_withCT, season_22, by=c('Treatment', 'Variety', 'VineID'), all.x = T)

overwinter_withCT_long<- overwinter_withCT %>% 
  dplyr::select(Treatment, Variety, VineID, Recovery_22, Recovery_23) %>% 
  pivot_longer( cols=4:5) 

overwinter_withCT_long<- merge(overwinter_withCT_long, overwinter_withCT, by= c("Treatment", "Variety", "VineID"))

overwinter_withCT_long<- overwinter_withCT_long %>%
  mutate(log = case_when(name == 'Recovery_22' ~ log_min.x, 
            name == 'Recovery_23' ~ log_min.y)) %>%
  dplyr::select(Treatment, Variety, log, value, name, VineID)

overwinterCT_excludeNever<- overwinter_withCT_long %>%
  filter(! value %in% c('never', 'rec22')) %>%
  filter(Treatment != 'control') %>%
  filter(log != 0 ) %>%
  mutate(year = name)

overwinterCT_excludeNever <- overwinterCT_excludeNever %>% 
  mutate(year_match = ifelse(name == 'Recovery_22', 2021, 2022)) %>% 
  mutate(Variety = factor(Variety), Variety = relevel(Variety, ref = "Ciliegiolo"))

# overwinter data for symptoms 
overwinter_corrected <- overwinter  %>% pivot_longer( cols=6:8, names_to= "name")  %>% 
  mutate(year = str_remove(name, "Pos"))

###symptoms
symptom<- read.csv('../data_files/Symptom_2021.csv')
symptom_aug22<- read.csv('../data_files/August_2022_symptoms.csv')
symptom_sept22<- read.csv('../data_files/September_2022_symptoms.csv')
symptom_2023<- read.csv('../data_files/Symptoms2023.csv')
symptom_2023<- symptom_2023[,-1]

symptom<- rbind(symptom, symptom_aug22, symptom_sept22) 
symptom$Dead<- 0

symptom<- rbind(symptom, symptom_2023) %>%mutate(Variety = case_when(Variety== 'Graciano'~ 'Tempranillo', TRUE ~ Variety))

symptom <- symptom %>%
  mutate(Direction = case_when(
    Direction %in% c('West', 'W') ~ 'West', 
    Direction %in% c('East', 'E') ~ 'East', 
    TRUE~ Direction
  ))

symptom$year<- format(as.Date(symptom$Sample.Date, format="%m/%d/%y"),"%Y")
symptom$month<- format(as.Date(symptom$Sample.Date, format="%m/%d/%y"),"%m")
symptom$VineID <- paste(symptom$Row, symptom$Vine)

##group by vine
symp_byvine<- symptom %>%
  filter(!is.na(Variety)) %>%
  group_by(year, Variety, Treatment, VineID, Direction) %>%
  summarise(no_symptom = sum(None, na.rm = T), n = n()) %>%
  mutate(symptomatic = ifelse(no_symptom ==n, 0, 1))

symp_byvine<- merge(symp_byvine, overwinter_corrected, by= c("Treatment", "Variety", "year", "VineID"), all.y = T)
symptom<- merge(symptom, overwinter_corrected, by= c("Treatment", "Variety", "year", "VineID"), all.y = T)

symptom_summary<- symptom %>%
  filter(Treatment != 'control') %>%
  filter(value == 1) %>%
  group_by(year, month, Direction, Treatment) %>% 
  summarise(no_symptom = sum(None), 
            bar_symptom = sum(BarrenShoot), 
            scorch_symptom = sum(Scorch), 
            match_symptom= sum(Match), 
            stunt_symptom = sum(Stunt), 
            SC_symptom = sum(SC), 
            UnevenL_symptom = sum(UnevenL),
            dead_symptom = sum(Dead),
            n = n()) %>%
  mutate(symptomatic = n -no_symptom, 
         percent_symptomatic = (symptomatic / n) * 100, 
        BarrenShoot = (bar_symptom / n) * 100, 
        Match = (match_symptom / n) * 100, 
        Scorch = (scorch_symptom / n) * 100, 
        Stunt = (stunt_symptom / n) * 100, 
        SC = (SC_symptom / n) * 100, 
        UnevenL = (UnevenL_symptom / n) * 100, 
        Dead = (dead_symptom / n) * 100) 

melted_symp<- melt(as.data.table(symptom_summary), id.vars = c('month', 'year', 'Direction', 'Treatment'), measure.vars = c('BarrenShoot', 'Match', 'Scorch', 'Stunt', 'SC', 'UnevenL', 'Dead'))

##symptomatic at either point
symp_sum<- symp_byvine %>%
  filter(value == 1) %>% 
   filter(Treatment != 'control') %>%
  group_by(year, Treatment, Direction) %>% 
  summarise(symptom_perc = round(sum(symptomatic / n()), 2), count= sum(symptomatic), n = n())

symptom_summary$Direction<- factor(symptom_summary$Direction, levels = c('West', 'WC', 'EC', 'East'))

symp_sum<- symp_sum %>% 
  mutate(Direction = case_when(
    Direction == 'West' ~ 'W', 
    Direction == 'East' ~ 'E', 
    TRUE~ Direction
  ))
symp_sum$Direction<- factor(symp_sum$Direction, levels = c('W', 'WC', 'EC', 'E'))

symp_sum_var <- symp_byvine %>%
  filter(value == 1) %>% 
  filter(Treatment != 'control') %>%
  group_by(year, Variety, Direction) %>% ## change to var
  summarise(symptom_perc = round(sum(symptomatic / n()), 2), count= sum(symptomatic), n = n())

## symptoms and recovery 
symp_sum_tot<- symp_byvine %>%
  filter(value == 1) %>% 
  filter(Treatment != 'control') %>%
  group_by(year, Variety, Direction) %>% #
  summarise(symptom_perc = round(sum(symptomatic / n()), 2), count= sum(symptomatic), n = n())%>% filter(year == 2021, Direction== 'West')

sum_ct <- merge(sum_ct, symp_sum_tot, by = "Variety")

## symptoms in year 1, grouped by month 
symptom_treated_2021<- symptom %>% filter(Treatment != 'control', value == 1, year ==2021) %>% group_by(VineID, Variety, Treatment) %>% summarize(symptoms_tot = sum(None), symptomatic = ifelse(symptoms_tot == 8, 0, 1))

## only positive vines - all years 
symptom_treated<- symptom %>% filter(Treatment != 'control', value == 1) %>% mutate(Direction = factor(Direction, levels = c("West", "WC", "EC", "East")), symptomatic = 1 - None)

## import stunting data 
stunt_2022 <- read.csv('../data_files/stunting_2022.csv') %>% 
  mutate(year = 2022)
stunt_april2023<- read.csv('../data_files/stunting_april2023.csv') %>% 
  mutate(year = 2023)
stunt_may2023<- read.csv('../data_files/stunting_may2023.csv') %>% 
  mutate(year = 2023)

stunt<- rbind(stunt_2022, stunt_april2023, stunt_may2023)

stunt_average<- stunt %>% 
  group_by(Row, Vine, Variety, `flag.color`, ind, year, month) %>%
  summarize(ave_height = mean(values), sd_height= sd(values), count= n())

stunt_merge<- merge(stunt_average, overwinter_corrected, by = c("Row", "Vine", "Variety", "year"), all.x = T )

stunt_merge_ave <- stunt_merge %>% 
  filter(value== 0) %>% 
  group_by(Variety, ind, month, year) %>% 
  summarize(average_uninfect = mean(ave_height, na.rm =T), count = n())

stunt_merge<- merge(stunt_merge, stunt_merge_ave, by= c("Variety", 'ind', 'month', 'year'), all.x = T) %>% mutate(normalized_height = ave_height / average_uninfect)

stunt_merge %>% filter(! is.na(Treatment)) %>%
  filter(Treatment != 'control') %>%
  filter(ind == 'W') %>% ggplot(aes(normalized_height)) + geom_histogram() + facet_wrap(month~year)

## data set for modeling
stunt_filter<- stunt_merge %>%
  filter(Treatment != 'control', ind == 'W')

### data of chill hours from Sanel Valley site
hrec_chill<- read.csv('../data_files/chill_hours.csv')

hrec_chill$Date <- as.Date(hrec_chill$Date, format = "%m/%d/%y")
hrec_chill$Hours<- as.numeric(hrec_chill$Hours)

## data from iButton loggers 
temps_bk <- read.csv("../data_files/Bakersfield_temps_22_23.csv") 
temps_hop <- read.csv('../data_files/Hopland_temps_22_23.csv') 
temps_hop<- temps_hop %>% select(! Index) %>% mutate(location = 'Hopland')
temps_bk <- temps_bk[30:2077, -3] %>% mutate(location = 'Bakersfield')
colnames(temps_bk)<- c('Date', 'Time', 'Temp', 'location')
colnames(temps_hop)<- c('Date', 'Time', 'Temp', 'location')

all_temps<- rbind(temps_bk, temps_hop) %>% 
  mutate(chill_hours_7 = ifelse(Temp < 7.22, 2, 0), 
         chill_hours_0 = ifelse(Temp< 0, 2, 0))

all_temps$Date <- as.Date(all_temps$Date, format = "%m/%d/%y")

chill_sum <- all_temps %>% 
  group_by(location, Date) %>%
  summarize(chill_sum_7 = sum(chill_hours_7), 
            chill_sum_0 = sum(chill_hours_0)) %>%
  mutate(totalhours_0 = base::cumsum(chill_sum_0), 
         totalhours_7 = base::cumsum(chill_sum_7)) %>% 
  pivot_longer(cols = c(totalhours_0, totalhours_7))

