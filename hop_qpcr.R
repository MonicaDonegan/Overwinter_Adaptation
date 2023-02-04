
library(dplyr)

treatments<- read.csv('/Users/monicadonegan/Downloads/inoculation_treatment.csv')

qpcr_setup<- read.csv('/Users/monicadonegan/Downloads/qpcr_record.csv')

metadata<- merge(qpcr_setup, treatments, by = c("Row", 'Vine'), all.x = TRUE)
metadata$Run.Date<- as.Date(metadata$Run.Date, format="%m/%d/%y")

setwd("/Users/monicadonegan/Downloads/hopland_qpcr")

library(stringr)
list<- list.files()
outputs<- c()

for (i in 1:length(list)) {
  name <- list[i]
  name<- str_remove(name, '.csv')
  print(name)
  plate<- read.csv(list[i])
  datestring<- as.POSIXct(str_split(name, "_")[[1]][3],format="%m-%d-%y")
  #print(plate[29,1])
  if (plate[28,1] == 'Well') {
    colnames(plate)<- plate[28,]
    plate<- plate[29:124,]

  } else if (plate[29,1] == 'Well') {
    colnames(plate)<- plate[29,]
    plate<- plate[30:125,]
  }
  colnames(plate)[2]<- 'Sample.Name'
  print(datestring)
  plate$Run.Date<- as.Date(datestring)
  print(plate$Run.Date)
  plate$Ct[plate$Ct == 'Undetermined'] <- 40
  plate$Ct <- as.numeric(plate$Ct)
  name<- str_split(name, "_")[[1]][2]
  assign(name, plate)
  outputs<- c(outputs, plate)
}

colnames(plate30) <- colnames(plate29)
all_data<-rbind(plate01,plate02,plate03,plate04, plate05,
                plate06,plate07, plate08, plate09, plate10, 
                plate11, plate12, plate13, plate14, plate15, 
                plate16, plate18, plate19, plate20, plate21,
                plate22, plate24, plate25, plate26,
                plate27, plate28, plate29, plate30, plate31,
                plate32, plate33, plate34, plate35,
                plate36, plate37, plate38, plate39, 
                plate40, plate41, plate42, plate43, plate44)

all_data<- merge(all_data, metadata, by=c('Run.Date', 'Sample.Name'))

summary<- all_data %>%
  group_by(Sample.Name,Sample.Date, Variety, Treatment, Row, Vine, Side, Run.Date) %>%
  summarise(mean = mean(Ct), n=n(), sd = sd(Ct))

large_st<- summary %>%
  filter(sd > 1 && mean < 37)

index_list <- c()
for(i in seq(1:dim(all_data)[1])) {
  slice<- all_data[i,] 
  for(j in seq(1:dim(large_st)[1])) {
    if(slice$Row == large_st$Row[j] & 
       slice$Vine == large_st$Vine[j] &
       slice$Side == large_st$Side[j] & 
       slice$Run.Date == large_st$Run.Date[j]) {
      print('HIT')
      print(i)
      if(slice$Ct == 40) {
        print('this is 40')
        index_list<- c(index_list, i)
      }
    } 
  }
}

all_data_2<- all_data[-index_list,]         

summary<- all_data_2 %>%
  group_by(Sample.Name,Sample.Date, Variety, Treatment, Row, Vine, Side, Run.Date) %>%
  summarise(mean = mean(Ct), n=n(), sd = sd(Ct))


summary$CFUperG <- (10 ^ ((log10(summary$mean) - 1.54) / -0.0454) ) *62500
### add an infected column 
summary$Status<- ''
summary$Response <- 0
for (i in 1:dim(summary)[1]) {
  if (summary$mean[i] < 36.5 ) {
    summary$Status[i]<- 'Xf+'
    summary$Response[i]<- 1
  } else {
    summary$Status[i]<- '(-)'
  }
}

summary$log <- log10(summary$CFUperG)
summary$log[summary$Status != 'Xf+'] <- 0 

summary$year<- format(as.Date(summary$Sample.Date, format="%m/%d/%y"),"%Y")
summary$month<- format(as.Date(summary$Sample.Date, format="%m/%d/%y"),"%m")
summary$VineID<- paste(summary$Row, summary$Vine)


movement<- summary %>% 
  filter(Variety %in% c('Tannat', 'Albarino', 'Tinta Francisca')) %>%
  filter(year == 2022) %>%
  group_by( month, Side, Variety) %>%
  summarize(count_pos = sum(Response), n = n())

movement<- summary %>% 
  filter(year == 2021, month == '09') %>%
  group_by( Side) %>%
  summarize(count_pos = sum(Response), n = n())

library(viridis)
library(ggpubr)
BlueRed <- colorRampPalette(c("slateblue2", "tomato1"), bias=1, space="rgb", interpolate="linear") 


summarize_allvar<- summary %>%
  group_by(year, month, Treatment, Side) %>% 
  summarize(infected = sum(Status == 'Xf+'), n = n(), ctmean = sum(CFUperG)) %>%
  mutate(percent_infected = infected / n * 100, ct = ctmean / n)
summarize_allvar$Treatment<- as.factor(summarize_allvar$Treatment)
summarize_allvar$Treatment <- relevel(summarize_allvar$Treatment, 'Hopland')
barplot_all<- summarize_allvar %>%
  filter(Treatment %in% c('Bakersfield', 'Hopland')) %>%
  filter(Side == 1) %>%
  #filter(year == 2021) %>%
  ggplot( aes(x=Treatment, y=infected, fill = Treatment)) +
  xlab('PD Strain') + 
  geom_col() +
  ylab('Count of infected vines') + 
  ylim(0.0, 100.0)+
  geom_text(aes(label = paste(infected, "/", round(n)), 
                vjust = -.25)) + 
  scale_colour_manual(aesthetics = c("colour", "fill"),values = BlueRed(2))+  
  facet_grid(month~year)+
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  )  + 
  theme_bw() + 
  guides(fill = "none", size = "none")

summary$Treatment<- as.factor(summary$Treatment)
summary$Treatment <- relevel(summary$Treatment, 'Hopland')

bp_log<- summary %>%
  #filter(Sample.Date == '6/7/21') %>%
  filter(Status == 'Xf+') %>%
  filter(year == 2022) %>%
  filter(month != '06') %>%
  filter(Treatment %in% c('Bakersfield', 'Hopland')) %>%
  ggplot( aes(x=Treatment, y=log, fill = Treatment)) +
  geom_boxplot(show.legend = F) +
  geom_point(show.legend = F) + 
  xlab('PD Strain') + 
  ylab('log(CFU/g)') + 
  facet_wrap(~month, ncol =1)+
  stat_compare_means(method = "wilcox.test", label.y=8) +
  #scale_colour_manual(aesthetics = c("colour", "fill"),values = BlueRed(2))+  
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_point(color="black", size=0.4, alpha=0.9) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) + 
  theme_classic() + 
  guides(fill = "none")

ggarrange(barplot_all, bp_log, 
          labels = c("A", "B"),
          widths = c(1.5,1),
          ncol = 2)


### look at rates in Var, Treatment
summarize_var<- summary %>%
  group_by(year, month, Variety, Treatment, Side) %>% 
  summarize(infected = sum(Status == 'Xf+'), n = n(), ctmean = sum(CFUperG)) %>%
  mutate(percent_infected = infected / n * 100, ct = ctmean / n)


library(ggplot2)
summarize_var %>%
  filter(year == 2021) %>%
  filter(month == '09') %>%
  filter(Side ==1 ) %>%
  filter(Variety %in% c('Albarino', 'Tinta Francisca', 'Tannat', 'Teroldego')) %>%
  filter(Treatment %in% c('Bakersfield', 'Hopland')) %>%
  ggplot( aes(x=Treatment, y=percent_infected, fill = Treatment)) +
  geom_col() +
  ylim(0,110) + 
  facet_wrap(~Variety, ncol=2)+
  geom_text(aes(label = paste(infected, "/", round(n)), 
                vjust = -.25,  size=1.5)) + 
  #stat_compare_means(method = "t.test", label.y=10) +
  scale_colour_manual(aesthetics = c("colour", "fill"),values = BlueRed(2))+  
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  )

sept21_sum <- summary  %>%
  filter(year == '2021', month == '09', Side == '1', 
         Treatment != 'control') %>%
  group_by(year, month, Variety, Side) %>% 
  summarize(infected = sum(Status == 'Xf+'), n = n()) %>%
  mutate(percent_infected = round(infected / n * 100, 1))



ggbarplot(filtered, x= "Treatment", y = "percent_infected",  ylab = "Percent Positive", fill ="Treatment", palette = BlueRed(2)) + 
  facet_grid(Variety~Sample.Date, switch = 'y') + guides(fill=guide_legend(title= "")) + ylim(0,100)+
  theme(text = element_text(size = 10), legend.position = "none") 


filtered %>%
  ggplot(aes(Treatment, percent_infected, fill = Treatment)) + 
  geom_col()+
  scale_colour_manual(aesthetics = c("colour", "fill"),values = BlueRed(2)) +
  facet_wrap(~Variety, ncol=3) + 
  ylim(0.0, 110.0)+
  ylab('Percent Positive')+
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = rel(1), color = 'black'),
    strip.background = element_rect(fill = "lightgray", colour = "black", size = 1)
  ) 
  theme(axis.ticks = element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  theme( legend.position = "none",axis.text.x = element_text(angle = 60,  hjust=1))+ 
  #geom_text(aes(label = round(percent_infected)), vjust = -0.05,  size=2.5) + 
 # theme(strip.text.y.left = element_text(angle = 0)) + 
filtered %>%
  ggplot(aes(Variety, percent_infected, fill = Variety)) + 
  geom_col()+
  scale_colour_manual(aesthetics = c("colour", "fill"),values = BlueRed(13)) +
  facet_wrap(~Treatment, nrow=2) + 
  ylim(0.0, 110.0)+
  ylab('Percent Positive')+
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = rel(1), color = 'black'),
    strip.background = element_rect(fill = "lightgray", colour = "black", size = 1)
  ) 
  
    
#budbreak<- read.csv('~/Downloads/budbreak.csv')
#summary<- merge(summary, budbreak, by = 'Variety', all.x = T)
#BudBreak - kinda significant 
  
##test with binomial glm 
library(glmmTMB)
alldata_2021<- summary[summary$year == '2021', ]
alldata_2021<- alldata_2021[alldata_2021$Treatment != 'control', ]
alldata_2021<- alldata_2021[alldata_2021$Side == '1', ]

all_data_side1<- summary %>%
  filter(Treatment != 'control' & Side == 1)

glmresult<- glmmTMB(Response~ Treatment + Variety + year + month + (1 | VineID),
    family = 'binomial', data = all_data_side1)

summary(glmresult)
library(car)
Anova(glmresult, test = "Chisq")

alldata_2021 <- alldata_2021[alldata_2021$log > 0, ]
logglm<- glmmTMB(log~ Treatment + Variety + month + (1 | VineID),
                  data = alldata_2021)

## Vine ID  as Random Effect 
library(lme4)
glmmres<- glmer(Response~Treatment + Variety + month  + year  +
      (1 | VineID), family='binomial', data=all_data_side1)
summary(glmmres)
anova(glmmres)

## chi square of Variety
counts<- summary %>%
  filter(Treatment != 'control') %>%
  filter(Sample.Date == '7/7/21') %>%
  group_by(Sample.Date, Variety) %>% 
  #group_by(Sample.Date, Treatment) %>% 
  summarize(infected = sum(Status == 'Xf+'), n = n()) %>%
  mutate(infected = infected, uninfected = n - infected)

vars<- counts$Variety
counts<- counts[,c(3,5)]
rownames(counts)<- vars

chisq.test(counts, correct = FALSE)

###symptoms

symptom<- read.csv('/Users/monicadonegan/Downloads/Symptom_2021.csv')

symptom_aug22<- read.csv('/Users/monicadonegan/Downloads/August19_22_symptoms.csv')
symptom_sept22<- read.csv('/Users/monicadonegan/Downloads/Sept29_symptom_2022.csv')
symptom<- rbind(symptom, symptom_aug22, symptom_sept22)

symptom$year<- format(as.Date(symptom$Sample.Date, format="%m/%d/%y"),"%Y")
symptom$month<- format(as.Date(symptom$Sample.Date, format="%m/%d/%y"),"%m")
symptom$VineID <- paste(symptom$Row, symptom$Vine)



symptom_summary<- symptom %>%
  filter(Treatment != 'control') %>%
  group_by(year, month, Treatment, Direction) %>% 
  summarise(no_symptom = sum(None), 
            bar_symptom = sum(BarrenShoot), 
            scorch_symptom = sum(Scorch), 
            match_symptom= sum(Match), 
            stunt_symptom = sum(Stunt), 
            SC_symptom = sum(SC), 
            UnevenL_symptom = sum(UnevenL),
            n = n()) %>%
  mutate(symptomatic = n -no_symptom) %>%
  mutate(percent_symptomatic = (symptomatic / n) * 100) %>%
  mutate(BarrenShoot = (bar_symptom / n) * 100) %>%
  mutate(Match = (match_symptom / n) * 100) %>%
  mutate(Scorch = (scorch_symptom / n) * 100) %>%
  mutate(Stunt = (stunt_symptom / n) * 100) %>%
  mutate(SC = (SC_symptom / n) * 100) %>%
  mutate(UnevenL = (UnevenL_symptom / n) * 100) 

BlueRed <- colorRampPalette(c("slateblue2", "tomato1"), bias=1, space="rgb", interpolate="linear") 

library(data.table)
melted_symp<- melt(as.data.table(symptom_summary), id.vars = c('month', 'year', 'Treatment', 'Direction'), measusre.vars = c('BarrentShoot', 'Match', 'Scorch', 'Stunt', 'SC', 'UnevenL'))
perc_symp<- melted_symp[181:320,]

perc_symp$Direction<- factor(perc_symp$Direction, levels = c('West', 'WC', 'EC', 'East'))

## graph by treatment 
perc_symp %>%
  filter(variable == 'percent_symptomatic') %>%
  #filter(Direction != 'West') %>%
  filter(year == '2022') %>%
  ggplot(aes(x= month, y = value, fill = Treatment)) +
  geom_col(position = 'dodge') + 
  facet_wrap(~Direction, ncol=4) + 
  ylab('Percentage') + 
  #ylim(0, 100) + 
  scale_colour_manual(aesthetics = c("colour", "fill"),values = BlueRed(2))
                      
### graph individual symptoms                
perc_symp %>%
  filter(variable != 'percent_symptomatic') %>%
  filter(Direction == 'West') %>%
 # filter(year == '2022') %>%
  ggplot(aes(x= Treatment, y = value, fill = variable)) +
  geom_col(position = 'dodge') + 
  facet_grid(year~month) + 
  ylab('Percentage') 
  #ylim(0, 100) + 

  
library(pals)
BlueRed <- colorRampPalette(c("slateblue2", "tomato1")) 



##group by months
symp_byvine<- symptom %>%
  filter(!is.na(Variety)) %>%
  group_by(year, Variety, Treatment, VineID, Direction) %>%
  summarise(no_symptom = sum(None), n = n()) %>%
  mutate(symptomatic = ifelse(no_symptom ==2, 0, 1))

symp_movement <- symp_byvine %>%
  filter(year == 2022) %>%
  filter(Treatment != 'control') %>%
 # filter(Variety %in% c('Tannat', 'Albarino', 'Tinta Francisca')) %>%
  group_by(Variety, Treatment, Direction) %>%
  summarize(count_s = sum(symptomatic), n = n())

##symptomatic at either point
symp_sum<- symp_byvine %>%
 # filter(!is.na(Variety)) %>% 
 # filter(Treatment != 'control') %>%
  group_by(year, Variety, Treatment, Direction) %>%
  summarise(symptom_perc = round(sum(symptomatic / n()), 2), count= sum(symptomatic), n = n())

symp_sum22<- symp_sum %>%
  #filter(year =='2022') %>%
  filter(Treatment != 'control') %>%
 filter(Variety %in% c('Falanghina', 'Periquita', 'Sagrantino')) %>%
  filter(Direction == 'West') %>%
  mutate(Variety = Variety, Treatment = Treatment, 
         West = paste(symptom_perc, " (", count, "/", n, ")"))



#summary_symp$Treatment<- as.factor(summary_symp$Treatment)
#summary_symp$Treatment <- relevel(summary_symp$Treatment, 'Hopland')
summary_symp %>%
  filter(Variety != '') %>%
  filter(Direction == 'West')
  ggplot(aes(Date, percent_symptomatic, fill = Variety)) + 
  geom_col( position = 'dodge') + 
  #geom_text(aes(label = round(percent_symptomatic)), vjust = -0.5) +
  ylim(0,100) #+
 # facet_wrap(~Variety)


count_bystrain<- 
  ggbarplot(symp_sum22,  x= "year", y = "count", palette = rev(BlueRed(2)), ylab = "# Symptomatic Vines",  fill ="Treatment", legend = "right") + 
  facet_wrap(~Variety)

##stat test
glmresult<- glm(None~ Treatment + Variety + month + year + Direction, 
                family = 'binomial', data = symptom)

summary(glmresult)
anova(glmresult, test = 'Chisq')


##group qpcr data by months
sum_byvine<- summary %>%
  filter(!is.na(Treatment)) %>%
  filter(Treatment != 'Control') %>%
  #filter(year == '2022') %>%
  group_by(Variety, Treatment, VineID, year) %>%
  summarise(pos_tot = sum(Response), n = n()) %>%
  mutate(pos_any = ifelse(pos_tot >0 , 1, 0))


sum_byvine_mov<- sum_byvine %>%
  filter(Variety %in% c('Falanghina', 'Sagrantino', 'Periquita')) %>%
  filter(Treatment != 'control') %>%
  # group_by(Treatment, both, year) %>% 
  group_by(Variety, Treatment, year) %>% 
  summarize(count_pos = sum(pos_any), n = n())

overwinter_count2<- 
  ggbarplot(sum_byvine_mov,  x= "year", y = "count_pos", fill = 'Treatment', palette  = rev(BlueRed(2)),ylab = "# Pos Vines", xlab = "year", legend = 'right') + 
  facet_wrap(~Variety)

#### merge with symptoms
symp_byvine_W <- symp_byvine %>%
  filter(Direction == 'West')

all_data<- merge(sum_byvine, symp_byvine_W,by = c('VineID', "Variety", "Treatment", "year"))
all_data<- all_data %>%
  filter(!is.na(Variety)) %>%
  filter(Treatment != 'control')
all_data$both = ''
all_data$either = 0
for (i in 1:dim(all_data)[1]) {
  if (all_data$symptomatic[i] ==1 & all_data$pos_any[i] == 1) {
    all_data$both[i] =  'Symptoms & qPCR'
    all_data$either[i] = 1
  } else if (all_data$pos_any[i] == 1) {
    all_data$both[i] = 'qPCR only'
    all_data$either[i] = 1
  } else if (all_data$symptomatic[i] == 1 ) {
    all_data$both[i] = 'Symptoms only'
    all_data$either[i] = 1
  } else {
    all_data$both[i] = 'none'
    all_data$either[i] = 0
  }
}


counts<- all_data %>%
  group_by(Variety, year)  %>%
  summarise(n = n())

all_data<- merge(all_data, counts, by = c('Variety', 'year'), all.x = T) 
all_data<- all_data %>% 
  mutate(percent_positive = (either / n) * 100)

all_data_sum<- all_data %>%
  filter(both != 'none') %>%
  filter(Treatment != 'control') %>%
 # group_by(Treatment, both, year) %>% 
  group_by(Variety, both, year) %>% 
  summarize(percent_positive = sum(percent_positive), count_pos = sum(either))

#palette = c("#00AFBB", "#E7B800", "#FC4E07")
#palette = c("#E7B800", "#FC4E07")
#subset<- all_data_sum %>%
  filter(Variety %in% c('Periquita', 'Sagrantino', 'Tannat', 'Albarino'))

all_data_sum_22$both <- as.factor(all_data_sum_22$both)
all_data_sum_22$both<- relevel(all_data_sum_22$both, 'Symptoms only')

all_data_sum_22<- all_data_sum %>%
  #filter(year == '2022') %>%
  filter(Variety %in% c('Tinta Francisca', 'Tannat', 'Albarino'))
positive_all<- ggbarplot(all_data_sum_22, x= "year", y = "count_pos",  ylab = "Count Positive", fill ="both") + 
  facet_wrap(~Variety)
          #,palette = c("#FC4E07","#E7B800","#00AFBB")) + guides(fill=guide_legend(title= "")) + ylim(0,100)+
  theme(text = element_text(size = 15),legend.position = c(.25,.8))
  #scale_fill_manual(name = "", values=c("red", "black"), labels=c("A","B"))

  
## change to count 
  
  
ggarrange( overwinter_count2, count_bystrain,
            labels = c("A", "B"),
           # widths = c(1.5,1),
            ncol = 1)
  
  
qpcr_sum <- summarize_var %>%
  filter(Treatment != 'control') %>%
  filter(Sample.Date != '7/7/21')
ggbarplot(qpcr_sum, x= "Treatment", y = "percent_infected",  ylab = "Percent Positive", fill ="Treatment", palette = BlueRed(2)) + 
  facet_wrap(vars(Variety), ncol =2) + guides(fill=guide_legend(title= "")) + ylim(0,100)+
  theme(text = element_text(size = 10), legend.position = "none") 

## GLM for symptomatic 

symptom_treated<- symptom[symptom$Treatment != 'control',]
glmresult<- glm(None~ Treatment + Variety + month, 
                family = 'binomial', data = symptom_treated)
summary(glmresult)
anova(glmresult, test = 'Chisq')


glmm_res_symp<- glmer(None~Treatment + Variety + month  +  
                  (1 | VineID), family='binomial', data=symptom_treated)

summary(glmm_res_symp)
anova(glmm_res_symp)

## culturing


culturing<- read.csv('/Users/mdonegan/Downloads/HREC_Vines_culturing.csv')

culturing %>% 
  filter(plate.pos. == 1) %>%
  count(Treatment)


## overwinter recovery 

##group by months
sum_byvine_2021<- summary %>%
  filter(!is.na(Treatment)) %>%
  filter(Treatment != 'Control') %>%
  filter(year == '2021') %>%
  group_by(Variety, Treatment, VineID) %>%
  summarise(pos_tot = sum(Response), n = n()) %>%
  mutate(pos_any = ifelse(pos_tot >0 , 1, 0))

##group by months
sum_byvine_2022<- summary %>%
  filter(!is.na(Treatment)) %>%
  filter(Treatment != 'Control') %>%
  filter(year == '2022') %>%
  group_by(Variety, Treatment, VineID) %>%
  summarise(pos_tot = sum(Response), n = n()) %>%
  mutate(pos_any = ifelse(pos_tot >0 , 1, 0))

overwinter<- merge(sum_byvine_2021, 
                   sum_byvine_2022, 
                   by = c('Variety', 'Treatment', 'VineID'))

overwinter$Recovery <- overwinter$pos_any.x - overwinter$pos_any.y
overwinter$Recovery<- ifelse(overwinter$pos_any.x == 0 & overwinter$pos_any.y == 0, 'never', overwinter$Recovery)

overwinter$Recovery[overwinter$Recovery== -1]<- 0

overwinter_sum <- overwinter %>%
  filter(Recovery != 'never') %>%
  group_by(Variety, Recovery) %>% 
  summarize(count = n())

cultivar_count <- overwinter %>%
  filter(Recovery != 'never') %>%
  group_by(Variety) %>% 
  summarize(count_var = n())


result<- merge(overwinter_sum, cultivar_count, by = c('Variety'), all.x = T)
result$Percent_Recovered <- paste(round(result$count / result$count_var , 2) * 100, " (", result$count, "/", result$count_var, ")")

result_curing<- result[result$Recovery == 1,]


### join with qPCR from sept 2021
sept21<- summary %>% 
  filter(year == '2021' & month == '09' & Side ==1)

overwinter_withCT<- merge(overwinter, sept21, by=c('Treatment', 'Variety', 'VineID'), all.x = T)

overwinter_withCT %>%
  filter(Recovery != 'never') %>%
  filter(Treatment != 'control') %>%
  filter(log != 0) %>%
  ggplot(aes(Recovery, log, fill = Recovery)) + 
  geom_boxplot() + 
  geom_point() #+ 
#  facet_wrap(~Variety) # + 
#stat_compare_means(method = "wilcox.test", label.y =0.9)


overwinterCT_excludeNever<- overwinter_withCT %>%
  filter(Recovery != 'never') %>%
  filter(Treatment != 'control') %>%
  filter(log != 0 )

model<- glm(as.factor(Recovery) ~ log + Treatment + Variety, family = 'binomial' ,data = overwinterCT_excludeNever)

summary(model)
anova(model, test = 'Chisq')

Anova(model)
####
overwinter_sum <- overwinter %>%
  filter(Recovery != 'never') %>%
  group_by(Treatment, Recovery) %>% 
  summarize(count = n())
overwinter_sum$n<- c(rep(112, 4), 
                     rep(38, 3), 
                     rep(106, 4))
  
#### Raster from caladapt

#gridMET Observed Meterological Data Derived Products

library(raster)
library(sf)
library(terra)
cal20_tmin<- rast("tmin_data/tmmn_year_gridmet_2020.CA_NV.tif")
cal20_tmax<- rast("tmax_data/tmmx_year_gridmet_2020.CA_NV.tif")

calcounties<- terra::vect("R-Geospatial-Fundamentals-master/notebook_data/california_counties/CaliforniaCounties.shp")
calcounties_proj<- project(calcounties, cal20_tmin)
#calcounties_WGS = st_transform(calcounties, st_crs(cal21_proj))

plot(cal20_tmin)
plot(calcounties_proj, add = T)

plot(crop(cal20_tmin -273.15,calcounties_proj, snap = 'out', mask = T, touches =F))
plot(calcounties_proj, add = T)

#calcounties_WGS = st_transform(calcounties, st_crs(cal21_trim))

## add max temp and dots for locations
library(colorspace)
#BlueRed <- colorRampPalette(c("slateblue2", "tomato1"), bias=1, space="rgb", interpolate="linear") 
plot(crop(cal20_tmin -273.15,calcounties_proj, snap = 'out', mask = T, touches =F), axes=F,
     col=diverge_hcl(20,'Blue-Red 3'))
plot(mendocino, border = 'black', add = T)
points(pts, add = T)
plot(crop(cal20_tmax -273.15,calcounties_proj, snap = 'out', mask = T, touches =F), axes=F,
     col=diverge_hcl(20,'Green-Brown'))

plot(mendocino, border = 'black', add = T)
points(pts, add = T)


mendocino<- calcounties_proj[calcounties_proj$NAME %in% c('Mendocino', 'Kern'),]

mendo_crop <- terra::crop(cal21, mendocino, mask=TRUE, touches = F)

plot(mendo_crop, axes=F, box = F, col = rev(brewer.pal(n = 30, name = "RdBu")))
plot(mendocino, border = 'blue', lwd = 7)


setwd("/Users/monicadonegan/Downloads/tmin_data") # setwd() doesn't carry over between chunks
rasterList<-list.files("/Users/monicadonegan/Downloads/tmin_data") 
AMminTs <- stack(rasterList)
setwd("/Users/monicadonegan/Downloads/tmax_data") # setwd() doesn't carry over between chunks
rasterList2<-list.files("/Users/monicadonegan/Downloads/tmax_data") 
AMmaxTs <- stack(rasterList2)# Make a raster stack from the list of files in the folder
years <- gsub("tmmn_year_gridmet_", "", rasterList)
years <- gsub(".CA_NV.tif", "", years)
names(AMminTs) <- years
names(AMmaxTs) <- years


#35.2387	-118.7615 Je115
#39.050783	-123.13742 D06
xy<- data.frame(x = c(-123.13742, -118.7615), y = c(39.050783, 35.2387))
pts <- matrix(nrow=2, ncol=2)
pts[,1] <- xy[,1]
pts[,2] <- xy[,2]

min_temps <- extract(AMminTs, pts)
max_temps <- extract(AMmaxTs, pts)

all_temps<- data.frame(t(min_temps))
max_temps<- data.frame(t(max_temps))
all_temps<- rbind(all_temps, max_temps)
all_temps<- all_temps - 273.15
all_temps$year<- as.numeric(rep(years,2))
all_temps$Temp_measure<- c(rep('Tmin', 42), rep('Tmax', 42))
colnames(all_temps)[1:2]<- c('Hopland', 'Bakersfield')

library(data.table)
melt_alltemp<- melt(as.data.table(all_temps), id.vars = c('Temp_measure', 'year'), measusre.vars = c('Hopland', 'Bakersfield'))
colnames(melt_alltemp)[4]<- 'degrees Celsius'
colnames(melt_alltemp)[3]<- 'Location'

library(ggplot2)
melt_alltemp %>%
  ggplot(aes(year, `degrees Celsius`, group = Location, color = Location)) + 
  geom_line() + 
  geom_point() + 
  facet_wrap(~Temp_measure, ncol = 1, scales = 'free') + 
  theme_bw()



#### looking at correlation with stunting 
#change name to stunt
stunt<- read.csv('~/Downloads/stunting_april.csv') # april 
stunt_may<- read.csv('~/Downloads/stunting_maydata.csv')


colnames(stunt)<- stunt[1,]
stunt<- stunt[2:478,]
colnames(stunt)[5:16] <- c('a', 'b','c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l')

stunt_melt<- cbind(stunt[1:4], stack(stunt[5:16]))
stunt_melt$ind<- as.character(stunt_melt$ind)
stunt_melt<- stunt_melt %>%
  mutate(ind = case_when(ind %in% c('a', 'b', 'c') ~ 'W',
                         ind %in% c('d', 'e', 'f') ~'WC', 
                         ind %in% c('g', 'h', 'i') ~ 'EC',
                         ind %in% c('j', 'k', 'l') ~ 'E', 
                         TRUE ~ 'problem'))
stunt_melt$values[stunt_melt$values == 'x'] <- -10
#stunt_melt$values[stunt_melt$values == 'n/a'] <- 0
stunt_melt$values[stunt_melt$values == 'N/A'] <- 0
stunt_melt$values[is.na(stunt_melt$values)] <- -10
stunt_melt$values<- as.numeric(stunt_melt$values)
stunt_melt <- stunt_melt[stunt_melt$values != -10,] 


stunt_average<- stunt_melt %>% 
  filter(ind == 'W') %>%
  filter(`flag color` == 'pink') %>%
  group_by(Variety, `flag color`) %>%
  summarize(ave_height = mean(values), sd_height= sd(values), count= n())

inoc_2021<- sum_byvine_2021 %>%
  filter(Treatment != 'control') %>%
  group_by(Variety) %>% 
  summarize(all_pos = sum(pos_any == 1), n = n()) %>%
  mutate(perc_pos = all_pos / n * 100)

merged_stunt<- merge(stunt_average, inoc_2021, by = 'Variety')
library("ggpubr")
ggscatter(merged_stunt, x = "ave_height", y = "perc_pos", 
          add = "reg.line", conf.int = TRUE, 
          label = "Variety",
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Average Height in control plants, April 22", ylab = "Inoculation Success, 2021")
