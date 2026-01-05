library(ggplot2)
library(ggpubr) 
library(patchwork)
library(ggpattern)
library(dplyr)
library(tidyverse)

setwd("/Users/monicadonegan/Downloads/TRADEOFF_MANUSCRIPT/HREC_analysis")
source("input_datasets.R")
source("../models_adaptation.R")

BlueRed <- colorRampPalette(c("slateblue2", "tomato1"), bias=1, space="rgb", interpolate="linear") 

# Figure 1A 
chill_sum %>%
  ggplot(aes(Date, value, color = location)) + 
  geom_line(linewidth = 1) + 
  scale_color_manual(values = rev(c(  'slateblue2', 'tomato1'))) +
  facet_grid(~name, scales = 'free') + 
  theme_classic() + 
  ylab('Cumulative Hours')+ theme(
  plot.title = element_text(size=11),
  strip.text = element_text(size = 8), 
  strip.text.y = element_text(angle = 0))
  
# Figure 1B - left
fig1b_left <- summarize_allvar_tr %>%
  filter(Treatment %in% c('Bakersfield', 'Hopland'), year == '2021') %>%
  mutate(strain = case_when(Treatment == 'Bakersfield'~ 'Je115', 
                            Treatment == 'Hopland' ~ 'D06', 
                            TRUE~ 'control')) %>%
  ggplot( aes(x=strain,  y = percent_infected, fill = strain)) +
  geom_col(position = 'dodge') +
  scale_colour_manual(aesthetics = c("colour", "fill"),values = BlueRed(2))+  
  annotate("text", label = "n.s.", x = 1.5, y = 75) +
  ylim(0,100) + 
  theme_bw() +
  labs(x = 'Strain', y ='Infection success', fill = "Positive" ) +
  theme(
    legend.position = 'none',
    axis.title.x = element_blank(),
    plot.title = element_text(size=11),
    strip.text = element_text(size = 8), 
    strip.text.y = element_text(angle = 0)
  ) 

# Figure 1B - center
fig1b_mid<- summary %>%
  mutate(strain = case_when(Treatment == 'Bakersfield'~ 'Je115', 
                            Treatment == 'Hopland' ~ 'D06', 
                            TRUE~ 'control')) %>%
  filter(Status == 'Xf+', Side ==1) %>%
  filter(year ==2021, month == '09') %>%
  filter(Treatment %in% c('Bakersfield', 'Hopland')) %>%
  mutate(month = case_when(month == '06'~'June', 
                           month == '07'~ 'July', 
                           month == '09' ~ 'Sept', 
                           TRUE ~ month), 
         month = factor(month, levels = c('June', 'July', 'Sept'))) %>%
  ggplot( aes(x=strain, y=log, fill = Treatment)) +
  geom_boxplot( outlier.size =0.5) +
  geom_point(aes(fill=Treatment), 
             position=position_dodge(width=0.75) , size =0.5) + 
  xlab('Strain') + 
  ylab('log(CFU/g)') + 
  scale_colour_manual(aesthetics = c("colour", "fill"),values = (BlueRed(2)))+
  annotate("text", label = "n.s.", x = 1.5, y = 8) +
  theme_bw() + 
  theme(
    axis.title.x = element_blank(),
    plot.title = element_text(size=11), 
    legend.position = 'none'
  )  

# Figure 1B - right
fig1b_right<- symp_sum %>% 
  mutate(strain = case_when(Treatment == 'Bakersfield'~ 'Je115', 
                            Treatment == 'Hopland' ~ 'D06', 
                            TRUE~ 'control')) %>% 
  filter(! is.na(Direction), year == 2021, Direction == 'W') %>% 
  ggplot(aes(x = strain, y = symptom_perc * 100, fill = strain)) +
  geom_col(position = 'dodge') + 
  annotate("text", label = "p=0.0004", x = 1.5, y = 75) +
  ylab('Symptom Incidence') + 
  xlab('Strain') + 
  scale_colour_manual(aesthetics = c("colour", "fill"),values = (BlueRed(2))) + 
  ylim(0,100) +
  theme_bw()+
  theme(
    axis.title.x = element_blank(),
    legend.position="none",
    plot.title = element_text(size=11), 
    text=element_text(size=11),strip.text.y = element_text(angle = 0))

## Figure 1C - left
fig1c_left<- infection_by_var %>% 
  mutate(Variety = fct_reorder(Variety, perc)) %>%
  ggplot(aes(Variety, perc)) + 
  geom_col() + 
  theme_bw() + 
  coord_flip()

# Figure 1C - center
fig1c_mid<- summary %>%
  mutate(strain = case_when(Treatment == 'Bakersfield'~ 'Je115', 
                            Treatment == 'Hopland' ~ 'D06', 
                            TRUE~ 'control')) %>%
  mutate(Variety = factor(Variety, levels = rev(c('Tinta Francisca', 'Sagrantino', 'Albarino', 'Tinta Amarella', 'Petit Manseng', 'Tempranillo', 'Greco di Tufo', 'Mencia', 'Periquita','Tannat',  'Falanghina', 'Ciliegiolo', 'Teroldego')))) %>%
  filter(Status == 'Xf+', Side ==1, year == 2021, month == '09') %>%
  filter(Treatment %in% c('Bakersfield', 'Hopland')) %>%
  ggplot( aes(x=Variety, y=log)) +
  geom_boxplot( outlier.size =0.5) +
  geom_point(position=position_dodge(width=0.75) , size =0.5) + 
  ylab('log(CFU/g)') + 
  theme_bw() + 
  theme(
    plot.title = element_text(size=11), 
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) + 
  coord_flip()

# Figure 1C - right
fig1c_right<- symp_sum_var %>% 
  mutate(Variety = factor(Variety, levels = rev(c('Tinta Francisca', 'Sagrantino', 'Albarino', 'Tinta Amarella', 'Petit Manseng', 'Tempranillo', 'Greco di Tufo', 'Mencia', 'Tannat', 'Periquita', 'Falanghina', 'Ciliegiolo', 'Teroldego')))) %>%
  filter(! is.na(Direction), year == '2021') %>%
  ggplot(aes(x = Variety, y = symptom_perc * 100)) +
  geom_col(position = 'dodge') + 
  ylab('Percentage of Xf+ Symptomatic Vines') + 
  theme_bw() + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_flip()

fig1<- '
AABBCC
DDEEFF
DDEEFF'

fig1b_left + fig1b_mid + fig1b_right + fig1c_left + fig1c_mid + fig1c_right + plot_layout(design = fig1)

# Figure 2B 
mv_graph<- movement %>%
  ggplot( aes(x=Side, y=infected, fill = Variety)) +
  geom_col(position = 'dodge') +
  scale_y_continuous(breaks = c( 4,8, 12, 16, 20, 24)) +
  facet_wrap(~year)+
  theme_bw() + 
  theme(
    plot.title = element_text(size=11), 
    text=element_text(size=11),strip.text.y = element_text(angle = 0)) + 
  ylab('Xf Positive Vines') + xlab('Vine Location')

# Figure 2C
symp_graph<- symp_sum %>% 
  filter(! is.na(Direction)) %>%
  filter(Treatment != 'control') %>%
  filter(Treatment != 'control') %>%
  ggplot(aes(x = Direction, y = symptom_perc, fill = Treatment)) +
  geom_col(position = 'dodge') + 
  facet_wrap(~year) + 
  ylab('% Symptomatic') + 
  xlab('Vine Location') + 
  scale_colour_manual(aesthetics = c("colour", "fill"),values = rev(BlueRed(2))) + 
  theme_bw()+
  theme(
    plot.title = element_text(size=11), 
    text=element_text(size=11),strip.text.y = element_text(angle = 0))

mv_graph / symp_graph

# Figure 3A 
strain_ow<- overwinter_sum_tr %>%
  mutate(strain = case_when(
    Treatment == 'Hopland' ~ 'D06', 
    Treatment == 'Bakersfield'~'Je115', 
    TRUE ~ 'control'
  )) %>% 
  ggplot( aes(factor(strain, levels = rev(levels(factor(strain)))), recover_perc, fill = strain)) +
  geom_col() + 
  scale_colour_manual(aesthetics = c("colour", "fill"),values = BlueRed(2))+  
  ylab("Overwinter Recovery")+ 
  xlab("Xf strain") +
  ylim(0,100) + 
  facet_wrap(~year) + theme_bw() + 
  theme(
    plot.title = element_text(size=11), 
    text=element_text(size=11),
    strip.text.y = element_text(angle = 0)
  ) + coord_flip()

# Figure 3B
var_ow<- overwinter_sum_var %>%
  ggplot( aes(Variety, recover_perc, group = Variety, fill= Variety)) +
  geom_col(position = "dodge") + 
  ylim(0,100) + 
  scale_colour_manual(aesthetics = c("colour", "fill"), values = c(rep('gray',4), "#619CFF", rep('gray', 5), "#00BA38", "gray", "#F8766D" ))+  
  theme_bw() + 
  theme(
    plot.title = element_text(size=11), 
    text=element_text(size=11),
    strip.text.y = element_text(angle = 0), 
    legend.position = "none"
  )+
  ylab("Overwinter Recovery Percentage")+ 
  facet_wrap(~year) + 
  coord_flip()

strain_ow / var_ow + plot_layout(heights  = c(1, 3))

# Figure 3C
hrec_chill %>%
  filter(Site == 'Hopland', year %in% c('2021-2022', '2022-2023')) %>%
  ggplot(aes(as.Date(Date), Hours, linetype = year)) + 
  geom_line() + 
  ylim(0, 2300) + 
  theme_bw()+
  theme(
    plot.title = element_text(size=11), 
    text=element_text(size=11),strip.text.y = element_text(angle = 0), 
    legend.position = 'inside',
    legend.position.inside = c(0.75, 0.25)
  ) + 
  xlab('Date') + ylab('Hours < 7.2 C')

# Figure 3D
overwinterCT_excludeNever %>%
  mutate(year = case_when(year == "Recovery_22"~"2022", 
                          year == "Recovery_23"~"2023")) %>%
  filter(! is.na(value)) %>%
  #  filter(Treatment== 'Hopland') %>%
  ggplot(aes(x = as.factor(value), y= as.numeric(log), fill = as.factor(value))) + 
  geom_boxplot(outlier.size =0.5) + 
  geom_point(aes(fill=factor(value)), 
             position=position_dodge(width=0.75) , size =0.5) +
  scale_colour_manual(aesthetics = c("colour", "fill"),values = c("red", "gray" ), labels = c("Xf+", "negative"))+
  labs( x = "Recovery", y = "log(CFU/g) previous fall", fill = "Vine Status") +
  facet_wrap(~year) + 
  theme_bw()+
  theme(
    plot.title = element_text(size=11), 
    text=element_text(size=11),
    strip.text.y = element_text(angle = 0)
  )

# Figure 3E
table_anova<- Anova(model_overwinter)
table_anova$`Chisq`<- round(table_anova$`Chisq`, 3)
table_anova$`Pr(>Chisq)`<- signif(table_anova$`Pr(>Chisq)`, 3)
rownames(table_anova)<- c('Host Variety', 'Pathogen Strain', 'log(CFU/g)', 'Year')
gg_table<- ggtexttable(table_anova, theme = ttheme("light"))
gg_table 

##Figure 4 
summary_filter %>%
  mutate(strain = case_when(Treatment == 'Bakersfield'~ 'Je115', 
                            Treatment == 'Hopland' ~ 'D06', 
                            TRUE~ 'control')) %>%
  filter(Status == 'Xf+', Side ==1) %>%
  filter(year %in% c(2022, 2023)) %>%
  #filter(month %in% c('07', '09')) %>%
  filter(Treatment %in% c('Bakersfield', 'Hopland')) %>%
  mutate(month = case_when(month == '06'~'June', 
                           month == '07'~ 'July', 
                           month == '09' ~ 'Sept', 
                           TRUE ~ month), 
         month = factor(month, levels = c('June', 'July', 'Sept'))) %>%
  ggplot( aes(x=month, y=log, fill = strain)) +
  geom_boxplot( outlier.size =0.5) +
  geom_point(aes(fill=strain), 
             position=position_dodge(width=0.75) , size =0.5) + 
  xlab('Month') + 
  ylab('log(CFU/g)') + 
  facet_wrap(~year, scales = 'free_x', space = 'free_x')+
  scale_colour_manual(aesthetics = c("colour", "fill"),values = (BlueRed(2)))+
  theme_bw() + 
  theme(
    plot.title = element_text(size=11)
  )  

## Figure S3 - CT cutoffs 
summary %>% 
  filter(VineID!= '14 33') %>%
  mutate(result = case_when(
    Treatment == 'control' ~ 'control', 
    Response == 1 ~ 'positive', 
    TRUE ~ 'negative')) %>% 
  ggplot(aes(mean, fill = result)) +   geom_histogram() + facet_wrap(~Response, scales = 'free') + xlab('CT') + theme_bw()

# Figure S4A 
false_negatives %>% 
  mutate(year = as.factor(ifelse(name == 'Pos2021', 2021, 2022))) %>%
  ggplot(aes(Variety, n, fill = year)) + 
  geom_col(position = 'stack') + 
  ggtitle('False negatives') + 
  coord_flip() + 
  ggtitle("(A) False negatives by Variety") +  theme(
    plot.title = element_text(size=11))

# Figure S4B
overwinter_fn %>%
  pivot_longer(Pos2021:Pos2023, values_to = 'pos', names_to = 'year') %>%
  mutate(year = gsub('Pos', '', year)) %>% 
  group_by(year, Treatment) %>%  
  dplyr::summarize(infected = sum(pos == 1, na.rm =T), inferred = sum(pos == 2, na.rm =T), n = n()) %>%
  mutate(percent_infected = infected / n * 100) %>%
  pivot_longer(cols = c('infected', 'inferred')) %>%
  filter(Treatment %in% c('Bakersfield', 'Hopland')) %>%
  mutate(strain = case_when(Treatment == 'Bakersfield'~ 'Je115', 
                            Treatment == 'Hopland' ~ 'D06', 
                            TRUE~ 'control')) %>%
  ggplot( aes(x=strain, y=value)) +
  geom_col_pattern(aes(fill =strain, pattern = name), position = 'stack') +
  scale_colour_manual(aesthetics = c("colour", "fill"),values = BlueRed(2))+  
  scale_pattern_manual(name = 'pattern', values = c( 'none', 'crosshatch'), labels = c('qPCR', 'Inferred'))+
  facet_wrap(~year) + 
  theme_bw() +
  labs(x = 'Year', y ='Count of infected vines', fill = "Strain" ) +
  theme(
    plot.title = element_text(size=11),
    strip.text = element_text(size = 8), 
    strip.text.y = element_text(angle = 0), 
    legend.key.size = unit(1.25, 'cm')
  ) + ggtitle('(B) Inferred Positives')

# Figure S9
corr_1 <- sum_ct %>% ggplot(aes(log_ct, perc)) + 
  geom_point() + geom_smooth(method = lm) + 
  xlab('Log(CFU/g)') + ylab('Infection Success') + theme_bw()

corr_2 <- sum_ct %>% ggplot(aes(perc, symptom_perc)) + 
  geom_point() + geom_smooth(method = lm) + 
  xlab('Infection Success') + ylab('Sympt. Incidence') + theme_bw()

corr_3 <- sum_ct %>% ggplot(aes(symptom_perc, log_ct)) + geom_point() + geom_smooth(method = lm) + 
  xlab('Sympt. Incidence') + ylab('Log(CFU/g)') + theme_bw()

corr_1 + corr_2 + corr_3

# Figure S11 
melted_symp %>%
  mutate(Strain = case_when(
    Treatment == 'Bakersfield' ~ 'Je115', 
    TRUE ~ 'D06'
  )) %>% 
  filter(Direction== 'West', !is.na(month)) %>%
  ggplot(aes(Strain, value, fill = variable)) + 
  geom_col(position = 'dodge') + 
  facet_grid(month~year) + 
  labs(y = "Percentage of Xf+ vines") + 
  ggtitle("Western Side")+ 
  theme_bw()

# Figure S12
stunt_merge %>%
  filter(! is.na(Treatment), Treatment != 'control', ind == 'W',!is.na(value)) %>%
  ggplot(aes(x = as.factor(month), y = normalized_height, fill = as.factor(value))) + 
  geom_boxplot( aes(fill=factor(value)), outlier.size =0.5) +
  geom_point( aes(fill=factor(value)), 
              position=position_dodge(width=0.75) , size =0.5) + 
  facet_grid(~year) +
  scale_colour_manual(aesthetics = c("colour", "fill"),values = c("gray", "red"), labels = c("uninfected", "Xf+"))+
  ylab('Normalized height') + 
  xlab('Month') + 
  theme_bw() + 
  theme(
    plot.title = element_text(size=11), 
    text=element_text(size=11),strip.text.y = element_text(angle = 0))

# Figure S13A
summarize_allvar_tr %>%
  filter(Treatment %in% c('Bakersfield', 'Hopland')) %>%
  mutate(strain = case_when(Treatment == 'Bakersfield'~ 'Je115', 
                            Treatment == 'Hopland' ~ 'D06', 
                            TRUE~ 'control')) %>%
  ggplot( aes(x=year, y=infected, fill = strain)) +
  geom_col(position = 'dodge') +
  scale_colour_manual(aesthetics = c("colour", "fill"),values = BlueRed(2))+  
  theme_bw() +
  labs(x = 'Year', y ='Count of infected vines', fill = "Strain" ) +
  theme(
    plot.title = element_text(size=11),
    strip.text = element_text(size = 8), 
    strip.text.y = element_text(angle = 0)
  ) 

# Figure S13B
movement_tr  %>%
  filter(Treatment %in% c('Bakersfield', 'Hopland')) %>%
  mutate(strain = case_when(Treatment == 'Bakersfield'~ 'Je115', 
                            Treatment == 'Hopland' ~ 'D06', 
                            TRUE~ 'control')) %>%
  ggplot( aes(x=Side, y=infected, fill = strain)) +
  geom_col(position = 'dodge') +
  scale_colour_manual(aesthetics = c("colour", "fill"),values = BlueRed(2))+  
  theme_bw() +
  facet_grid(Variety~year) +
  labs(x = 'Quadrant', y ='Xf Positive vines') +
  theme(
    plot.title = element_text(size=11),
    strip.text = element_text(size = 8)
  ) 

# Figure S13C
summarize_tr_var  %>%
  filter(Treatment %in% c('Bakersfield', 'Hopland')) %>%
  mutate(strain = case_when(Treatment == 'Bakersfield'~ 'Je115', 
                            Treatment == 'Hopland' ~ 'D06', 
                            TRUE~ 'control')) %>%
  ggplot( aes(x=year, y=infected, fill = strain)) +
  geom_col(position = 'dodge') +
  scale_colour_manual(aesthetics = c("colour", "fill"),values = BlueRed(2))+  
  theme_bw() +
  facet_wrap(~Variety) +
  labs(x = 'Year', y ='Count of infected vines', fill = "Strain" ) +
  theme(
    plot.title = element_text(size=11),
    strip.text = element_text(size = 8), 
    strip.text.y = element_text(angle = 0)
  ) 

## Figure S14
sus_1<- sum_ct %>% ggplot(aes(perc, recover_perc)) + 
  geom_point() + geom_smooth(method = lm) + 
  xlab('Infection Success') + ylab('Overwinter Recovery') + theme_bw()

sus_2<- sum_ct %>% ggplot(aes(log_ct, recover_perc)) + 
  geom_point() + geom_smooth(method = lm) + 
  xlab('log(CFU/g)') + ylab('Overwinter Recovery') + theme_bw()

sus_3<- sum_ct %>% ggplot(aes(symptom_perc, recover_perc)) + 
  geom_point() + geom_smooth(method = lm) + 
  xlab('Sympt. Incidence') + ylab('OverwinterRecovery') + theme_bw()

sus_1 + sus_2 + sus_3

