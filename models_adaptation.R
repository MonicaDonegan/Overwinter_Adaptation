library(glmmTMB)
library(car)
library(lme4)

setwd("/Users/monicadonegan/Downloads/TRADEOFF_MANUSCRIPT/HREC_analysis")
source("input_datasets.R")

## Model 1 - Infection success in year 1
infection_success<- glm(Pos2021~ Treatment + Variety,
                family = 'binomial', data = year1_no_control)

summary(infection_success)
Anova(infection_success, test = "Wald")

## Model 2 - Bacterial population sizes before the winter
prewinter_mod<- lm(log~ Variety  + Treatment , data = all_nocontrols)

summary(prewinter_mod)
Anova(prewinter_mod, test = "Chisq")

## Model 3 - Bacterial populations after winter in 2022

postwinter_mod<- lmer(log~ Variety+ month * Treatment  + (1 | VineID), data = all_nocontrols_2022)

summary(postwinter_mod)
Anova(postwinter_mod, test = "Chisq")

## Model 4 - movement across vine 
movement_mod<- glmmTMB(Response~ Treatment + Variety  + as.factor(Side) *as.factor(year)  + (1 | VineID),
    family = 'binomial', data = movement_full_sum)

summary(movement_mod)
Anova(movement_mod, test = "Chisq")

## model 5 - Overwinter pathogen survival
model_overwinter<- glmmTMB((1-as.numeric(value)) ~  Variety + Treatment +  log + as.factor(year)  + (1 | VineID)  
                , family = 'binomial' ,data = overwinterCT_excludeNever)

summary(model_overwinter)
Anova(model_overwinter)

## model 6 - Symptoms before winter, grouped by month 
prewinter_symptom<- glm(symptomatic~   Variety  + Treatment  ,
                 family = 'binomial', data = symptom_treated_2021)
summary(prewinter_symptom)
Anova(prewinter_symptom, test = 'Wald')

## model 7 - Symptoms all years across quadrants and months 
full_symptom<- glmmTMB(symptomatic~   Variety  + Treatment +  month + year *  Direction + (1|VineID),
                     family = 'binomial', data = symptom_treated)
summary(full_symptom)
Anova(full_symptom, test = 'Chisq')

## model 8 - Stunting after the winter
stunt_model<- lmer(normalized_height~value+ Treatment  + year + month  + Variety+ (1| VineID), data = stunt_filter)
summary(stunt_model)
Anova(stunt_model)

## Correlations - susceptibility 
cor.test(sum_ct$log_ct,sum_ct$perc)
cor.test(sum_ct$symptom_perc,sum_ct$perc)
cor.test(sum_ct$log_ct,sum_ct$symptom_perc)

## Correlations - recovery and susceptibility
cor.test(sum_ct$perc,sum_ct$survival_perc)
cor.test(sum_ct$log_ct,sum_ct$survival_perc)
cor.test(sum_ct$symptom_perc,sum_ct$survival_perc)

