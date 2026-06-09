library(dplyr)
library(stringr)
library(tidyverse)

input_qpcr_files<- function(setup_path, raw_data_path){
  qpcr_setup<- read.csv(setup_path)
  qpcr_setup$Run.Date<- as.Date(qpcr_setup$Run.Date, format="%m/%d/%y")
  setwd(raw_data_path)
  list<- list.files()
  outputs<- c()
  for (i in 1:length(list)) {
    name <- list[i]
    name<- str_remove(name, '.csv')
    plate<- read.csv(list[i])
    datestring<- as.POSIXct(str_split(name, "_")[[1]][3],format="%m-%d-%y")    
    if (plate[28,1] == 'Well') {
      colnames(plate)<- plate[28,]
      plate<- plate[29:124,]
      
    } else if (plate[29,1] == 'Well') {
      colnames(plate)<- plate[29,]
      plate<- plate[30:125,]
    }
    colnames(plate)[2]<- 'Sample.Name'
    plate$Run.Date<- as.Date(datestring)
    plate$Ct[plate$Ct == 'Undetermined'] <- 40
    plate$Ct <- as.numeric(plate$Ct)
    name<- str_split(name, "_")[[1]][2]
    assign(name, plate)
    outputs<- c(outputs, plate)
  }
  rm(plate)
  local_env<- environment()
  Pattern1<-grep("plate",names(local_env),value=TRUE)
  Pattern1_list<-do.call("list",mget(Pattern1))
  all_data<- bind_rows(Pattern1_list, .id = "column_label")
  all_data2<- merge(all_data, qpcr_setup, by=c('Run.Date', 'Sample.Name'))
  return(all_data2)
}

determine_positives<- function(summary, ct_cutoff){ 
  summary$CFUperG <- (summary$mean * -0.245 + 12.8)
  ### add an infected column 
  summary$Status<- ''
  summary$Response <- 0
  for (i in 1:dim(summary)[1]) {
    if (summary$mean[i] < ct_cutoff ) {
      summary$Status[i]<- 'Xf+'
      summary$Response[i]<- 1
    } else {
      summary$Status[i]<- '(-)'
    }
  }
  
  summary$log <- (summary$CFUperG)
  summary$log[summary$Status != 'Xf+'] <- 0 
  
  summary$year<- format(as.Date(summary$Sample.Date, format="%m/%d/%y"),"%Y")
  summary$month<- format(as.Date(summary$Sample.Date, format="%m/%d/%y"),"%m")
  return(summary)
}


check_culture<- function(sum, culture, summary){ 
  culture <- culture %>% filter(year == sum$year[1])
  for(i in 1:dim(culture)[1]) {
    index<- which(sum$VineID == culture$VineID[i])
    if(identical(index, integer(0))) {
      newrow = data.frame(Variety = summary$Variety[summary$VineID==culture$VineID[i]][1], Treatment = summary$Treatment[summary$VineID==culture$VineID[i]][1], VineID = culture$VineID[i], Row = culture$Row[i], Vine = culture$Vine[i], year = sum$year[1], 
                          pos_tot=0, n=0, pos_any=1)
      sum <- rbind(newrow, sum)
    } else {
      sum$pos_any[index]<- 1
    }
  }
  return(sum)
}

sum_by_vine_general<- function(df, year_filter ) {
  new_sum<- df %>%
    filter(year == year_filter, 
           Treatment != 'Control', 
           !is.na(Treatment)) %>%
    group_by(Variety, Treatment, VineID, Row, Vine, year) %>%
    summarise(pos_tot = sum(Response), n = n()) %>%
    mutate(pos_any = ifelse(pos_tot >0 , 1, 0))
  return(new_sum)
}
