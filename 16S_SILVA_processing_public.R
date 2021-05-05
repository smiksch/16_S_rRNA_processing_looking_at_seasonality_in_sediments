#load needed libraries#####
library(tidyverse)
library(reshape2)
library(vegan)
library(ggdendro)
library(ggpubr)
library(zoo)
library(dendextend)
library(ggplotify)
library(cowplot)
#set and check my wd#####
getwd()
setwd("PATH/TO/DATA/AND/WORKSPACE")
getwd()
#load my self created functions#####
#return tabel with relativ abundances
get.pct <- function(COLUMN) {
  (COLUMN/sum(COLUMN))*100
}

#return "others" value
get.others <- function(DATAFRAM) {
  x <- dplyr::select_if(DATAFRAM, is.numeric)
  y <- rbind(x, 100 - colSums(x))
  withothers <- DATAFRAM%>% 
    add_row(taxonomicpath= 'ZZOther')
  z <- cbind(y, taxonomicpath = withothers$taxonomicpath)
  return(z)
}

#sums up data of small Taxa on phylum, or in case of Proteobacteria on class, level (Bacteria)
sum.up.everything.with.low.abundance <-
  function(
    DATA, #your datafram
    MINIMUM #cut of in % e.g. 1 for 1% abundance
  ){
    over_min <- 
      DATA %>%
      filter_at(.vars = vars(which(sapply(.,class)=="numeric")), .vars_predicate = any_vars(. > MINIMUM))
    
    under_min <- 
      anti_join(DATA, over_min, by = "taxonomicpath")%>% #gets all rows that are not in over_1
      separate(taxonomicpath, 
               into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Specis", "Strain"),
               sep = ";")%>% #to analyse all taxalevels
      select(which(sapply(.,class)=="numeric"), Phylum, Class)%>% #only use the levels I want to sum up
      mutate(Phylum = ifelse(Phylum == "Proteobacteria", Class, Phylum))%>% #use the Class name for Proteobacteria
      select(which(sapply(.,class)=="numeric"), Phylum)%>%
      group_by(Phylum)%>% #Grouping needed for the summary 
      summarise_at(.vars = vars(which(sapply(.,class)=="numeric")), .funs = sum)%>% #finally creates subtotals
      filter_at(.vars = vars(which(sapply(.,class)=="numeric")), .vars_predicate = any_vars(. > MINIMUM))
    under_min$Phylum <- 
      sub("$", "_zother", under_min$Phylum)#renames vars
    under_min$Phylum <- 
      sub("^", "Bacteria_", under_min$Phylum)#renames vars
    under_min$Phylum <- 
      sub("Alphaproteobacteria", "Proteobacteria_Alphaproteobacteria", under_min$Phylum)#renames vars
    under_min$Phylum <- 
      sub("Deltaproteobacteria", "Proteobacteria_Deltaproteobacteria", under_min$Phylum)#renames vars
    under_min$Phylum <- 
      sub("Gammaproteobacteria", "Proteobacteria_Gammaproteobacteria", under_min$Phylum)#renames vars
    
    OUTPUT  <- 
      under_min %>%
      separate(Phylum, 
               into = c("Domain", "Phylum", "Class", "Family"),
               sep = "_")%>% 
      tidyr::unite(taxonomicpath, Domain, Phylum, Class, Family, sep = ";")%>%
      full_join(over_min)%>%
      arrange(taxonomicpath)%>% 
      get.others()%>%
      mutate(
        taxonomicpath = str_replace(
          taxonomicpath, "No Relative;NA;NA;NA;NA;NA;NA;NA", "ZNoRelative"
        )
      )
    
    return(OUTPUT)
  }

#calculates mean and SD from ASV data based on seasonal groups
get.mean.SD.of.group<- 
  function(DF=.){
    SEASON<-
      DF$season%>%
      as.vector()%>%
      unique()
    SUMMARY<-
      DF%>%
      select(
        tax,
        ASV,
        value
      )%>%
      drop_na()%>%
      group_by(ASV)%>%
      dplyr::summarise(
        mean=mean(value),
        sd=sd(value),
        tax=unique(tax)
      )%>%
      mutate(
        season=SEASON
      )
    return(SUMMARY)
  }

#function to make ASV plots for sediment and water from lists of groups 
get.asw.plots3<-
  function(DF, TOP=100){
    TOPDF<-
      DF%>%
      filter(str_detect(season,"spring"))%>%
      top_n(TOP,mean)%>%
      mutate(
        ordering=mean+max
      )
    
    ORDER<-
      TOPDF%>%
      select(
        ASV,
        ordering
      )
    
    FORPLOT<-
      DF%>%
      filter(
        ASV %in% TOPDF$ASV
      )%>%
      full_join(ORDER)%>%
      mutate(
        min = if_else(min < 0, 0, min)
      )
    
    FIRST<-"Top"
    SECOND<-"ASVs in sediments"
    PLOTNAME<-paste(FIRST,TOP,SECOND)
    
    PLOT<-
      FORPLOT%>%
      mutate(
        ASV=fct_reorder(
          ASV, ordering, .desc = T)
      )%>%
      ggplot()+
      geom_point(
        aes(
          x=ASV,
          y=mean,
          color=season
        ),size=2
      )+
      geom_errorbar(
        aes(
          x=ASV,
          ymin=min,
          ymax=max,
          color=season
        )
      )+
      scale_color_manual(values = c(
        "winter"="black",
        "spring"="green"
      )
      )+
      scale_y_continuous(
        expand = c(0, 0),
        limits = c(0, NA)
      )+
      theme_minimal()+
      theme(
        axis.text.x = element_text(
          angle = 90,
          size = 6,
          hjust=0.95,
          vjust=0.2
        ),
        strip.placement = "inside",
        panel.grid = element_blank(),
        panel.background = element_rect(fill=NA,color = "black",size = 1),
        legend.position="bottom",
        legend.text = element_text(size = 8),
        axis.ticks.y = element_line(color = "black", size = 1.5),
        axis.ticks.x = element_line(color = "black", size = 1.5),
        aspect.ratio = 9/21
      )+
      xlab("ASVs")+
      ylab("relativ ASV abundance[%]")+
      ggtitle(PLOTNAME)
    
    return(PLOT)
  }

get.asw.plots4<-
  function(DF, TOP=100){
    TOPDF<-
      DF%>%
      filter(str_detect(season,"spring"))%>%
      top_n(TOP,mean)%>%
      mutate(
        ordering=mean+max
      )
    
    ORDER<-
      TOPDF%>%
      select(
        ASV,
        ordering
      )
    
    FORPLOT<-
      DF%>%
      filter(
        ASV %in% TOPDF$ASV
      )%>%
      full_join(ORDER)%>%
      mutate(
        min = if_else(min < 0, 0, min)
      )
    
    FIRST<-"Top"
    SECOND<-"ASVs in water"
    PLOTNAME<-paste(FIRST,TOP,SECOND)
    
    PLOT<-
      FORPLOT%>%
      mutate(
        ASV=fct_reorder(
          ASV, ordering, .desc = T)
      )%>%
      ggplot()+
      geom_point(
        aes(
          x=ASV,
          y=mean,
          color=season,
          shape=season
        ),size=2
      )+
      geom_errorbar(
        aes(
          x=ASV,
          ymin=min,
          ymax=max,
          color=season
        )
      )+
      scale_color_manual(values = c(
        "winter"="black",
        "spring2018"="green",
        "spring2019"="green4"
      )
      )+
      scale_y_continuous(
        expand = c(0, 0),
        limits = c(0, NA)
      )+
      theme_minimal()+
      theme(
        axis.text.x = element_text(
          angle = 90,
          size = 6,
          hjust=0.95,
          vjust=0.2
        ),
        strip.placement = "inside",
        panel.grid = element_blank(),
        panel.background = element_rect(fill=NA,color = "black",size = 1),
        legend.position="bottom",
        legend.text = element_text(size = 8),
        axis.ticks.y = element_line(color = "black", size = 1.5),
        axis.ticks.x = element_line(color = "black", size = 1.5),
        aspect.ratio = 9/21
      )+
      xlab("ASVs")+
      ylab("relativ ASV abundance[%]")+
      ggtitle(PLOTNAME)
    
    return(PLOT)
  }

#to plot specific taxa separately
get.asw.plots.onlisttaxa<-
  function(DF, TOP=100){
    TAXON<-
      unique(DF$group)
    TOPDF<-
      DF%>%
      filter(str_detect(season,"spring"))%>%
      top_n(TOP,mean)%>%
      mutate(
        ordering=mean+max
      )
    
    ORDER<-
      TOPDF%>%
      select(
        ASV,
        ordering
      )
    
    FORPLOT<-
      DF%>%
      filter(
        ASV %in% TOPDF$ASV
      )%>%
      full_join(ORDER)%>%
      mutate(
        min = if_else(min < 0, 0, min)
      )
    
    FIRST<-"Top"
    SECOND<-"ASVs in sediments"
    PLOTNAME<-paste(TAXON,FIRST,TOP,SECOND)
    
    PLOT<-
      FORPLOT%>%
      mutate(
        ASV=fct_reorder(
          ASV, ordering, .desc = T)
      )%>%
      ggplot()+
      geom_point(
        aes(
          x=ASV,
          y=mean,
          color=season
        ),size=1.5
      )+
      geom_errorbar(
        aes(
          x=ASV,
          ymin=min,
          ymax=max,
          color=season
        )
      )+
      scale_color_manual(values = c(
        "winter"="black",
        "spring"="green"
      )
      )+
      scale_y_continuous(
        expand = c(0, 0),
        limits = c(0, NA)
      )+
      theme_minimal()+
      theme(
        axis.text.x = element_text(
          angle = 90,
          size = 6,
          hjust=0.95,
          vjust=0.2
        ),
        strip.placement = "inside",
        panel.grid = element_blank(),
        panel.background = element_rect(fill=NA,color = "black",size = 1),
        legend.position="bottom",
        legend.text = element_text(size = 8),
        axis.ticks.y = element_line(color = "black", size = 1.5),
        axis.ticks.x = element_line(color = "black", size = 1.5),
        aspect.ratio = 9/21
      )+
      xlab("ASVs")+
      ylab("relativ ASV abundance[%]")+
      ggtitle(PLOTNAME)
    
    return(PLOT)
  }

#data import and normalization#####
#for own use imput comunity data in read_tsv() function or addaped further depending on input
SILVAngs_out <- read_tsv('seasonpaper_v138---ssu---fingerprint----Total---sim_93---tax_silva---td_20.csv',
                         col_names = TRUE, col_types = NULL
)%>%
  separate(`taxonomic path`, 
           into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Specis", "Strain"),
           sep = ";") %>%
  filter(str_detect(Domain, "Bacteria") |
           str_detect(Domain, "No Relative")) %>%
  unite(taxonomicpath,
        Domain, Phylum, Class, Order, Family, Genus, Specis, Strain,
        sep = ";"
  )%>%
  select(
    which(sapply(.,class)=="numeric"),taxonomicpath
  )%>%#following part removes absolut singeltons
  mutate(sumVars = rowSums(.[1:ncol(.)-1])) %>%
  filter(sumVars >1) %>%
  select(which(sapply(.,class)=="numeric"),
         taxonomicpath)

#normalise as relative abundance
pct<- 
  SILVAngs_out %>%
  filter(
    !str_detect(taxonomicpath,"Chloroplast")
  )%>%
  mutate_at(.vars = vars(which(sapply(.,class)=="numeric")),.funs = get.pct)%>%
  select(
    !contains("_mock")&
      !sumVars
  )

#write all values in a table ousite R 
write_excel_csv(pct,"pct.csv")

#for PANGEA upload
pct_sval_export<-
  pct%>%
  select(
    contains("s_")|
      contains("w_")|
      taxonomicpath
  )%>%
  melt()%>%
  mutate(
    variable=str_remove(
      variable,
      "sub20000\\."
    ),
    taxonomicpath=str_remove_all(
      taxonomicpath,
      ";NA"
    ),
    taxonomicpath=str_remove(
      taxonomicpath,
      ";$"
    )
  )%>%
  separate(
    variable,
    into = c(
      "type",
      "station",
      "date"
    )
  )%>%
  mutate(
    adding=case_when(
      type=="s" ~ "S",
      type=="w" ~ "W"
    ),
    station=str_remove(
      station,
      "st"
    )
  )%>%
  unite(
    st_type,
    station,
    adding
  )%>%
  mutate(
    st_type=str_remove(
      st_type,
      "_"
    ),
    type="ISF"
  )%>%
  unite(
    tmp,
    type,
    st_type,
    sep = "/"
  )%>%
  unite(
    ID,
    tmp,
    date
  )%>%
  dcast(
    taxonomicpath~...,
    value.var = "value"
  )%>%
  filter_at(
    .vars = vars(which(sapply(.,class)=="numeric")),
    .vars_predicate = any_vars(. > 0)
    )

pct_Hel_export<-
  pct%>%
  select(
    contains("Hel_")|
      taxonomicpath
  )%>%
  melt()%>%
  mutate(
    variable=str_remove(
      variable,
      "sub20000\\.Nr[:digit:]{2}_"
    ),
    taxonomicpath=str_remove_all(
      taxonomicpath,
      ";NA"
    ),
    taxonomicpath=str_remove(
      taxonomicpath,
      ";$"
    )
  )%>%
  dcast(
    taxonomicpath~...,
    value.var = "value"
  )%>%
  filter_at(
    .vars = vars(which(sapply(.,class)=="numeric")),
    .vars_predicate = any_vars(. > 0)
  )
  
write_excel_csv(pct_Hel_export,"Hel_pangaea_16S_tabel.csv")
write_excel_csv(pct_sval_export,"Sval_pangaea_16S_tabel.csv")
#data processing for Plotting######
#filter out low abundant taxa to reduce plot size 
#also here we manually summed up several closely related taxa
over2<-
  pct%>%
  sum.up.everything.with.low.abundance(.,2)%>%
  mutate(
    taxonomicpath=str_remove(taxonomicpath,"Clade II;"),
    taxonomicpath=str_remove(taxonomicpath,"Clade I;"),
    taxonomicpath=str_remove(taxonomicpath,"Clade Ia;"),
    taxonomicpath=str_remove(taxonomicpath,"Clade Ib;"),
    taxonomicpath=str_replace(taxonomicpath,"SAR11 clade;uncultured","SAR11 clade"),
    taxonomicpath=str_replace(
      taxonomicpath,
      "Gammaproteobacteria Incertae Sedis;Unknown Family;uncultured",
      "zother"),
    taxonomicpath=str_replace(taxonomicpath,"Bacteria;Cyanobacteria;zother;NA","ZZOther"),
    taxonomicpath=str_replace(taxonomicpath,
                              "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Colwelliaceae;Colwellia",
                              "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Colwelliaceae.unc"),
    taxonomicpath=str_replace(taxonomicpath,
                              "Bacteria;Proteobacteria;Alphaproteobacteria;uncultured;;NA;NA;NA",
                              "Bacteria;Proteobacteria;Alphaproteobacteria;zother"),
    taxonomicpath=str_replace(taxonomicpath,
                              "Bacteria;Proteobacteria;Gammaproteobacteria;uncultured;;NA;NA;NA",
                              "Bacteria;Proteobacteria;Gammaproteobacteria;zother"),
    taxonomicpath=str_replace(taxonomicpath,";uncultured",".unc"),
    taxonomicpath=str_remove_all(taxonomicpath,";NA"),
    taxonomicpath=str_remove(taxonomicpath,";$")
    )%>%
  group_by(taxonomicpath)%>% #Grouping needed for the summary 
  summarise_at(.vars = vars(which(sapply(.,class)=="numeric")), .funs = sum)

#export table to read exact values
write_excel_csv(over2,"over2.csv")

#now we edit and rearrange some parts of the table to plot it nicely with ggplot
for_bubble_plot<-
  over2%>%
  melt%>%
  rename(
    Sample=variable,
    Taxa=taxonomicpath
  )%>%
  mutate(
    Sample=str_remove(Sample,"sub20000\\."),
    Sample=str_remove(Sample,"combind\\."),
    Sample=str_remove(Sample,"Nr[:digit:]{2}_"),
    Sample=str_remove(Sample,"\\.merged"),
    Sample=str_remove(Sample,"_Svalbard")
  )%>%
  mutate(
    Taxa=str_remove_all(Taxa,";NA"),
    Taxa=str_remove(Taxa,"Bacteria;"),
    Taxa=str_remove(Taxa,";$"),
    Taxa=str_replace_all(Taxa,"[:blank:]","\\."),
    Taxa=str_replace_all(Taxa,"\\(","\\."),
    Taxa=str_replace_all(Taxa,"\\)","\\."),
    Taxa=str_replace_all(Taxa,"_","\\."),
    Taxa=str_replace_all(Taxa,"-","\\."),
    Taxa=str_replace(Taxa,";uncultured",".unc"),
    Phy=str_extract(Taxa,"^[:alpha:]+;"),
    Gen=str_extract(Taxa,";[[:alnum:]\\.]+$"),
    Phy=str_remove(Phy,";$"),
    Gen=str_remove(Gen,"^;"),
    Phy=case_when(
      str_detect(Taxa,";Gammaproteobacteria;") ~ "Gammaproteobacteria",
      str_detect(Taxa,";Deltaproteobacteria;") ~ "Deltaproteobacteria",
      str_detect(Taxa,";Alphaproteobacteria;") ~ "Alphaproteobacteria",
      str_detect(Taxa,"NoRelative") ~ "No Relative",
      str_detect(Taxa,"ZZOther") ~ "Other",
      T ~ Phy
    ),
    station=str_extract(
      Sample,
      "_st[:digit:]{2}"
    ),
    Sample=str_remove(
      Sample,
      "_st[:digit:]{2}"
    ),
    season=case_when(
      str_detect(Sample,"_201[78]12") ~ "winter",
      str_detect(Sample,"_201[89]01") ~ "winter",
      str_detect(Sample,"_201[89]0[7-9]") ~ "summer/fall",
      str_detect(Sample,"_201[89]1[0-1]") ~ "summer/fall",
      str_detect(Sample,"_20180228") ~ "twilight",
      T ~ "spring",
    )
  )%>%
  unite(Sample_tmp,Sample,station,sep = "%")%>%
  mutate(
    Sample_tmp=str_remove(
      Sample_tmp,
      "%"
      ),
    Sample_tmp=str_remove(
      Sample_tmp,
      "NA"
    )
  )%>%
  rename(
    Sample=Sample_tmp
  )%>%
  mutate(
    Sample=str_replace(
      Sample,
      "s_20190913_st07",
      "s_20190913_Ast07"
      )
  )%>%
  unite(label_y,Phy,Gen,remove = F,sep = " ")%>%
  mutate(
    label_y=case_when(
      label_y=="NA NA" ~ Taxa,
      T ~ label_y
    ),
    label_y=str_remove(label_y," NA$")
  )%>%
  unite(label_x,Sample,remove = F,sep = " ")%>%
  arrange(
    Taxa
  )

#now we finally create the bubble plot 
p_bubble<-
  ggplot(for_bubble_plot)+
  geom_point(aes(
    x=Sample,
    y=Taxa,
    size=value,
    color=season
  )
  )+
  scale_size_continuous(range = c(-1, 10),breaks = c(0,2,5,10,15,20,40,60))+
  scale_y_discrete(
    labels=unique(for_bubble_plot$label_y)
  )+
  scale_color_manual(
    values = c(
      "winter"="#1A0E34", 
      "twilight"="#929489", 
      "spring"="#A7CE9D", 
      "summer/fall"="#6E3E67"
    )
  )+
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 90,size = 7,vjust = 0.3),
    axis.text.y = element_text(size = 7, vjust = 0.2),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    strip.text = element_text(size = 8),
    strip.placement = "inside",
    strip.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.border = element_rect(fill=NA,color = "black"),
    legend.position="bottom",
    legend.text = element_text(size = 8),
    axis.ticks.y = element_line(color = "black", size = 0.5),
    axis.ticks.x = element_line(color = "black", size = 0.5),
    aspect.ratio = 0.5/1
  )

p_bubble

pdf("test.pdf",paper = "a4r",width = 16.5,height = 15)
p_bubble
dev.off()

#here we determine the frequency of Chloroplasts in relation to all reads 
pct_chlor<- 
  SILVAngs_out %>%
  mutate_at(.vars = vars(which(sapply(.,class)=="numeric")),.funs = get.pct)%>%
  select(
    !contains("_mock")&
      !sumVars
  )%>%
  filter(
    str_detect(taxonomicpath,"Chloroplast")
  )

#creation of data table used to create chloroplast bar plot 
for_chloro_plot<-
  pct_chlor%>%
  melt%>%
  rename(
    Sample=variable,
    Taxa=taxonomicpath
  )%>%
  mutate(
    Sample=str_remove(Sample,"sub20000\\."),
    Sample=str_remove(Sample,"combind\\."),
    Sample=str_remove(Sample,"Nr[:digit:]{2}_"),
    Sample=str_remove(Sample,"\\.merged"),
    Sample=str_remove(Sample,"_Svalbard")
  )%>%
  mutate(
    Taxa=str_remove_all(Taxa,";NA"),
    Taxa=str_remove(Taxa,"Bacteria;"),
    Taxa=str_remove(Taxa,";$"),
    Taxa=str_replace_all(Taxa,"[:blank:]","\\."),
    Taxa=str_replace_all(Taxa,"\\(","\\."),
    Taxa=str_replace_all(Taxa,"\\)","\\."),
    Taxa=str_replace_all(Taxa,"_","\\."),
    Taxa=str_replace_all(Taxa,"-","\\."),
    Taxa=str_replace(Taxa,";uncultured",".unc"),
    station=str_extract(
      Sample,
      "_st[:digit:]{2}"
    ),
    Sample=str_remove(
      Sample,
      "_st[:digit:]{2}"
    ),
    season=case_when(
      str_detect(Sample,"_201[78]12") ~ "winter",
      str_detect(Sample,"_201[89]01") ~ "winter",
      str_detect(Sample,"_201[89]0[7-9]") ~ "summer/fall",
      str_detect(Sample,"_201[89]1[0-1]") ~ "summer/fall",
      str_detect(Sample,"_20180228") ~ "twilight",
      T ~ "spring",
    )
  )%>%
  unite(Sample_tmp,Sample,station,sep = "%")%>%
  mutate(
    Sample_tmp=str_remove(
      Sample_tmp,
      "%"
    ),
    Sample_tmp=str_remove(
      Sample_tmp,
      "NA"
    )
  )%>%
  rename(
    Sample=Sample_tmp
  )%>%
  mutate(
    Sample=str_replace(
      Sample,
      "s_20190913_st07",
      "s_20190913_Ast07"
    )
  )%>%
  unite(label_x,Sample,remove = F,sep = " ")%>%
  arrange(
    Taxa
  )

#create chloroplast bar plot 
p_bar_chloro<-
  ggplot(for_chloro_plot)+
  geom_bar(aes(
    x=Sample,
    y=value
  ),
  fill="gray69",
  stat="identity"
  )+
  coord_trans()+
  theme_minimal()+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    strip.text = element_text(size = 8),
    strip.placement = "inside",
    strip.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.border = element_rect(fill=NA,color = "black"),
    legend.position="bottom",
    legend.text = element_text(size = 8),
    axis.ticks.y = element_line(color = "black", size = 0.5),
    axis.ticks.x = element_line(color = "black", size = 0.5),
    aspect.ratio = 1/19.96
  )

#merge plots 
P_arranged<-
  ggarrange(
    p_bar_chloro,
    p_bubble,
    nrow = 2,
    align = "v",
    heights = c(0.1,1)
    )

P_arranged

#final output for manual beauty tweaking in AI 
pdf(
  "bubble_p_Bac2_Sval_all_forRef.pdf",
  paper = "a4r",width = 16.5,height = 15
  )
P_arranged
dev.off()
#NMDS plot####
#import SILVA out with Probandt Data
SILVA_out2 <- read_tsv('comparison_v132---ssu---fingerprint----Total---sim_93---tax_silva---td_20.csv',
                      col_names = TRUE, col_types = NULL
)%>%
  separate(`taxonomic path`, 
           into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Specis", "Strain"),
           sep = ";") %>%
  filter(str_detect(Domain, "Bacteria") |
           str_detect(Domain, "No Relative")) %>%
  unite(taxonomicpath,
        Domain, Phylum, Class, Order, Family, Genus, Specis, Strain,
        sep = ";"
  )%>%
  select(
    which(sapply(.,class)=="numeric"),
    taxonomicpath
  )%>%#following part removes absolut singeltons
  mutate(sumVars = rowSums(.[1:ncol(.)-1])) %>%
  filter(sumVars >1) %>%
  select(!sumVars)

pct2<-
  SILVA_out2%>%
  filter(
    !str_detect(taxonomicpath,"Chloroplast")
  )%>%
  mutate_at(.vars = vars(which(sapply(.,class)=="numeric")),.funs = get.pct)%>%
  select(
    which(sapply(.,class)=="numeric"),
    taxonomicpath,
  )

#run clustering
cluster_all<-
  pct2%>%
  select(
    which(sapply(.,class)=="numeric")
  )%>%
  data.matrix()%>%
  t()%>%
  metaMDS(distance = "bray")

plot(cluster_all)

#extract score for ggplot 
scores_all<-
  as.data.frame(scores(cluster_all))

#make plotting nicer 
for_nmds_plotting_all <-
  scores_all%>%
  mutate(
    ID=row.names(scores_all),
    Ano_ID=row.names(scores_all),
    ID=str_remove(ID,".10000reads"),
    ID=str_replace_all(ID,"C_D","_CCP-D"),
    ID=str_replace_all(ID,"C_G","_CCP-G"),
    ID=str_replace_all(ID,"_00_10","0-10cm_201409"),
    ID=str_replace_all(ID,"_10_20","10-20cm_201409")
  )%>%
  separate(
    ID, c("1", "2", "3"), sep = "_"
  )%>%
  select(
    !"1"
  )%>%
  rename(
    station="2",
    date="3"
  )%>%
  mutate(
    location=case_when(
      station=="Hel" ~ "Hel",
      str_detect(station,"CCP-[D,G]") ~ "Hel",
      TRUE ~ "Sval"
    )
  )%>%
  mutate(
    season=case_when(
      str_detect(date,"201[78]12") ~ "winter",
      str_detect(date,"201[89]01") ~ "winter",
      str_detect(date,"201[89]0[7-9]") ~ "summer/fall",
      str_detect(date,"201[89]1[0-1]") ~ "summer/fall",
      str_detect(date,"201409") ~ "summer/fall",
      str_detect(date,"20180228") ~ "twilight",
      T ~ "spring",
    ),
    ID=str_remove(station, "01"),
    ID=str_remove(ID, "02"),
    ID=str_remove(ID, "_3"),
    ID=str_remove(ID, "_5")
  )

#run ANOSIM Stats
matrix_for_ano_all<-
  pct2%>%
  select(
    which(sapply(.,class)=="numeric")
  )%>%
  data.matrix()%>%
  t()

ano_location_all <- 
  anosim(
    matrix_for_ano_all, for_nmds_plotting_all$location, distance = "bray", permutations = 9999
  )

ano_location_all

ano_season_all <- 
  anosim(
    matrix_for_ano_all, for_nmds_plotting_all$season, distance = "bray", permutations = 9999
  )

ano_season_all

#create polygones
for_circal<-
  for_nmds_plotting_all%>%
  select(NMDS1,
         NMDS2,
         location
  )

circal_Hel <- for_circal[for_circal$location == "Hel", ][chull(for_circal[for_circal$location == 
                                                                            "Hel", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
circal_Sval <- for_circal[for_circal$location == "Sval", ][chull(for_circal[for_circal$location == 
                                                                              "Sval", c("NMDS1", "NMDS2")]), ]  # hull values for grp B

circal.data <- rbind(circal_Hel, circal_Sval)%>%
  mutate(
    location=str_replace_all(location,"Hel","Helgoland"),
    location=str_replace_all(location,"Sval","Svalbard")
  )#combine grp.a and grp.b
circal.data

#make the plot 
nmds_p<-
ggplot(for_nmds_plotting_all)+
  geom_polygon(data=circal.data,
               mapping = aes(x=NMDS1,y=NMDS2,linetype = location,group=location),
               fill=NA,
               color="black",
               size=1
  ) +
  geom_point(aes(x = NMDS1,
                 y = NMDS2,
                 color=season,
                 shape=ID),
             size=4.5
  )+
  scale_color_manual(
    values = c(
      "winter"="#1A0E34", 
      "twilight"="#929489", 
      "spring"="#A7CE9D", 
      "summer/fall"="#6E3E67"
    )
  )+
  scale_shape_manual(values = c("CCP-D0-10cm"=15,
                                "CCP-D10-20cm"=0,
                                "CCP-G0-10cm"=16,
                                "CCP-G10-20cm"=1,
                                "Hel"=7,
                                "st05"=17,
                                "st06"=18,
                                "st07"=25,
                                "st23"=8
  )
  )+
  geom_text(aes(x=-0.73,y=-0.5,label="stress:"))+
  geom_text(aes(x=-0.605,y=-0.55,label=cluster_all$stress))+
  geom_text(aes(x=-0.7,y=-0.6,label="location R:"))+
  geom_text(aes(x=-0.63,y=-0.65,label=ano_location_all$statistic))+
  geom_text(aes(x=-0.61,y=-0.7,label="location Significance:"))+
  geom_text(aes(x=-0.75,y=-0.75,label=ano_location_all$signif))+
  geom_text(aes(x=-0.7,y=-0.8,label="season R:"))+
  geom_text(aes(x=-0.605,y=-0.85,label=ano_season_all$statistic))+
  geom_text(aes(x=-0.61,y=-0.9,label="season Significance:"))+
  geom_text(aes(x=-0.605,y=-0.95,label=ano_season_all$signif))+
  theme_minimal_grid()

nmds_p

#add to file 
pdf("nmds_plot.pdf",paper = "a4r",width = 16.5,height = 15)
nmds_p
dev.off()

#plot to display changes and cluster in a tree based on dissimilarity#####
#select values to compare all others too in this case the Winter in sediments
my_basline<-
  for_bubble_plot%>%
  select(
    label_y,
    value,
    label_x,
    season
  )%>%
  filter(
    str_detect(label_x,"s_")&
    str_detect(season,"winter")
  )

#now lets make the means and SD for all taxa in the baseline
mean_StD<-
  my_basline%>% #dataframe in typical ggplot formate
  group_by(label_y)%>%
  summarise(mean=mean(value),sd=sd(value))

#add the baseline to the rest to process with ggplot later
combin<-
  full_join(for_bubble_plot, mean_StD)

#lets see who changes from baseline
updown_regulation<-
  combin%>%
  filter(
    str_detect(label_x,"s_")
  )%>%
  filter(
    !str_detect(label_y,"zother")&
      !str_detect(label_y,"Other")&
      !str_detect(label_y,"Relative")
  )%>%
  mutate(
    change=case_when(
      value > mean+sd ~ "up",
      value < mean-sd ~ "down",
      TRUE ~ "no"
    ),
    change_value=case_when(
      change == "up" ~ value-(mean+sd),
      change == "down" ~ value-(mean-sd),
      TRUE ~ 0
    )
  )

#to look at it manually
display_change<-
  updown_regulation%>%
  filter(
    change!="no"
  )%>%
  filter(
    change_value>=1 |
      change_value<=-1
  )

for_change_plot<-
  updown_regulation%>%
  filter(
    label_y %in% unique(display_change$label_y)
  )%>%
  mutate(
    label_x=str_replace(
      label_x,
      "Ast",
      "st"
      )
  )

#make change plot
change_p<-
  ggplot(for_change_plot)+
  geom_point(
    aes(
      x=label_x,
      y=label_y,
      shape=change,
      color=season,
      fill=season
      )
    )+
  scale_shape_manual(values = c("up"=24,"down"=25,"no"=16))+
  scale_color_manual(
    values = c(
      "winter"="#1A0E34", 
      "twilight"="#929489", 
      "spring"="#A7CE9D", 
      "summer/fall"="#6E3E67"
    )
  )+
  scale_fill_manual(
    values = c(
      "winter"="#1A0E34", 
      "twilight"="#929489", 
      "spring"="#A7CE9D", 
      "summer/fall"="#6E3E67"
    )
  )+
  scale_size_continuous()+
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 90,size = 7,vjust = 0.3),
    axis.text.y = element_text(size = 7, vjust = 0.2),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    strip.text = element_text(size = 8),
    strip.placement = "inside",
    strip.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.border = element_rect(fill=NA,color = "black"),
    legend.position="bottom",
    legend.text = element_text(size = 8),
    axis.ticks.y = element_line(color = "black", size = 0.5),
    axis.ticks.x = element_line(color = "black", size = 0.5),
    aspect.ratio = 0.5/1
  )

change_p
  
#now I also want to cluster in a tree so I need a distance matrix 
dist_bray<-
  pct%>%
  select(
    contains("s_")
  )%>%
  t()%>%
  vegdist(method = "bray")

#compute cluster
tree_cluster<-
  hclust(dist_bray)

plot(tree_cluster)

#following is to make the tree look nice
season_c<-
  tree_cluster%>%
  labels()%>%
  str_remove(.,"sub20000")%>%
  str_remove(.,"[s,w]_st.._")%>%
  str_remove(.,"_[5,3]")%>%
  str_replace(.,"201[7,8]12[1,2][0,7]","#1A0E34")%>%
  str_replace(.,"201[8,9]0[4,5][0,2][1,5]","#A7CE9D")%>%
  str_replace(.,"20180228","#929489")%>%
  str_replace(.,"20190913","#6E3E67")%>%
  str_remove(.,"\\.")


label_plot1<-
  labels(tree_cluster)%>%
  str_remove(.,"sub20000\\.")%>%
  str_remove(.,"_st[0,2][3,5,6,7]")

label_plot2<-
  labels(tree_cluster)%>%
  str_remove(.,"sub20000")%>%
  str_extract(.,"_st[0,2][3,5,6,7]")

label_plot3<-
  cbind(
    label_plot1,
    label_plot2
    )%>%
  as.data.frame()%>%
  unite(
    label_plot,
    label_plot1,
    label_plot2
  )%>%
  mutate(
    lab=str_replace(label_plot,"__","_")
  )
  
label_plot<- as.vector(label_plot3$lab)

#nice plotting of the tree
ggd2<- 
  tree_cluster%>%
  as.dendrogram()%>%
  set("labels",
      label_plot
  )%>%
  set("branches_k_color", 
      value = c("brown"), k = 1
  )%>%
  set("leaves_pch", 19) %>%  # node point type
  set("leaves_cex", 1) %>%  # node point size
  set("leaves_col", season_c)%>% # node point color
  set("labels_cex", NA)

tree_p <- 
  as.ggplot(function() plot(ggd2)) 

tree_p

#now we reorder the changeplot based on the tree clustering
change_p<-
  ggplot(for_change_plot)+
  geom_point(
    aes(
      x=factor(
        label_x,
        levels = label_plot
          ),
      y=label_y,
      shape=change,
      color=season,
      fill=season
    ),size=2.5
  )+
  scale_shape_manual(values = c("up"=24,"down"=25,"no"=16))+
  scale_color_manual(
    values = c(
      "winter"="#1A0E34", 
      "twilight"="#929489", 
      "spring"="#A7CE9D", 
      "summer/fall"="#6E3E67"
    )
  )+
  scale_fill_manual(
    values = c(
      "winter"="#1A0E34", 
      "twilight"="#929489", 
      "spring"="#A7CE9D", 
      "summer/fall"="#6E3E67"
    )
  )+
  scale_size_continuous()+
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 90,size = 7,vjust = 0.3),
    axis.text.y = element_text(size = 7, vjust = 0.2),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    strip.text = element_text(size = 8),
    strip.placement = "inside",
    strip.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.border = element_rect(fill=NA,color = "black"),
    legend.position="bottom",
    legend.text = element_text(size = 8),
    axis.ticks.y = element_line(color = "black", size = 0.5),
    axis.ticks.x = element_line(color = "black", size = 0.5),
    aspect.ratio = 0.5/0.9
  )

change_p

#now combine bothe, still needs some manual editing since both are made by same plotting functions but as good as it gets 
treeandchange<-
  ggarrange(
    tree_p,
    change_p,
    nrow = 2,
    ncol = 1,
    align = "v",
    heights = c(0.6,1)
  )

treeandchange

pdf("changeplot.pdf",paper = "a4r")
treeandchange
dev.off()


#plotting ASV data####
#read in ASV data from DADA2
imput_asv<-
  read_delim("SVA_ASVs_rawcounts.csv",delim = ";")%>%
  drop_na()

#calculate relative abundance
pct_asvs<-
  imput_asv%>%
  mutate_at(.vars = vars(which(sapply(.,class)=="numeric")),.funs = get.pct)%>%
  mutate(
    tax=make.unique(tax)
  )%>%
  as.data.frame()%>%
  melt()

pct_asvs_context<-
  pct_asvs%>%
  mutate(
    season=case_when(
      str_detect(variable,"201[78]12") ~ "winter",
      str_detect(variable,"201[89]0[45]") ~ "spring",
      str_detect(variable,"201909") ~ "summer/fall",
      str_detect(variable,"201802") ~ "twilight",
      T ~ "error"
    ),
    station=str_extract(variable,"st[:digit:]{2}"),
    biom=case_when(
      str_detect(variable,"w\\.") ~ "water",
      T ~ "sediment"
    ),
    year=str_extract(variable,"201[789]"),
    ASV=case_when(
      str_detect(tax, "_[:digit:]+$") ~ str_extract(tax,"[:digit:]+$"),
      str_detect(tax, "_[:alpha:]+$") ~ str_extract(tax,"[:alpha:]+_[:alpha:]+$"),
      TRUE ~ str_extract(tax,"[:alpha:]+_[:alpha:]+\\.[:digit:]+$")
    )
  )

sediment<-
  pct_asvs_context%>%
  filter(
    biom=="sediment"
  )

water<-
  pct_asvs_context%>%
  filter(
    biom=="water"
  )%>%
  mutate(
    seasonforwater=case_when(
      season=="winter" ~ "winter",
      season=="spring"&year=="2018" ~ "spring2018",
      season=="spring"&year=="2019" ~ "spring2019",
      season=="summer/fall" ~ "summer/fall",
      season=="twilight" ~ "twilight",
      T ~ "error"
    )
  )

list_sediment<-
  sediment%>%
  group_by(season)%>%
  group_split()

list_water<-
  water%>%
  mutate(
    season=seasonforwater
  )%>%
  group_by(season)%>%
  group_split()

test<-
  list_sediment[[1]]

get.mean.SD.of.group(test)%>%view()

sediment_seasonal_ASV_means<-
  map(list_sediment,get.mean.SD.of.group)%>%
  reduce(full_join)%>%
  mutate(
    min=mean-sd,
    max=mean+sd,
    group=case_when(
      str_detect(tax, "_Gammaproteobacteria") ~ "Gammaproteobacteria",
      str_detect(tax, "_Deltaproteobacteria") ~ "Deltaproteobacteria",
      str_detect(tax, "Bacteroidota_") ~ "Bacteroidetes",
      str_detect(tax, "Planctomycetota_") ~ "Planctomycetes",
      str_detect(tax, "Verrucomicrobiota_") ~ "Verrucomicrobia",
      str_detect(tax, "Actinobacteriota_") ~ "Actinobacteria",
      T ~ "other"
    )
  )


water_seasonal_ASV_means<-
  map(list_water,get.mean.SD.of.group)%>%
  reduce(full_join)%>%
  mutate(
    min=mean-sd,
    max=mean+sd,
    group=case_when(
      str_detect(tax, "_Gammaproteobacteria") ~ "Gammaproteobacteria",
      str_detect(tax, "_Deltaproteobacteria") ~ "Deltaproteobacteria",
      str_detect(tax, "Bacteroidota_") ~ "Bacteroidetes",
      str_detect(tax, "Planctomycetota_") ~ "Planctomycetes",
      str_detect(tax, "Verrucomicrobiota_") ~ "Verrucomicrobia",
      str_detect(tax, "Actinobacteriota_") ~ "Actinobacteria",
      T ~ "other"
    )
  )

top100pctsed_p<-
  sediment_seasonal_ASV_means%>%
  filter(
    str_detect(season,"spring")|
      str_detect(season,"winter")
  )%>%
  get.asw.plots3()

top100pctwat_p<-
  water_seasonal_ASV_means%>%
  filter(
    str_detect(season,"spring")|
      str_detect(season,"winter")
  )%>%
  get.asw.plots4()

pdf("pct_ASVs_top100_sed_wat.pdf")
top100pctsed_p
top100pctwat_p
dev.off()

list_taxa_sed<-
  sediment_seasonal_ASV_means%>%
  filter(
    str_detect(season,"spring")|
      str_detect(season,"winter")
  )%>%
  group_by(group)%>%
  group_split()

test<-
  list_taxa_sed[[1]]

get.asw.plots.onlisttaxa(DF = test)

taxa_plot_list<-
  map(list_taxa_sed,get.asw.plots.onlisttaxa)

pdf("taxon_likeFISH_top100.pdf")
taxa_plot_list
dev.off()

for_bubble_plot<-
  list_sediment%>%
  reduce(full_join)%>%
  mutate(
    variable=str_remove(variable, "neb\\."),
    plotTax=str_remove(tax,"_[:digit:]+$")
  )%>%
  separate(
    variable, 
    into = c("U1","U2"),
    sep = "\\.",
    remove = F
  )%>%
  unite(
    plotID,
    biom,
    U2,
    U1,
    remove = F
  )

#now just to see if the ASV annotation shows a similar picture than the SILVA one
for_bubble_plot_reduced<-
  for_bubble_plot%>%
  select(
    plotID,
    plotTax,
    value
  )%>%
  unite(
    ID,
    plotID,plotTax,
    sep = ":"
  )%>%
  group_by(ID)%>%
  dplyr::summarise_all(.funs = sum)%>%
  group_by(ID)%>%
  filter(any(value>2))%>%
  separate(
    ID,
    into = c("plotID","plotTax"),
    sep = ":"
  )%>%
  mutate(
    season=case_when(
      str_detect(plotID,"201[78]12") ~ "winter",
      str_detect(plotID,"201[89]0[45]") ~ "spring",
      str_detect(plotID,"201909") ~ "summer/fall",
      str_detect(plotID,"201802") ~ "twilight",
      T ~ "error"
    )
  )

p_bubble<-
  ggplot(for_bubble_plot_reduced)+
  geom_point(aes(
    x=plotID,
    y=plotTax,
    size=value,
    color=season
  )
  )+
  scale_size_continuous(range = c(-1, 10),breaks = c(0,2,5,10,15,20,40,60),name = "[% relativ abundance]")+
  scale_color_manual(
    values = c("winter"="#111E2C", 
               "twilight"="#E6E5FF", 
               "spring"="#8D8556", 
               "summer/fall"="#575134"
    )
  )+
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 90,size = 7),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    strip.text = element_text(size = 8),
    strip.placement = "inside",
    strip.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.border = element_rect(fill=NA,color = "black"),
    legend.position="bottom",
    legend.text = element_text(size = 8),
    axis.ticks.y = element_line(color = "black", size = 1),
    axis.ticks.x = element_line(color = "black", size = 1),
    aspect.ratio = 1/1.8
  )

p_bubble

pdf("bubble_ASV_taxa_biger2pct_lisbon.pdf",width = 16.5,height = 11.7 )
p_bubble
dev.off()
