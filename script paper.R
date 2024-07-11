library(ChAMP)
library(MethylResolver)
library(survminer)
library(survival)
library(biomaRt)
library(data.table)
library(openxlsx)
library(readxl)


testDir<-paste0(getwd(), "/copiati")
myLoad_v1<-champ.load(directory = testDir, method = "minfi", force=TRUE, arraytype = "EPIC")
source("S:/SARA/IFH/champ.load.v2.R")
#this function is an adaptation of champ.load that handles EPICv2 arrays
testDir<-paste0(getwd(), "/copiativ2")
myLoad_v2<-champ.load.v2(directory = testDir, filterXY = FALSE)
load("S:/SARA/IFH/sonde_comuni_v2_v1.R")
#this file contains probes witht the same genomic location in both arraytypes
beta_v1 <- myLoad_v1$beta
beta_v2 <- myLoad_v2$beta
row.names(beta_v2) <- sub("^_*(.*?)_.*","\\1",row.names(beta_v2))
merge_betas <- merge(beta_v1, beta_v2, by.x=0, by.y=0)
rownames(merge_betas) <- merge_betas$Row.names
merge_betas$Row.names <- NULL
merge_def <- merge(sonde_comuni_v2_v1, merge_betas, by.x="Probe_ID", by.y=0)
rownames(merge_def) <- merge_def$Probe_ID
merge_def <- merge_def[,c(-1:-6)]
merge_def <- as.matrix(merge_def)
myNorm <- champ.norm(beta=merge_def, arraytype = "EPIC", cores = 8)
myNorm<-na.omit(myNorm)
MethylResolver(methylMix = myNorm)
methy<-read.table("MethylResolver.txt")
methy<-methy[,5:27]

km<-read_xlsx("casistica Ewing aggiornata aprile 2024.xlsx")
km<-merge(km,methy,by.x="ID",by.y="row.names")
res<-t(methy)
res<-res[,1:3]
View(res)
res<-as.data.frame(res)
colnames(res)<-c("cutoff", "bonferroni_pvalue", "better_prognosis")
res$cutoff<-NA
res$bonferroni_pvalue<-NA
res$better_prognosis<-NA
res$types<-row.names(res)
km$os_years<-km$`Overall survival (months)`/12
km$event<-gsub("ALIVE", FALSE, km$`Alive/Dead`)
km$event<-gsub("DOD", TRUE, km$event)
km$event<-as.logical(km$event)

types<-rownames(res)
#some throw errors, try them all and eliminate them before the cycle
types<-types[-c(3,6,11,14,17,22)]

for (i in types) {
#find expressio cutoff
  s<-surv_cutpoint(km,time = "os_years",event = "event",variables=i,minprop = 0.1,progressbar = TRUE)
  res$cutoff[res$types==i]<-as.numeric(s$cutpoint$cutpoint)
#assign new categories
  res.cat <- surv_categorize(s)
  colnames(res.cat)[3]<-"type"
#calculate survival pvalue between the two groups
  p<-pairwise_survdiff(Surv(os_years, event) ~type, data = res.cat, p.adjust.method = "bonferroni")
  res$bonferroni_pvalue[res$types==i]<-as.numeric(p$p.value)
#fit the curve and compute the mean survival probability for each group
 fit <- survfit(Surv(os_years, event) ~type, data = res.cat)
 sum<-surv_summary(fit)
 mean_p_high<-mean(sum$surv[sum$strata %like% "high"])
 mean_p_low<-mean(sum$surv[sum$strata %like% "low"])
 if(mean_p_high>mean_p_low) {res$better_prognosis[res$types==i]<-"high"}
 else {if(mean_p_high==mean_p_low){res$better_prognosis[res$types==i]<-"equal"}
   else{res$better_prognosis[res$types==i]<-"low"}
 }
#Only draw the curves where the comparisons are significant
if (as.numeric(p$p.value)<=0.05) {
  svg(filename = paste0("./destinazione grafici KM/",i,".svg"), width = 5, height = 5)
  plot<-ggsurvplot(fit, data = res.cat, xlab="Overall survival (years)", ylab= "Survival probability", pval = TRUE, conf.int = TRUE, ggtheme = theme_classic2(base_size=12, base_family = "serif"),font.family = "serif")
  print(plot)
  dev.off()
  }
}
write.xlsx(km, file = "casistica con risultati deconvoluzione.xlsx", rowNames=TRUE) #Supplementary Table 1
write.xlsx(res, file = "risultati KM scan con strati.xlsx", rowNames=TRUE) #Supplementary Table 2


#prepare data for boxplots
library(tidyr)
df<-km[,c(1,18:41)]
reshaped_df <- df %>%
  gather(key = "numeric_column", value = "numeric_value", -ID, -Translocation_type)
colnames(reshaped_df)[2:4]<-c("Translocation", "Type", "Abundance")
data_rel<-reshaped_df[1:352,]
data_abs<-reshaped_df[353:736,]

svg(filename = "boxplot_abs_whole.svg", width = 8, height = 4)
ggplot(data_abs, aes(x=Type, y=Abundance)) + geom_boxplot() + scale_y_continuous(limits = c(0, 1))+ scale_y_break(c(0.25, 0.75)) +theme_light()+ theme(text=element_text(size=16,  family="serif"), axis.text.x = element_text(angle=45, vjust=1, hjust=1))
dev.off()

svg(filename = "boxplot_abs_split.svg", width = 8, height = 4)
ggplot(data_abs, aes(x=Type, y=Abundance, fill=Translocation)) + geom_boxplot() + scale_y_continuous(limits = c(0, 1))+ scale_y_break(c(0.25, 0.75)) +scale_fill_grey() + theme_light()+ theme(text=element_text(size=16,  family="serif"), axis.text.x = element_text(angle=45, vjust=1, hjust=1))+ theme(legend.position = "none")
dev.off()


svg(filename = "boxplot_rel_whole.svg", width = 8, height = 3)
ggplot(data_rel, aes(x=Type, y=Abundance)) + geom_boxplot() + scale_y_continuous(limits = c(0, 1)) +theme_light()+ theme(text=element_text(size=16,  family="serif"), axis.text.x = element_text(angle=45, vjust=1, hjust=1))
dev.off()

svg(filename = "boxplot_rel_split.svg", width = 8, height = 3)
ggplot(data_rel, aes(x=Type, y=Abundance, fill=Translocation)) + geom_boxplot() + scale_y_continuous(limits = c(0, 1)) +scale_fill_grey() + theme_light()+ theme(text=element_text(size=16,  family="serif"), axis.text.x = element_text(angle=45, vjust=1, hjust=1))+ theme(legend.position = "none")
dev.off()