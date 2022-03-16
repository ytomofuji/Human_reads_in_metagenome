library(tidyverse)

argv=commandArgs(T)
DATA=argv[1]
FREQ_DIR=argv[2]
ID=argv[3]

METAGE_data<-read_tsv(gzfile(DATA))

FREQ_TEMP<-dplyr::select(METAGE_data,CHROM,POS,REF,ALT)%>%
rename(CHR=CHROM)
FREQ<-tibble()
for(i in seq(1,22)){
tmp_FREQ<-read_tsv(gzfile(str_c(FREQ_DIR,"/ALL_POP.chr",i,".freq.chr.pos.gz")))%>%
inner_join(FREQ_TEMP,by=c("CHR","POS","REF","ALT"))
FREQ<-bind_rows(FREQ,tmp_FREQ)
gc()
}

METAGE_data<-inner_join(METAGE_data,FREQ,by=c("CHROM"="CHR","POS","REF","ALT"))

USED_BASE_NUM<-METAGE_data%>%
mutate(NUM=REF_CT+ALT_CT)%>%
pull(NUM)
USED_BASE_NUM<-sum(USED_BASE_NUM)


E_LIST<-c(1e-2,1e-3,1e-4,1e-5,1e-6,1e-7)
ENAME_LIST<-c("1e_2","1e_3","1e_4","1e_5","1e_6","1e_7")

total_summary<-tibble()

for(i in seq(1,6)){
E<-E_LIST[i]
data<-METAGE_data%>%
filter(REF_CT+ALT_CT>=1)%>%
mutate(AMR_LIK=log10(
    AMR*AMR*(E**REF_CT)*((1-E)**ALT_CT)+
    2*(1-AMR)*AMR*((1/2)**(REF_CT+ALT_CT))+
    (1-AMR)*(1-AMR)*(E**ALT_CT)*((1-E)**REF_CT)
))%>%
mutate(EUR_LIK=log10(
    EUR*EUR*(E**REF_CT)*((1-E)**ALT_CT)+
    2*(1-EUR)*EUR*((1/2)**(REF_CT+ALT_CT))+
    (1-EUR)*(1-EUR)*(E**ALT_CT)*((1-E)**REF_CT)
))%>%
mutate(AFR_LIK=log10(
    AFR*AFR*(E**REF_CT)*((1-E)**ALT_CT)+
    2*(1-AFR)*AFR*((1/2)**(REF_CT+ALT_CT))+
    (1-AFR)*(1-AFR)*(E**ALT_CT)*((1-E)**REF_CT)
))%>%
mutate(EAS_LIK=log10(
    EAS*EAS*(E**REF_CT)*((1-E)**ALT_CT)+
    2*(1-EAS)*EAS*((1/2)**(REF_CT+ALT_CT))+
    (1-EAS)*(1-EAS)*(E**ALT_CT)*((1-E)**REF_CT)
))%>%
mutate(SAS_LIK=log10(
    SAS*SAS*(E**REF_CT)*((1-E)**ALT_CT)+
    2*(1-SAS)*SAS*((1/2)**(REF_CT+ALT_CT))+
    (1-SAS)*(1-SAS)*(E**ALT_CT)*((1-E)**REF_CT)
))

AMR_LIK_SUM=sum(data$AMR_LIK)
EUR_LIK_SUM=sum(data$EUR_LIK)
AFR_LIK_SUM=sum(data$AFR_LIK)
EAS_LIK_SUM=sum(data$EAS_LIK)
SAS_LIK_SUM=sum(data$SAS_LIK)

LIK_VEC<-c(AMR_LIK_SUM,EUR_LIK_SUM,AFR_LIK_SUM,EAS_LIK_SUM,SAS_LIK_SUM)
POPS<-c("AMR","EUR","AFR","EAS","SAS")
TOP_POP<-POPS[which(LIK_VEC==max(LIK_VEC))]

summary<-tibble(Sample_ID=ID,
THR=ENAME_LIST[i],USED_BASE_NUM=USED_BASE_NUM,
AMR_LIK=sum(data$AMR_LIK),
EUR_LIK=sum(data$EUR_LIK),
AFR_LIK=sum(data$AFR_LIK),
EAS_LIK=sum(data$EAS_LIK),
SAS_LIK=sum(data$SAS_LIK),
TOP_POP=TOP_POP
)
total_summary<-bind_rows(total_summary,summary)
}

write_tsv(total_summary,str_c(ID,"_population_likelihood_result.txt"))