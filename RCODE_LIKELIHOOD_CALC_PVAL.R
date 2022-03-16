library(tidyverse)

argv=commandArgs(T)
DATA=argv[1]
GENO=argv[2]
SG_ID=argv[3]
ID=argv[4]
FREQ_FILE=argv[5]
E<-as.numeric(argv[6])
ENAME<-argv[7]
NPERM<-as.numeric(argv[8])

geno<-read_tsv(gzfile(GENO))
FREQ<-as_tibble(read.table(gzfile(FREQ_FILE),header=F))
METAGE_data<-read_tsv(DATA)
sample_list<-colnames(geno[5:(ncol(geno))])
result<-tibble(ID=sample_list)

FREQ<-dplyr::select(FREQ,V2,V4,V5,V6)
colnames(FREQ)<-c("SNP","ALT","REF","ALT_FREQ")

data<-METAGE_data%>%
filter(REF_CT+ALT_CT>=1)%>%
mutate(REF_REF=log10(1-E)*REF_CT+log10(E)*ALT_CT)%>%
mutate(REF_ALT=log10(1/2)*(REF_CT+ALT_CT))%>%
mutate(ALT_ALT=log10(E)*REF_CT+ log10(1-E)*ALT_CT)

#calculate per sample likelihood score
#sum up genotype probability
res<-c()
for(sample in sample_list){
    WGS<-dplyr::select(geno,all_of(c("CHROM","POS","REF","ALT",sample)))
    colnames(WGS)<-c("CHROM","POS","REF","ALT","GT")
    P<-data %>%
    dplyr::select(CHROM,POS,REF,ALT,REF_REF,REF_ALT,ALT_ALT)%>%
    inner_join(WGS,by=c("CHROM","POS","REF","ALT"))%>%
    mutate(P=case_when(GT=="0/0"~REF_REF,
    GT=="0/1"~REF_ALT,
    GT=="1/1"~ALT_ALT))%>%
    pull(P)
    Score<-sum(P)
    res<-c(res,Score)
}

result$Score<-res

if(SG_ID %in% result$ID){
    ANS<-filter(result,ID==SG_ID)%>%
    pull(Score)
}else{
    ANS=NA
}

NUM1<-nrow(data)
freq_data<-inner_join(data,FREQ,by=c("SNP","REF","ALT"))
NUM2<-nrow(freq_data)

if(NUM1==NUM2){
    
freq_data<-mutate(freq_data,EXP=
(ALT_FREQ*ALT_FREQ*ALT_ALT+
2*ALT_FREQ*(1-ALT_FREQ)*REF_ALT+
(1-ALT_FREQ)*(1-ALT_FREQ)*REF_REF
))%>%
mutate(VAR=
(ALT_FREQ*ALT_FREQ*(ALT_ALT-EXP)*(ALT_ALT-EXP)+
2*ALT_FREQ*(1-ALT_FREQ)*(REF_ALT-EXP)*(REF_ALT-EXP)+
(1-ALT_FREQ)*(1-ALT_FREQ)*(REF_REF-EXP)*(REF_REF-EXP)
))

EXP<-sum(freq_data$EXP)
SD<-sqrt(sum(freq_data$VAR))

p_result<-c()
for(i in seq(1,NPERM)){
freq_data$PERM<-runif(nrow(freq_data),0,1)
freq_data<-mutate(freq_data,Score=case_when(
    (PERM<=ALT_FREQ*ALT_FREQ)~ALT_ALT,
    ((ALT_FREQ*ALT_FREQ<PERM)&(PERM<=ALT_FREQ*(2-ALT_FREQ)))~REF_ALT,
    (ALT_FREQ*(2-ALT_FREQ)<PERM)~REF_REF))

p_result<-c(p_result,sum(freq_data$Score))}

EMP_P<-(length(p_result[p_result>ANS])+1)/(NPERM+1)
EMP_E<-mean(p_result)
EMP_SD<-sd(p_result)

write_csv(
tibble(Empirical_P=EMP_P,
Empirical_EXP=EMP_E,
Empirical_SD=EMP_SD,
Analytic_P=pnorm((ANS-EXP)/SD,lower.tail=F),
Analytic_EXP=EXP,
Analytic_SD=SD,
USED_BASE_NUM=(sum(freq_data$REF_CT)+sum(freq_data$ALT_CT))
),str_c(ID,"_likelihood_p_val_summary.txt")
)

write_tsv(tibble(PERM_result=p_result),
str_c(ID,"_permutation_result.txt")
)

}

empirical<-function(x){
    return((length(p_result[p_result>x])+1)/(NPERM+1))
}

result$EMP_P<-unlist(lapply(result$Score,empirical))
result$ANA_P<-pnorm((result$Score-EXP)/SD,lower.tail=F)

tmp<-arrange(result,-Score)
TOP_Score<-tmp$Score[1]
tmp$RANK<-seq(1,nrow(tmp))
tmp<-dplyr::select(tmp,ID,RANK)
result<-left_join(result,tmp,by="ID")
result<-mutate(result,ID_MATCH=ifelse(ID==SG_ID,"Yes","No"))

write_tsv(result,
str_c(ID,"_likelihood_p_val_result.txt")
)
