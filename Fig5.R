library("ggplot2")
library("readxl")
library(reshape2)
library("ggpubr")
library(jcolors)
library(gridExtra)
library(FEAST)

setwd("/Users/xu-wenwang/Dropbox/Projects/FEAST/code")

OTU_fmt_yd = read.table("../data/Lu/OTU.csv",sep = ',',row.names = 1,header = T)

meta_table = read_excel("../data/Lu/meta.xlsx")

# replace some source names
meta_table$source1[meta_table$source1=='0901D'] = "S0901D"
meta_table$source2[meta_table$source2=='0901D'] = "S0901D"
meta_table$source3[meta_table$source3=='0901D'] = "S0901D"
meta_table$source4[meta_table$source4=='0901D'] = "S0901D"
meta_table$source1[meta_table$source1=='0901F'] = "S0901F"
meta_table$source2[meta_table$source2=='0901F'] = "S0901F"
meta_table$source3[meta_table$source3=='0901F'] = "S0901F"
meta_table$source4[meta_table$source4=='0901F'] = "S0901F"
meta_table$source1[meta_table$source1=='0904A'] = "S0904A"
meta_table$source2[meta_table$source2=='0904A'] = "S0904A"
meta_table$source3[meta_table$source3=='0904A'] = "S0904A"
meta_table$source4[meta_table$source4=='0904A'] = "S0904A"
meta_table$source1[meta_table$source1=='0904B'] = "S0904B"
meta_table$source2[meta_table$source2=='0904B'] = "S0904B"
meta_table$source3[meta_table$source3=='0904B'] = "S0904B"
meta_table$source4[meta_table$source4=='0904B'] = "S0904B"

unique_source = rownames(OTU_fmt_yd)[1:24]
overlapped = intersect(rownames(OTU_fmt_yd),c(unique_source,meta_table$sinkname))
OTU_fmt_yd = OTU_fmt_yd[overlapped,]
meta_table = meta_table[meta_table$sinkname%in%overlapped,]

Env = paste('Donor_',1:length(unique_source),sep = "")
SourceSink = rep('Source',length(unique_source))
id = c(1:length(unique_source))

# excluded sinks mixed by same sources and numerical souces
keeped = c()
for (i in 1:nrow(meta_table)){
  N_s1 = sum(!is.na(c(meta_table$source1[i],meta_table$source2[i],meta_table$source3[i],meta_table$source4[i])))
  N_s2 = sum(c(meta_table$source1[i],meta_table$source2[i],meta_table$source3[i],meta_table$source4[i])%in%unique_source)
  N_s3 = length(unique(c(meta_table$source1[i],meta_table$source2[i])))
  if (N_s1==2&(N_s2==2)&(N_s3==2)){keeped = c(keeped,i)}
}
meta_table = meta_table[keeped,]

Env = c(Env, paste('Receipent_',1:nrow(meta_table),sep = ""))
SourceSink = c(SourceSink, rep('Sink',nrow(meta_table)))
id = c(id, 1:nrow(meta_table))

metadata = data.frame(Env = Env, SourceSink = SourceSink, id = id)
rownames(metadata) = c(unique_source,meta_table$sinkname)

overlaped = intersect(rownames(metadata), rownames(OTU_fmt_yd))
metadata = metadata[overlaped,]
OTU_fmt_yd = OTU_fmt_yd[overlaped,]
OTU_fmt_yd = OTU_fmt_yd[rownames(metadata),]
OTU_fmt_yd[is.na(OTU_fmt_yd)] <- 0

metadata$id[25:nrow(metadata)] = 1:(nrow(metadata)-24)

##
FEAST_output <- FEAST(C = as.matrix(OTU_fmt_yd), metadata = metadata, different_sources_flag = 0, dir_path = "../results/FMT_YD/",
                      outfile="FMT")
setwd("/Users/xu-wenwang/Dropbox/Projects/FEAST/code")

rownames(FEAST_output) = sub("_Receipent_.*", "", rownames(FEAST_output))
colnames(FEAST_output) = sub("_Donor_.*", "", colnames(FEAST_output))
FEAST_output = FEAST_output[rownames(metadata)[25:nrow(metadata)],]

write.table(FEAST_output, file = "../results/Lu/FEAST_merged.csv", row.names = TRUE, col.names = TRUE,sep=",")

predicted_contribution = c()
true_contribution = c()
source_num = c()
types = c()
for (i in 1:nrow(FEAST_output))
{
 m1 = which(meta_table$sinkname==rownames(FEAST_output)[i])
 m2 = which(colnames(FEAST_output)==meta_table$source1[m1])
 m3 = which(colnames(FEAST_output)==meta_table$source2[m1])
 m4 = which(colnames(FEAST_output)==meta_table$source3[m1])
 m5 = which(colnames(FEAST_output)==meta_table$source4[m1])
 
 m2 = colnames(FEAST_output)[m2]
 m3 = colnames(FEAST_output)[m3]
 m4 = colnames(FEAST_output)[m4]
 m5 = colnames(FEAST_output)[m5]
 
 m2 = m2[m2%in%unique_source]
 m3 = m3[m3%in%unique_source]
 m4 = m4[m4%in%unique_source]
 m5 = m5[m5%in%unique_source]
 
 m_all = unique(c(m2,m3,m4,m5))
 source_num = c(source_num, length(m_all))
 
 predicted_contribution = c(predicted_contribution,as.numeric(FEAST_output[i,m_all]),as.numeric(FEAST_output[i,25]))
 true_contribution = c(true_contribution,rep(1/length(m_all),length(m_all)),0)
 types = c(types,'true source','true source','unknown')
}
predicted_contribution = data.frame(contribution = predicted_contribution,ture = true_contribution,types=types)
data_sub = lm(contribution~ture, data = predicted_contribution)
corr_R = summary(data_sub)$r.squared 


g1 = ggplot(predicted_contribution, aes(x= contribution,fill=types,color=types)) + 
  geom_density()+scale_color_manual(values = c('#009988','#BB5566'))+scale_fill_manual(values = c('#009988','#BB5566'))+
  theme_bw()+xlab("Predicted contribution") + ylab('Density')+
  theme(axis.text.x = element_text(color="black", size=10),axis.text.y = element_text(color="black", size=10),legend.position = "none")+
  theme(strip.background = element_blank(),strip.text.x = element_blank())

####
FEAST_output_1 = FEAST_output[rownames(FEAST_output)[source_num==2],1:(ncol(FEAST_output)-1)]
for (i in 1:nrow(FEAST_output_1)){
  largeset_two = which((FEAST_output_1[i,]>=as.numeric(sort(FEAST_output_1[i,]))[23])==TRUE)
  FEAST_output_1[i,] = 0
  FEAST_output_1[i,largeset_two] = 1
}

melted_cormat = melt(as.matrix(FEAST_output_1))
melted_cormat['ground'] = 0
for (i in 1:nrow(melted_cormat)){
  m1 = which(meta_table$sinkname==melted_cormat$Var1[i])
  if (melted_cormat$Var2[i]%in%as.character(meta_table[m1,])){melted_cormat$ground[i]=1}
}

############## split 
unique_sinks = unique(melted_cormat$Var1)
g2 = ggplot(data = melted_cormat[melted_cormat$Var1%in%unique_sinks[1:64],]) + 
  geom_point(aes(x=Var1, y=Var2,size=factor(value)),color="#228833",shape = 15) +
  geom_point(aes(x=Var1, y=Var2,size=factor(ground)),color="#EE6677",shape=20) + scale_size_manual(values=c('1'=2, '0'=NA), guide="none")+
  theme_bw()+theme(axis.text.x = element_text(color="black", size=6, angle = 90),axis.text.y = element_text(color="black", size=6),legend.position = "none")+
  theme(strip.background = element_blank(),strip.text.x = element_blank())

g3 = ggplot(data = melted_cormat[melted_cormat$Var1%in%unique_sinks[65:128],]) + 
  geom_point(aes(x=Var1, y=Var2,size=factor(value)),color="#228833",shape = 15) +
  geom_point(aes(x=Var1, y=Var2,size=factor(ground)),color="#EE6677",shape=20) + scale_size_manual(values=c('1'=2, '0'=NA), guide="none")+
  theme_bw()+theme(axis.text.x = element_text(color="black", size=6, angle = 90),axis.text.y = element_text(color="black", size=6),legend.position = "none")+
  theme(strip.background = element_blank(),strip.text.x = element_blank())

g4 = ggplot(data = melted_cormat[melted_cormat$Var1%in%unique_sinks[129:192],]) + 
  geom_point(aes(x=Var1, y=Var2,size=factor(value)),color="#228833",shape = 15) +
  geom_point(aes(x=Var1, y=Var2,size=factor(ground)),color="#EE6677",shape=20) + scale_size_manual(values=c('1'=2, '0'=NA), guide="none")+
  theme_bw()+theme(axis.text.x = element_text(color="black", size=6, angle = 90),axis.text.y = element_text(color="black", size=6),legend.position = "none")+
  theme(strip.background = element_blank(),strip.text.x = element_blank())

g5 = ggplot(data = melted_cormat[melted_cormat$Var1%in%unique_sinks[193:256],]) + 
  geom_point(aes(x=Var1, y=Var2,size=factor(value)),color="#228833",shape = 15) +
  geom_point(aes(x=Var1, y=Var2,size=factor(ground)),color="#EE6677",shape=20) + scale_size_manual(values=c('1'=2, '0'=NA), guide="none")+
  theme_bw()+theme(axis.text.x = element_text(color="black", size=6, angle = 90),axis.text.y = element_text(color="black", size=6),legend.position = "none")+
  theme(strip.background = element_blank(),strip.text.x = element_blank())

ggsave(g1,file="../figures/Lu_proportion.pdf",width=2.55, height=2.5, dpi=300)
ggsave(g2,file="../figures/Lu_relation_2_1.pdf",width=7.3, height=2.5, dpi=300)
ggsave(g3,file="../figures/Lu_relation_2_2.pdf",width=7.3, height=2.5, dpi=300)
ggsave(g4,file="../figures/Lu_relation_2_3.pdf",width=7.3, height=2.5, dpi=300)
ggsave(g5,file="../figures/Lu_relation_2_4.pdf",width=7.3, height=2.5, dpi=300)
