library("ggplot2")
library("readxl")
library("ggalluvial")
library(reshape2)
library("ggpubr")
library(jcolors)
library(gridExtra)
library(FEAST)
library(vegan)

setwd("/Users/xu-wenwang/Dropbox/Projects/FEAST/code")

OTU_fmt_yd = read_excel("../data/FMT_YD/OTU_YD.xlsx")
meta_fmt_yd = read_excel("../data/FMT_YD/meta_YD.xlsx")
OTU_fmt_yd <- as.data.frame(t(as.matrix(OTU_fmt_yd)))

# merge the meta
sampleID = c()
Env = c()
SourceSink = c()
id = c()
types = c()
ground_truth = c()
respond = c()

################## source #######################
unique_source = unique(meta_fmt_yd$Donor)
for (i in 1:nrow(meta_fmt_yd)){
  if (meta_fmt_yd$Donor[i]!="'NaN'"){
    if (meta_fmt_yd$Donor[i]%in%sampleID==FALSE){
      sampleID = c(sampleID, meta_fmt_yd$Donor[i])
      Env = c(Env, paste('Donor_',i,sep = ""))
      SourceSink = c(SourceSink, 'Source')
      id = c(id, which(unique_source==meta_fmt_yd$Donor[i]))
      types = c(types, 'Donor')
    }
  }
}

unique_pre = unique(meta_fmt_yd$Before)
for (i in 1:nrow(meta_fmt_yd)){
  if (meta_fmt_yd$Before[i]!="'NaN'"){
    if (meta_fmt_yd$Before[i]%in%sampleID==FALSE){
      sampleID = c(sampleID, meta_fmt_yd$Before[i])
      Env = c(Env, paste('Donor_',i+length(unique_source),sep = ""))
      SourceSink = c(SourceSink, 'Source')
      id = c(id, which(unique_pre==meta_fmt_yd$Before[i])+7)
      types = c(types, 'Pre')
    }
  }
}

############### sink of days ##############
id_label = 1
for (i in 1:nrow(meta_fmt_yd)){
  if (meta_fmt_yd$Days[i]!="'NaN'"){
    sampleID = c(sampleID, meta_fmt_yd$Days[i])
    Env = c(Env, 'Recipient')
    SourceSink = c(SourceSink, 'Sink')
    id = c(id, id_label)
    types = c(types, 'Days')
    ground_truth = c(ground_truth, meta_fmt_yd$Donor[i])
    respond = c(respond, meta_fmt_yd$Res[i])
    id_label = id_label + 1
  }
}

############### sink of weeks ##############
for (i in 1:nrow(meta_fmt_yd)){
  if (meta_fmt_yd$Weeks[i]!="'NaN'"){
    sampleID = c(sampleID, meta_fmt_yd$Weeks[i])
    Env = c(Env, 'Recipient')
    SourceSink = c(SourceSink, 'Sink')
    id = c(id, id_label)
    types = c(types, 'Weeks')
    ground_truth = c(ground_truth, meta_fmt_yd$Donor[i])
    respond = c(respond, meta_fmt_yd$Res[i])
    id_label = id_label + 1
  }
}

############### sink of months ##############
for (i in 1:nrow(meta_fmt_yd)){
  if (meta_fmt_yd$Months[i]!="'NaN'"){
    sampleID = c(sampleID, meta_fmt_yd$Months[i])
    Env = c(Env, 'Recipient')
    SourceSink = c(SourceSink, 'Sink')
    id = c(id, id_label)
    types = c(types, 'Months')
    ground_truth = c(ground_truth, meta_fmt_yd$Donor[i])
    respond = c(respond, meta_fmt_yd$Res[i])
    id_label = id_label + 1
  }
}

############### sink of Longterm ##############
for (i in 1:nrow(meta_fmt_yd)){
  if (meta_fmt_yd$Longterm[i]!="'NaN'"){
    sampleID = c(sampleID, meta_fmt_yd$Longterm[i])
    Env = c(Env, 'Recipient')
    SourceSink = c(SourceSink, 'Sink')
    id = c(id, id_label)
    types = c(types, 'Longterm')
    ground_truth = c(ground_truth, meta_fmt_yd$Donor[i])
    respond = c(respond, meta_fmt_yd$Res[i])
    id_label = id_label + 1
  }
}

metadata = data.frame(Env = Env, SourceSink = SourceSink, id = id, ground_truth = c(rep('NA',7+88),ground_truth), types = types, 
                      response = c(rep('NA',7+88),respond))
rownames(metadata) = sampleID
overlaped = intersect(rownames(metadata), rownames(OTU_fmt_yd))
metadata = metadata[overlaped,]
OTU_fmt_yd = OTU_fmt_yd[overlaped,]
OTU_fmt_yd = OTU_fmt_yd[rownames(metadata),]
OTU_fmt_yd[is.na(OTU_fmt_yd)] <- 0

#
FEAST_output_merged = matrix(0,nrow = 259,ncol = 95+1)
rownames(FEAST_output_merged) = rownames(metadata)[96:354]
colnames(FEAST_output_merged) = c(rownames(metadata)[1:95],'Unknown')
for (i in 96:nrow(metadata)){
  if(metadata$types[i]=='Days'){
    m2 = which(meta_fmt_yd$Days==rownames(metadata)[i])
    pre_name = meta_fmt_yd$Before[m2]
  }
  if(metadata$types[i]=='Weeks'){
    m2 = which(meta_fmt_yd$Weeks==rownames(metadata)[i])
    pre_name = meta_fmt_yd$Before[m2]
  }
  if(metadata$types[i]=='Months'){
    m2 = which(meta_fmt_yd$Months==rownames(metadata)[i])
    pre_name = meta_fmt_yd$Before[m2]
  }
  if(metadata$types[i]=='Longterm'){
    m2 = which(meta_fmt_yd$Longterm==rownames(metadata)[i])
    pre_name = meta_fmt_yd$Before[m2]
  }
  ############## feast
  OTU_fmt_yd_sub = OTU_fmt_yd[c(rownames(metadata)[1:7],pre_name,rownames(metadata)[i]),]
  metadata_sub = metadata[c(rownames(metadata)[1:7],pre_name,rownames(metadata)[i]),]
  metadata_sub$id=c(1:8,1)

  FEAST_output <- FEAST(C = as.matrix(OTU_fmt_yd_sub), metadata = metadata_sub, different_sources_flag = 0, dir_path = "../results/FMT_YD/",
                        outfile="FMT")
  FEAST_output = as.data.frame(FEAST_output[1:9])
  setwd("/Users/xu-wenwang/Dropbox/Projects/FEAST/code")
  colnames(FEAST_output) = sub("X.", "", colnames(FEAST_output))
  colnames(FEAST_output) = sub("._Donor_.*", "", colnames(FEAST_output))
  colnames(FEAST_output)[1:8] = paste0("'",colnames(FEAST_output)[1:8],"'")
  FEAST_output_merged[i-95,colnames(FEAST_output)] = as.numeric(FEAST_output[1,])
}

# write.table(FEAST_output_merged, file = "../results/FMT_YD/FEAST_seprate.csv", row.names = TRUE, col.names = TRUE,sep=",")

FEAST_output = read.table(file = "../results/FMT_YD/FEAST_seprate.csv", row.names = 1, header = TRUE,sep=",",check.names = F)
## add annotation of ground truth

##################### pcoa #########################
OTU_fmt_yd_rela = OTU_fmt_yd
for (i in 1:nrow(OTU_fmt_yd_rela)){OTU_fmt_yd_rela[i,] = OTU_fmt_yd_rela[i,]/sum(OTU_fmt_yd_rela[i,])}
endo.pcoa<-wcmdscale(vegdist(OTU_fmt_yd_rela,method='bray'))
data_pcoa = data.frame(PC1=endo.pcoa[,1],PC2=endo.pcoa[,2],types=metadata$types)
g_pcoa = ggplot(data = data_pcoa, aes(x=PC1, y = PC2,color=types))+
  geom_point(alpha=0.8)+scale_color_manual(values = c("#0072BD","gray70","#77ab30","#7d2e8d","#ff0000","#d99058"))+
  theme_bw() +
  theme(axis.text.x = element_text(color="black", size=10),axis.text.y = element_text(color="black", size=10),legend.position = 'none')+
  theme(strip.background = element_blank(),strip.text.x = element_blank())+
  xlab("PC1") + ylab("PC2")
ggsave(g_pcoa,file="../figures/pcoa_fmt_yd.pdf",width=2.7, height=2.5, dpi=300)

######################################################################################
predicted_contribution = c()
pre_contribution = c()
types = c()
sink_name = c()
for (i in 1:nrow(FEAST_output))
{
  m1 = which(rownames(metadata)==rownames(FEAST_output)[i])
  source_name = metadata$ground_truth[m1]
  if(metadata$types[m1]=='Days'){
    m2 = which(meta_fmt_yd$Days==rownames(metadata)[m1])
    pre_name = meta_fmt_yd$Before[m2]
  }
  if(metadata$types[m1]=='Weeks'){
    m2 = which(meta_fmt_yd$Weeks==rownames(metadata)[m1])
    pre_name = meta_fmt_yd$Before[m2]
  }
  if(metadata$types[m1]=='Months'){
    m2 = which(meta_fmt_yd$Months==rownames(metadata)[m1])
    pre_name = meta_fmt_yd$Before[m2]
  }
  if(metadata$types[m1]=='Longterm'){
    m2 = which(meta_fmt_yd$Longterm==rownames(metadata)[m1])
    pre_name = meta_fmt_yd$Before[m2]
  }
  predicted_contribution = c(predicted_contribution,as.numeric(FEAST_output[i,source_name])+as.numeric(FEAST_output[i,pre_name]),as.numeric(FEAST_output[i,96]))
  sink_name = c(sink_name,rep(rownames(FEAST_output)[i],2))
  types = c(types,'known','unknown')
}
predicted_contribution = data.frame(contribution = predicted_contribution,types=types,sink=sink_name)
#predicted_contribution$types = factor(predicted_contribution$types,levels = c("true source","unknown","Pre"))

g1 =  ggplot(predicted_contribution, mapping=aes_string(y = "contribution", x = "types"))+
  geom_boxplot(aes_string(colour="types", fill="types"))+
  stat_summary(geom = "crossbar", width=0.65, fatten=0, color="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x)))})+
  geom_line(aes(group=sink),colour = "gray",alpha=0.2)+ylim(0,1)+
  scale_color_manual(values = c('#009988','#BB5566','#CCBB44'))+scale_fill_manual(values = c('#009988','#BB5566','#CCBB44'))+
  stat_compare_means(paired = TRUE)+
  theme_bw()+xlab("") + ylab('Contribution')+
  theme(axis.text.x = element_text(color="black", size=10,angle = 45,hjust = 1),axis.text.y = element_text(color="black", size=10),legend.position = "none")+
  theme(strip.background = element_blank(),strip.text.x = element_blank())

  

####
FEAST_output_1 = FEAST_output[,rownames(metadata)[1:7]]
for (i in 1:nrow(FEAST_output_1)){
  largeset_two = which((FEAST_output_1[i,]>=as.numeric(sort(FEAST_output_1[i,]))[7])==TRUE)
  FEAST_output_1[i,] = 0
  FEAST_output_1[i,largeset_two] = 1
}

melted_cormat = melt(as.matrix(FEAST_output_1))
melted_cormat['ground'] = 0
for (i in 1:nrow(melted_cormat)){
  m1 = which(rownames(metadata)==melted_cormat$Var1[i])
  source_name = metadata$ground_truth[m1]
  if (melted_cormat$Var2[i]%in%as.character(source_name)){melted_cormat$ground[i]=1}
}

############## split 
unique_sinks = unique(melted_cormat$Var1)
g2 = ggplot(data = melted_cormat[melted_cormat$Var1%in%unique_sinks[1:65],]) + 
  geom_point(aes(x=Var1, y=Var2,size=factor(value)),color="#228833",shape = 15) +
  geom_point(aes(x=Var1, y=Var2,size=factor(ground)),color="#EE6677",shape=20) + scale_size_manual(values=c('1'=2, '0'=NA), guide="none")+
  theme_bw()+theme(axis.text.x = element_text(color="black", size=6, angle = 90),axis.text.y = element_text(color="black", size=6),legend.position = "none")+
  theme(strip.background = element_blank(),strip.text.x = element_blank())

g3 = ggplot(data = melted_cormat[melted_cormat$Var1%in%unique_sinks[66:130],]) + 
  geom_point(aes(x=Var1, y=Var2,size=factor(value)),color="#228833",shape = 15) +
  geom_point(aes(x=Var1, y=Var2,size=factor(ground)),color="#EE6677",shape=20) + scale_size_manual(values=c('1'=2, '0'=NA), guide="none")+
  theme_bw()+theme(axis.text.x = element_text(color="black", size=6, angle = 90),axis.text.y = element_text(color="black", size=6),legend.position = "none")+
  theme(strip.background = element_blank(),strip.text.x = element_blank())

g4 = ggplot(data = melted_cormat[melted_cormat$Var1%in%unique_sinks[131:195],]) + 
  geom_point(aes(x=Var1, y=Var2,size=factor(value)),color="#228833",shape = 15) +
  geom_point(aes(x=Var1, y=Var2,size=factor(ground)),color="#EE6677",shape=20) + scale_size_manual(values=c('1'=2, '0'=NA), guide="none")+
  theme_bw()+theme(axis.text.x = element_text(color="black", size=6, angle = 90),axis.text.y = element_text(color="black", size=6),legend.position = "none")+
  theme(strip.background = element_blank(),strip.text.x = element_blank())

g5 = ggplot(data = melted_cormat[melted_cormat$Var1%in%unique_sinks[196:259],]) + 
  geom_point(aes(x=Var1, y=Var2,size=factor(value)),color="#228833",shape = 15) +
  geom_point(aes(x=Var1, y=Var2,size=factor(ground)),color="#EE6677",shape=20) + scale_size_manual(values=c('1'=2, '0'=NA), guide="none")+
  theme_bw()+theme(axis.text.x = element_text(color="black", size=6, angle = 90),axis.text.y = element_text(color="black", size=6),legend.position = "none")+
  theme(strip.background = element_blank(),strip.text.x = element_blank())

ggsave(g1,file="../figures/FMT_YD_proportion_seprate.pdf",width=2.55, height=2.8, dpi=300)
ggsave(g2,file="../figures/YD_relation_1_seprate.pdf",width=7.8, height=2, dpi=300)
ggsave(g3,file="../figures/YD_relation_2_seprate.pdf",width=7.8, height=2, dpi=300)
ggsave(g4,file="../figures/YD_relation_3_seprate.pdf",width=7.8, height=2, dpi=300)
ggsave(g5,file="../figures/YD_relation_4_seprate.pdf",width=7.8, height=2, dpi=300)



