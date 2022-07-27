library("seqtime")
library("FEAST")
library("ggplot2")
library(scales)
library(jcolors)
library(randomForest)
setwd("/Users/xu-wenwang/Dropbox/Projects/FEAST/code")

################### parameter of glv ########################
realization = 3
M_s = 100
depth = 1000
method = c()
corr_final = c()
strength_final = c()
sd_final = c()

# vary the identity
for (sigma in c(0,0.00001,0.0001,0.001,0.01,0.1)){
  set.seed(123)
  N = 90
  M = 3
  C = 0.5
  A = matrix(rnorm(N*N,mean=0,sd=sigma), N, N)
  Connect = sample(0:1, N^2, prob = c(1-C, C), replace = TRUE)
  A = A * matrix(Connect, N, N)
  diag(A) = -sigma*5
  b = runif(N)
  b[b>0] = 0.5
  
  ################### construct sources ################
  corr_feast = c()
  corr_sourcetracker = c()
  corr_randomforest = c()
  JSD_all = c()
  for (rea in 1:realization){
    steady_source = matrix(0, M, N)
    meta_source = c()
    
    collection_all = sample(1:N,N)
    collection_1 = collection_all[1:30];
    collection_2 = collection_all[31:60];
    collection_3 = collection_all[61:90];
    
    
    ## first source
    y_1 = runif(N)
    y_1[setdiff(1:N,collection_1)] = 0
    x = glv(N = N, A, b = b, y = y_1, tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
    steady_source[1,] = x[,ncol(x)]
    meta_source = rbind(meta_source, c(paste("Adult",1, sep = ""), 'Source', 1))
    
    ## second source
    y_2 = runif(N)
    y_2[setdiff(1:N,collection_2)] = 0
    x = glv(N = N, A, b = b, y = y_2, tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
    steady_source[2,] = x[,ncol(x)]
    meta_source = rbind(meta_source, c(paste("Adult",2, sep = ""), 'Source', 2))
    
    ## third source
    y_3 = runif(N)
    y_3[setdiff(1:N,collection_3)] = 0
    x = glv(N = N, A, b = b, y = y_3, tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
    steady_source[3,] = x[,ncol(x)]
    meta_source = rbind(meta_source, c(paste("Adult",3, sep = ""), 'Source', 3))
    
    
    steady_source[steady_source<0] = 0
    steady_source[1,] = steady_source[1,] / sum(steady_source[1,])
    steady_source[2,] = steady_source[2,] / sum(steady_source[2,])
    steady_source[3,] = steady_source[3,] / sum(steady_source[3,])
    
    ################### construct the sinks ################
    steady_sink = matrix(0, M_s, N)
    true_source = matrix(0, M_s, M)
    meta_sink = c()
    
    for (i in 1:M_s){
      m = runif(M)
      m = m/sum(m)
      y = rep(0, N)
      for (j in 1:M){ y = y + m[j]*steady_source[j,]}
      x = glv(N = N, A, b = b, y = y, tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
      x[x<0]=0
      steady_sink[i,] = x[,ncol(x)]/sum(x[,ncol(x)])
      true_source[i,] = m
      meta_sink = rbind(meta_sink, c(paste("infant",i, sep = ""), 'Sink', i))
    }
    
    count_table = as.data.frame(rbind(steady_source*depth, steady_sink*depth))
    colnames(count_table) = dput(paste('taxa_', seq(1,N,1), sep = ""))
    rownames(count_table) = dput(paste('ERR', seq(1,M+M_s,1), sep = ""))
    
    meta_table = data.frame(rbind(meta_source, meta_sink))
    colnames(meta_table) = c('Env','SourceSink','id')
    rownames(meta_table) = dput(paste('ERR',seq(1,M+M_s,1),sep = ""))
    
    # FEAST
    FEAST_output <- FEAST(C = as.matrix(ceiling(count_table)), metadata = meta_table, different_sources_flag = 0, dir_path = "../results/GLV_temp/0.5_0.5/",
                          outfile="GLV_vary_step")
    FEAST_output = as.matrix(as.data.frame(as.data.frame(FEAST_output[1:4])))
    rownames(FEAST_output) = paste(rownames(meta_table)[which(meta_table$SourceSink=='Sink')],meta_table$Env[which(meta_table$SourceSink=='Sink')],sep='_')
    setwd("/Users/xu-wenwang/Dropbox/Projects/FEAST/code")
    
    rownames(FEAST_output) = sub("_.*", "", rownames(FEAST_output))
    colnames(FEAST_output) = sub("_.*", "", colnames(FEAST_output))
    FEAST_output = FEAST_output[rownames(count_table)[(M+1):nrow(count_table)],]
    FEAST_output = FEAST_output[,c(rownames(count_table)[1:M],'Unknown')]
    
    corr_R = R2_Score(melt(true_source)$value,  melt(as.numeric(FEAST_output[,1:M]))$value)
    corr_feast = c(corr_feast, corr_R)
    for (row in 1:nrow(steady_source)) {steady_source[row,]<-steady_source[row,]/sum(steady_source[row,])}
    JSD_all = c(JSD_all, mean(jsdmatrix(steady_source)[upper.tri(jsdmatrix(steady_source))]))
    
    ############ Sourcetracker ##########################
    train.ix <- which(meta_table$SourceSink=='Source')
    test.ix <- which(meta_table$SourceSink=='Sink')
    envs <- meta_table$Env
    if(is.element('Description',colnames(meta_table))) desc <- meta_table$Description
    # load SourceTracker package
    source('./sourcetracker/src/SourceTracker.r')
    alpha1 <- alpha2 <- 0.001
    st <- sourcetracker(ceiling(count_table[train.ix,]), envs[train.ix])
    results <- predict(st,ceiling(count_table[test.ix,]), alpha1=alpha1, alpha2=alpha2)
    results_pro = results$proportions
    results_pro = results_pro[rownames(count_table)[(M+1):nrow(count_table)],]
    results_pro = results_pro[,c(meta_table$Env[1:M],'Unknown')]
    
    corr_R = R2_Score(melt(true_source)$value,  melt(as.numeric(results_pro[,1:M]))$value)
    corr_sourcetracker = c(corr_sourcetracker, corr_R)
    
    ############ random forest ########################
    train.ix <- which(meta_table$SourceSink=='Source')
    test.ix <- which(meta_table$SourceSink=='Sink')
    count_table = ceiling(count_table)
    count_table['sourcesink'] = factor(meta_table$Env)
    train_data = count_table[train.ix,]
    
    rf <- randomForest(droplevels(train_data$sourcesink)~., data=train_data)
    p2 <- predict(rf, count_table[test.ix,],type = "prob")
    corr_R = R2_Score(melt(true_source)$value,  melt(as.numeric(p2[,1:M]))$value)
    corr_randomforest = c(corr_randomforest, corr_R)
  }
  corr_final = c(corr_final, mean(corr_feast), mean(corr_sourcetracker), mean(corr_randomforest))
  strength_final = c(strength_final, sigma, sigma, sigma)
  sd_final = c(sd_final, sd(corr_feast), sd(corr_sourcetracker), sd(corr_randomforest))
  method = c(method, "FEASE", 'Source tracker','Random forest')
}
dat = data.frame(corr_final = corr_final, sd_final = sd_final/sqrt(realization), strength_final = (strength_final), method = method)
dat1 = dat
dat1$corr_final[dat1$corr_final<0] = 0
dat1$sd_final[dat$corr_final<0] = 0

g1<- ggplot(dat1, aes(x=strength_final+0.000001, y=corr_final,color=method)) +
  geom_line() + geom_point(aes(shape=method),size=3)+
  scale_shape_manual(values=c(15, 16, 17)) +
  geom_errorbar(aes(ymin = corr_final-sd_final, ymax = corr_final+sd_final, color=method), width=0.001)+
  scale_color_jcolors(palette = "pal5") + theme_bw()+xlab("Interaction strength") + ylab(R^2~'of proportion estimates')+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+annotation_logticks() + 
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.background = element_rect(fill = "white", size = 4, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = 'none',
    axis.ticks = element_line(colour = "black", size = 0.2),
    panel.grid.major = element_line(colour = "grey70", size = 0.2),
    panel.grid.minor = element_blank(),axis.text.x = element_text(color="black", size=10),
    axis.text.y = element_text(color="black", size=10)
  )

ggsave(g1,file="../figures/GLV_linear_interaction_sigma_relative.pdf",width=3, height=2.8)
write.table(dat, file = "../results/GLV/R_linear_interaction_sigma_relative.csv", row.names = FALSE, col.names = FALSE,sep=",")

