library("seqtime")
library("vegan")
library("ggplot2")
library("ggpubr")
setwd("/Users/xu-wenwang/Dropbox/Projects/FEAST/code")

################### parameter of glv ########################
M = 5
dis = c()
types = c()
sink_source = c()

for (unique_species in c(1,3,5,10)){
  set.seed(123)
  C = 0.5
  N = 90
  N_sub = N/M
  A = matrix(rnorm(N*N,mean=0,sd=1), N, N)
  set.seed(123)
  Connect = sample(0:1, N^2, prob = c(1-C, C), replace = TRUE)
  A = A * matrix(Connect, N, N)
  diag(A) = -5*(C>0)
  set.seed(123)
  b = runif(N)
  b[b>0] = 0.5
  
  steady_source = matrix(0, M, N)
  meta_source = c()
  
  # five sources
  set.seed(123)
  collection_all = sample(1:N,N)
  for (i in 1:M){
    set.seed(i)
    collection_1 = collection_all[((i-1)*unique_species+1):(i*unique_species)]
    collection_1 = c(collection_1, collection_all[(M*unique_species+1):N])
    #collection_1 = sample(1:N,60)
    y_1 = runif(N)
    y_1[setdiff(1:N,collection_1)] = 0
    x = glv(N = N, A, b = b, y = y_1, tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
    steady_source[i,] = x[,ncol(x)]/sum(x[,ncol(x)])
  }
  steady_source[steady_source<0] = 0
  
  ################### construct the sinks ################
  M_s = 100
  steady_sink = matrix(0, M_s, N)
  
  i = 1
  while (i<=M_s){
    print(i)
    mix_order = sample(1:M,M)
    m = rep(1,M)/M
    y = steady_source[mix_order[1],]*m[1]
    y = glv(N = N, A, b = b, y = y, tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
    y[y<0] = 0
    y = y[,ncol(y)]
    success = 0
    for (j in 2:M){
      tryCatch({x = glv(N = N, A, b = b, y = y, tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)},warning=function(w) success <<-1)
      x[x<0] = 0
      y = steady_source[mix_order[j],]*m[j] + x[,ncol(x)]
    }
    if (success==0){
      steady_sink[i,] = x[,ncol(x)]/sum(x[,ncol(x)])
      i = i + 1}
  }
  
  dis_sub_1 = vegdist(steady_sink)
  dis_sub_2 = vegdist(steady_source)

  dis=c(dis,dis_sub_1,dis_sub_2)
  types = c(types,rep(unique_species,length(dis_sub_1)),rep(unique_species,length(dis_sub_2)))
  sink_source = c(sink_source,rep('sink',length(dis_sub_1)),rep('source',length(dis_sub_2)))
}
dis_merged = data.frame(dissimilarity=dis,types=factor(types),sink_source=sink_source)

g1 = ggboxplot(dis_merged,x='types',y='dissimilarity',
               color = "sink_source",width = 0.5)+
  scale_color_manual(values = c("#A20056FF", "#631879FF"))+
  stat_compare_means(aes(group = sink_source), label = "p.signif",label.y=0.45)+ylim(0,0.5)+
  theme_bw()+ylab('Bray-Curtis dissimilarity')+xlab("Numer unique species")+
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.background = element_rect(fill = "white", size = 4, colour = "white"),
    legend.justification = c(0, 1),
    axis.ticks = element_line(colour = "black", size = 0.2),
    axis.text.x = element_text(color="black", size=10),
    axis.text.y = element_text(color="black", size=10)
  )

ggsave(g1,file="../figures/glv_priority_dissmilarity.pdf",width=3.5, height=2.8)

