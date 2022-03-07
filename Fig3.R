library("seqtime")
library("FEAST")
library("vegan")
library("ggplot2")
library("ggpubr")
setwd("/Users/xu-wenwang/Dropbox/Projects/FEAST/code")
library(jcolors)
library(randomForest)
library(scales)
library(Rtsne)
library(reshape2)
################### parameter of glv ########################

M_s = 7
depth = 1000

# vary the identity
set.seed(123)
C = 0.5
N = 90
M = 3
A = matrix(rnorm(N*N,mean=0,sd=1), N, N)
Connect = sample(0:1, N^2, prob = c(1-C, C), replace = TRUE)
A = A * matrix(Connect, N, N)
diag(A) = -5*(C>0)
b = runif(N)
b[b>0] = 0.5

steady_source = matrix(0, M, N)
meta_source = c()

collection_all = sample(1:N,N)
collection_1 = collection_all[1:30]
collection_2 = collection_all[31:60]
collection_3 = collection_all[61:90]


## first source
y_1 = runif(N)
y_1[setdiff(1:N,collection_1)] = 0
x = glv(N = N, A, b = b, y = y_1, tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
steady_source[1,] = x[,ncol(x)]/sum(x[,ncol(x)])
meta_source = rbind(meta_source, c(paste("Adult",1, sep = ""), 'Source', 1))

## second source
y_2 = runif(N)
y_2[setdiff(1:N,collection_2)] = 0
x = glv(N = N, A, b = b, y = y_2, tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
steady_source[2,] = x[,ncol(x)]/sum(x[,ncol(x)])
meta_source = rbind(meta_source, c(paste("Adult",2, sep = ""), 'Source', 2))

## third source
y_3 = runif(N)
y_3[setdiff(1:N,collection_3)] = 0
x = glv(N = N, A, b = b, y = y_3, tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
steady_source[3,] = x[,ncol(x)]/sum(x[,ncol(x)])
meta_source = rbind(meta_source, c(paste("Adult",3, sep = ""), 'Source', 3))

steady_source[steady_source<0] = 0
################### construct the sinks ################
steady_sink = matrix(0, M_s, N)
true_source = matrix(0, M_s, M)
meta_sink = c()

m = c(1/3,1/3,1/3)
# order 1-2-3
for (i in 1:1){
  set.seed(i)
  y = steady_source[1,]*m[1]
  x_1 = glv(N = N, A, b = b, y = y, tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
  x_1[x_1<0]=0
  x_2 = glv(N = N, A, b = b, y = x_1[,ncol(x_1)]+steady_source[2,]*m[2], tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
  x_2[x_2<0]=0
  x_3 = glv(N = N, A, b = b, y = x_2[,ncol(x_2)]+steady_source[3,]*m[3], tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
  x_3[x_3<0]=0
  if (i==1){order_1 = cbind(x_1,x_2,x_3)}
  steady_sink[i,] = x_3[,ncol(x_3)]/sum(x_3[,ncol(x_3)])
  true_source[i,] = m
  meta_sink = rbind(meta_sink, c(paste("infant",i, sep = ""), 'Sink', i))
}
# order 1-3-2
for (i in 1:1){
  set.seed(i)
  y = steady_source[1,]*m[1]
  x_1 = glv(N = N, A, b = b, y = y, tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
  x_1[x_1<0]=0
  x_2 = glv(N = N, A, b = b, y = x_1[,ncol(x_1)]+steady_source[3,]*m[3], tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
  x_2[x_2<0]=0
  x_3 = glv(N = N, A, b = b, y = x_2[,ncol(x_2)]+steady_source[2,]*m[2], tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
  x_3[x_3<0]=0
  if (i==1){order_2 = cbind(x_1,x_2,x_3)}
  steady_sink[i+1,] = x_3[,ncol(x_3)]/sum(x_3[,ncol(x_3)])
  true_source[i+1,] = m
  meta_sink = rbind(meta_sink, c(paste("infant",i+1, sep = ""), 'Sink', i+1))
}
# order 2-1-3
for (i in 1:1){
  set.seed(i)
  y = steady_source[2,]*m[2]
  x_1 = glv(N = N, A, b = b, y = y, tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
  x_1[x_1<0]=0
  x_2 = glv(N = N, A, b = b, y = x_1[,ncol(x_1)]+steady_source[1,]*m[1], tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
  x_2[x_2<0]=0
  x_3 = glv(N = N, A, b = b, y = x_2[,ncol(x_2)]+steady_source[3,]*m[3], tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
  x_3[x_3<0]=0
  if (i==1){order_3 = cbind(x_1,x_2,x_3)}
  steady_sink[i+2,] = x_3[,ncol(x_3)]/sum(x_3[,ncol(x_3)])
  true_source[i+2,] = m
  meta_sink = rbind(meta_sink, c(paste("infant",i+2, sep = ""), 'Sink', i+2))
}
# order 2-3-1
for (i in 1:1){
  set.seed(i)
  y = steady_source[2,]*m[2]
  x_1 = glv(N = N, A, b = b, y = y, tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
  x_1[x_1<0]=0
  x_2 = glv(N = N, A, b = b, y = x_1[,ncol(x_1)]+steady_source[3,]*m[3], tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
  x_2[x_2<0]=0
  x_3 = glv(N = N, A, b = b, y = x_2[,ncol(x_2)]+steady_source[1,]*m[1], tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
  x_3[x_3<0]=0
  if (i==1){order_4 = cbind(x_1,x_2,x_3)}
  steady_sink[i+3,] = x_3[,ncol(x_3)]/sum(x_3[,ncol(x_3)])
  true_source[i+3,] = m
  meta_sink = rbind(meta_sink, c(paste("infant",i+3, sep = ""), 'Sink', i+3))
}
# order 3-1-2
for (i in 1:1){
  set.seed(i)
  y = steady_source[3,]*m[3]
  x_1 = glv(N = N, A, b = b, y = y, tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
  x_1[x_1<0]=0
  x_2 = glv(N = N, A, b = b, y = x_1[,ncol(x_1)]+steady_source[1,]*m[1], tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
  x_2[x_2<0]=0
  x_3 = glv(N = N, A, b = b, y = x_2[,ncol(x_2)]+steady_source[2,]*m[2], tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
  x_3[x_3<0]=0
  if (i==1){order_5 = cbind(x_1,x_2,x_3)}
  steady_sink[i+4,] = x_3[,ncol(x_3)]/sum(x_3[,ncol(x_3)])
  true_source[i+4,] = m
  meta_sink = rbind(meta_sink, c(paste("infant",i+4, sep = ""), 'Sink', i+4))
}
# order 3-2-1
for (i in 1:1){
  set.seed(i)
  y = steady_source[3,]*m[3]
  x_1 = glv(N = N, A, b = b, y = y, tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
  x_1[x_1<0]=0
  x_2 = glv(N = N, A, b = b, y = x_1[,ncol(x_1)]+steady_source[2,]*m[2], tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
  x_2[x_2<0]=0
  x_3 = glv(N = N, A, b = b, y = x_2[,ncol(x_2)]+steady_source[1,]*m[1], tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
  x_3[x_3<0]=0
  if (i==1){order_6 = cbind(x_1,x_2,x_3)}
  steady_sink[i+5,] = x_3[,ncol(x_3)]/sum(x_3[,ncol(x_3)])
  true_source[i+5,] = m
  meta_sink = rbind(meta_sink, c(paste("infant",i+5, sep = ""), 'Sink', i+5))
}

# simutanous mixing
for (i in 1:1){
  set.seed(i)
  y = steady_source[1,]*m[1]+steady_source[2,]*m[2]+steady_source[3,]*m[3]
  #x_1 = glv(N = N, A, b = b, y = y, tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
  #x_1[x_1<0]=0
  steady_sink[i+6,] = y/sum(y)
  true_source[i+6,] = m
  meta_sink = rbind(meta_sink, c(paste("infant",i+6, sep = ""), 'Sink', i+6))
}

steady_source[1,] = steady_source[1,] / sum(steady_source[1,])
steady_source[2,] = steady_source[2,] / sum(steady_source[2,])
steady_source[3,] = steady_source[3,] / sum(steady_source[3,])

count_table = as.data.frame(rbind(steady_source*depth, steady_sink*depth))
colnames(count_table) = dput(paste('taxa_', seq(1,N,1), sep = ""))
rownames(count_table) = dput(paste('ERR', seq(1,M+M_s,1), sep = ""))

meta_table = data.frame(rbind(meta_source, meta_sink))
colnames(meta_table) = c('Env','SourceSink','id')
rownames(meta_table) = dput(paste('ERR',seq(1,M+M_s,1),sep = ""))

# FEAST
FEAST_output <- FEAST(C = as.matrix(ceiling(count_table)), metadata = meta_table, different_sources_flag = 0, dir_path = "../results/GLV_temp/0.5_0.5/",
                      outfile="GLV_vary_step")
setwd("/Users/xu-wenwang/Dropbox/Projects/FEAST/code")

rownames(FEAST_output) = sub("_.*", "", rownames(FEAST_output))
colnames(FEAST_output) = sub("_.*", "", colnames(FEAST_output))
FEAST_output = FEAST_output[rownames(count_table)[(M+1):nrow(count_table)],]
FEAST_output = FEAST_output[,c(rownames(count_table)[1:M],'Unknown')]

colnames(FEAST_output) = c('Source 1','Source 2', 'Source 3','Unknown')
rownames(FEAST_output) = c('Sink 1','Sink 2', 'Sink 3','Sink 4','Sink 5','Sink 6','Sink 7')

melt_dat = melt(as.matrix(FEAST_output))
melt_ground = data.frame(Var1=c('Ground truth','Ground truth','Ground truth','Ground truth'),
                         Var2 = c('Source 1','Source 2', 'Source 3','Unknown'),
                         value = c(1/3,1/3,1/3,0))

melt_dat = rbind(melt_ground,melt_dat)
g1 = ggplot(melt_dat, aes(x = Var1, Var2, fill= value))+
  geom_tile(color = "white",size=1.5)+
  scale_fill_gradient2(low = "#44AA99", high = "#CC6677", mid = "white", 
                         midpoint = 0.2, limit = c(0,0.4), space = "Lab")+
  geom_text(aes(x = Var1, Var2, label= round(value,3)), color = "black", size = 3) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1),axis.text.y = element_text(size = 10))+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none")


tsne_out <- Rtsne(rbind(steady_source,steady_sink),pca=FALSE,perplexity=2,theta=0.0) 
df_tsne_out<-as.data.frame(tsne_out$Y)
df_tsne_out$types=c('source','source','source','sink','sink','sink','sink','sink','sink','sink')

g2 = ggplot(df_tsne_out, aes(x = V1, V2,fill = factor(types),shape=factor(types)))+
  geom_point(size=3, stroke=1,color='white')+scale_shape_manual(values=c(21, 23))+theme_bw()+xlab('tSNE_1')+ylab('tSNE_2')+
  scale_fill_manual(values = c("#1E88E5",'#D81B60'))+
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.background = element_rect(fill = "white", size = 4, colour = "white"),
    legend.justification = c(0, 1),
    axis.ticks = element_line(colour = "black", size = 0.2),
    axis.text.x = element_text(color="black", size=10),
    axis.text.y = element_text(color="black", size=10)
  )

p1 = grid.arrange(g2,g1,nrow = 1)
ggsave(p1,file="../figures/glv_priority_sigma_1.pdf",width=8.5, height=2.8)

