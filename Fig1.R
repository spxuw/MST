library("seqtime")
library("ggplot2")
setwd("/Users/xu-wenwang/Dropbox/Projects/FEAST/code")
library(jcolors)
library(reshape2)

################### parameter of glv ########################
set.seed(234)
N = 6
C = 0.5
sigma = 1
A = matrix(rnorm(N*N,mean=0,sd=sigma), N, N)
Connect = sample(0:1, N^2, prob = c(1-C, C), replace = TRUE)
A = A * matrix(Connect, N, N)
diag(A) = -1
b = runif(N)

A[1,2] = -2
A[1,3] = 1.5

## first source
set.seed(456)
y_0 = runif(N)
dis_sink = 0
while (dis_sink<0.1){
  collection = sample(1:N,N)
  y_1 = y_0
  y_2 = y_0
  y_3 = y_0
  y_1[collection[3:6]] = 0
  y_2[collection[c(1,2,5,6)]] = 0
  y_3[collection[1:4]] = 0

  y_1 = glv(N = N, A, b = b, y = y_1, tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
  y_2 = glv(N = N, A, b = b, y = y_2, tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
  y_3 = glv(N = N, A, b = b, y = y_3, tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)

  y_1 = y_1[,ncol(y_1)]/sum(y_1[,ncol(y_1)])
  y_2 = y_2[,ncol(y_2)]/sum(y_2[,ncol(y_2)])
  y_3 = y_3[,ncol(y_3)]/sum(y_3[,ncol(y_3)])
  
  # simutanonsly
  x_01_1 = glv(N = N, A, b = b, y = y_1*0.33+y_2*0.33+y_3*0.33, tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
  x_01_0 = y_1*0.33+y_2*0.33+y_3*0.33
  x_01_1 = x_01_1[,ncol(x_01_1)]/sum(x_01_1[,ncol(x_01_1)])
  
  # order 1-2-3
  x_11_1 = glv(N = N, A, b = b, y = y_1*0.33, tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
  x_11_2 = glv(N = N, A, b = b, y = x_11_1[,ncol(x_11_1)]+y_2*0.33, tstart = 100.1, tend = 200.1, tstep = 0.1, perturb = NULL)
  x_11_3 = glv(N = N, A, b = b, y = x_11_2[,ncol(x_11_2)]+y_3*0.33, tstart = 200.1, tend = 300.1, tstep = 0.1, perturb = NULL)
  x_11_3 = x_11_3[,ncol(x_11_3)]/sum(x_11_3[,ncol(x_11_3)])
  
  # order 1-3-2
  x_12_1 = glv(N = N, A, b = b, y = y_1*0.33, tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
  x_12_2 = glv(N = N, A, b = b, y = x_12_1[,ncol(x_12_1)]+y_3*0.33, tstart = 100.1, tend = 200.1, tstep = 0.1, perturb = NULL)
  x_12_3 = glv(N = N, A, b = b, y = x_12_2[,ncol(x_12_2)]+y_2*0.33, tstart = 200.1, tend = 300.1, tstep = 0.1, perturb = NULL)
  x_12_3 = x_12_3[,ncol(x_12_3)]/sum(x_12_3[,ncol(x_12_3)])
  
  # order 2-1-3
  x_21_1 = glv(N = N, A, b = b, y = y_2*0.33, tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
  x_21_2 = glv(N = N, A, b = b, y = x_21_1[,ncol(x_21_1)]+y_1*0.33, tstart = 100.1, tend = 200.1, tstep = 0.1, perturb = NULL)
  x_21_3 = glv(N = N, A, b = b, y = x_21_2[,ncol(x_21_2)]+y_3*0.33, tstart = 200.1, tend = 300.1, tstep = 0.1, perturb = NULL)
  x_21_3 = x_21_3[,ncol(x_21_3)]/sum(x_21_3[,ncol(x_21_3)])
  
  # order 2-3-1
  x_22_1 = glv(N = N, A, b = b, y = y_2*0.33, tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
  x_22_2 = glv(N = N, A, b = b, y = x_22_1[,ncol(x_22_1)]+y_3*0.33, tstart = 100.1, tend = 200.1, tstep = 0.1, perturb = NULL)
  x_22_3 = glv(N = N, A, b = b, y = x_22_2[,ncol(x_22_2)]+y_1*0.33, tstart = 200.1, tend = 300.1, tstep = 0.1, perturb = NULL)
  x_22_3 = x_22_3[,ncol(x_22_3)]/sum(x_22_3[,ncol(x_22_3)])
  
  # order 3-1-2
  x_31_1 = glv(N = N, A, b = b, y = y_3*0.33, tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
  x_31_2 = glv(N = N, A, b = b, y = x_31_1[,ncol(x_31_1)]+y_1*0.33, tstart = 100.1, tend = 200.1, tstep = 0.1, perturb = NULL)
  x_31_3 = glv(N = N, A, b = b, y = x_31_2[,ncol(x_31_2)]+y_2*0.33, tstart = 200.1, tend = 300.1, tstep = 0.1, perturb = NULL)
  x_31_3 = x_31_3[,ncol(x_31_3)]/sum(x_31_3[,ncol(x_31_3)])
  
  # order 3-2-1
  x_32_1 = glv(N = N, A, b = b, y = y_3*0.33, tstart = 0, tend = 100, tstep = 0.1, perturb = NULL)
  x_32_2 = glv(N = N, A, b = b, y = x_32_1[,ncol(x_32_1)]+y_2*0.33, tstart = 100.1, tend = 200.1, tstep = 0.1, perturb = NULL)
  x_32_3 = glv(N = N, A, b = b, y = x_32_2[,ncol(x_32_2)]+y_1*0.33, tstart = 200.1, tend = 300.1, tstep = 0.1, perturb = NULL)
  x_32_3 = x_32_3[,ncol(x_32_3)]/sum(x_32_3[,ncol(x_32_3)])
  
  steady_sink = cbind(x_11_3,x_12_3,x_21_3,x_22_3,x_31_3,x_32_3)
  dis_sink = mean(vegdist(t(steady_sink)))
}

steady_source = cbind(y_1/sum(y_1),y_2/sum(y_2),y_3/sum(y_3))
abundance = cbind(steady_source,steady_sink,x_01_0/sum(x_01_0),x_01_1)
colnames(abundance) = c('s1','s2','s3','s4','s5','s6','s7','s8','s9','s10','s11')

dat_abundance = melt(abundance)

g1 <- ggplot(data=dat_abundance, aes(x=Var2, y=abundance, fill=factor(Var1))) +
  geom_bar(stat="identity")+  scale_fill_brewer(palette = "Dark2")+ theme_classic()+ 
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.background = element_rect(fill = "white", size = 4, colour = "white"),
    legend.justification = c(0, 1),
    axis.ticks = element_line(colour = "black", size = 0.2),
    panel.grid.major = element_line(colour = "grey70", size = 0.2),
    panel.grid.minor = element_blank(),axis.text.x = element_text(color="black", size=10),
    axis.text.y = element_text(color="black", size=10)
  )
ggsave(g1, file="../figures/Demo_interaction_priority.pdf", width=4.5, height=3)

