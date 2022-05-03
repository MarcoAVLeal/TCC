rm(list=ls())
source(file = "funcoes.R",encoding = "UTF-8")
library(tidyverse)
library(binsmooth)
library(knitr)
library(kableExtra)
library(additive.models)
library(locfit)
library(rms)
library(splines)
library(foreach)
library(doParallel)
cat("\14")
qtd.amostras      <- 1000
amostra           <- c(150,250,350)
var               <- c(0.05,0.10,0.15)
seed              <- 103159
dados             <- list()
dados_x           <- list()
i                 <- 1
for( n in 1:length(amostra)){
  mt <- matrix(data = NA,nrow = amostra[n],ncol = qtd.amostras)
  mt_x <- matrix(data = NA,nrow = amostra[n],ncol = qtd.amostras)
  for( v in 1:length(var)){
    
    for (j in 1:qtd.amostras){
      
      
      x           <- seq(0.1,2,length.out = amostra[n])
      set.seed(103159+j)
      norms       <- rnorm(length(x),0,var[v])
      dens_gamma  <- dgamma(x = x,shape = 6,rate = 10)
      y           <- dens_gamma + norms
      
      mt[,j]      <- y
      
      #dados       <- data.frame(x=x,y=y,variable = "Gamma(6,10)",value = dens_gamma)
  
    }
    dados[[ i ]] <- as.data.frame(mt)
    dados_x[[ i ]] <- x
    i                    <- i + 1  
    
  }
  
}




### Loess LOOCV
i = 1
dados_eqm <- list()



for( n in 1:3){
  mt <-c()
  
  
  for (k in 1:qtd.amostras){
    cat(n,"\\",k,"\n")
    tr               = 1:length(dados[[n]][,k])
    df               = cbind(x= dados_x[[n]], y = dados[[n]][,k]) %>% as.data.frame
    dftr             = df[tr,]
    number_of_bins   = seq(0.05,0.25,0.01)
    number_of_bins_sp= 1:length(number_of_bins) 
    cv.errors.loess  = matrix( NA,nrow =  nrow(df), ncol = length(number_of_bins) )
    cv.errors.kernel = matrix( NA,nrow =  nrow(df), ncol = length(number_of_bins) )
    cv.errors.sp1 = matrix( NA,nrow =  nrow(df), ncol = length(number_of_bins_sp) )
    cv.errors.sp3 = matrix( NA,nrow =  nrow(df), ncol = length(number_of_bins_sp) ) 
    
    
    for( i in 1:length(number_of_bins)){ # for each number of knots to test
      
      for( j in tr ){ # for each fold
        
        
        fit.loess              <- loess(y ~ x,degree=1, span = number_of_bins[i], data=df[tr!=j,])
        loess.pred             = predict( fit.loess, newdata=data.frame(x = df$x[tr==j]) )
        cv.errors.loess[j,i]   = mean( ( df$y[tr==j] - loess.pred )^2,na.rm = TRUE ) 
        
        fit.kernel             <- locfit(y ~ x,deg=1, alpha = number_of_bins[i],kern="gauss", data=df[tr!=j,])
        kernel.pred            <- predict( fit.kernel, newdata=data.frame(x = df$x[tr==j]) )
        cv.errors.kernel[j,i]  <- mean( ( df$y[tr==j] - kernel.pred )^2,na.rm = TRUE ) 
        
        
        p              <-  seq(1,(i),1)/(i+1)
        knots          <-  quantile(df$x[tr!=j]  , p = p)
        fit.spg1       <-  ols(as.formula(glue::glue("y ~ lsp(x, knots)")), data=df[tr!=j,])
        spg1.pred    <- predict( fit.spg1, newdata=data.frame(x = df$x[tr==j]) )
        cv.errors.sp1[j,i] <- mean( ( df$y[tr==j] - spg1.pred )^2,na.rm = TRUE )
        
        #p              <-  seq(1,i-1,1)/i
        knots          <-  quantile(df$x[tr!=j]  , p = p)
        fit.spg3       <- lm(y ~ bs(x, knots = knots), data=df[tr!=j,] )
        spg3.pred      <- predict( fit.spg3, newdata=data.frame(x = df$x[tr==j]) )
        cv.errors.sp3[j,i] <- mean( ( df$y[tr==j] - spg3.pred )^2,na.rm = TRUE )
        
      }}
    
    cv.errors.mean.sp1   = apply(cv.errors.sp1,2,mean,na.rm = TRUE)
    min.cv.index.sp1     = which.min( cv.errors.mean.sp1 )
    cv.min.sp1           = cv.errors.mean.sp1[min.cv.index.sp1]
    par.sp1              = number_of_bins_sp[min.cv.index.sp1]  + 1
    
    cv.errors.mean.sp3   = apply(cv.errors.sp3,2,mean,na.rm = TRUE)
    min.cv.index.sp3     = which.min( cv.errors.mean.sp3 )
    cv.min.sp3           = cv.errors.mean.sp3[min.cv.index.sp3]
    par.sp3              = number_of_bins_sp[min.cv.index.sp3] + 1
    
    
    cv.errors.mean.loess   = apply(cv.errors.loess,2,mean,na.rm = TRUE)
    min.cv.index.loess     = which.min( cv.errors.mean.loess )
    cv.min.loess           = cv.errors.mean.loess[min.cv.index.loess]
    par.loess              = number_of_bins[min.cv.index.loess]
    ### Kernel (Nadaraya-Watson) LOOCV
    
    cv.errors.mean.kernel  = apply(cv.errors.kernel,2,mean,na.rm = TRUE)
    min.cv.index.kernel    = which.min( cv.errors.mean.kernel )
    cv.min.kernel          = cv.errors.mean.kernel[min.cv.index.kernel]
    par.kernel              = number_of_bins[min.cv.index.kernel]
    

 
    mt <- rbind(mt,c(cv.min.kernel,cv.min.loess,cv.min.sp1,cv.min.sp3,par.loess,par.kernel,par.sp1,par.sp3))
  }
  
  dados_eqm[[ n ]] <- as.data.frame(mt)
  cat("Fim",n)  
    
}


dados_copia = dados_eqm 

dados_copia[[1]]$n <- 150
dados_copia[[2]]$n <- 150
dados_copia[[3]]$n <- 150

dados_copia[[1]]$var <- 0.05
dados_copia[[2]]$var <- 0.1
dados_copia[[3]]$var  <- 0.15

dados_save <- rbind(dados_copia[[1]],dados_copia[[2]],dados_copia[[3]])

write.csv(x = dados_save,file = "dados_sim/dados_sim_1_3.csv",row.names = FALSE)





dados1 <- read.csv(file = "dados_sim/dados_sim_1_3.csv",header = TRUE)
dados2 <- read.csv(file = "dados_sim/dados_sim_4_5.csv",header = TRUE)
dados3 <- read.csv(file = "dados_sim/dados_sim_6_7.csv",header = TRUE)
dados4 <- read.csv(file = "dados_sim/dados_sim_8_9.csv",header = TRUE)

dados_final <- rbind(dados1,dados2,dados3,dados4)
dados_final$Tipo <- paste0(dados_final$n,dados_final$var)

dados_final$EQM_MIN <- apply(X = dados_final[,1:4],MARGIN = 1,FUN = which.min) 
write.csv(x = dados_final,file = "dados_sim_final.csv",row.names = FALSE)

#dados_final <- read.csv("dados_sim_final.csv",header = TRUE)


table_comp <- tapply(X = dados_final$EQM_MIN,INDEX = dados_final$Tipo,FUN = table)

table_mt <- matrix(data = 0,nrow = 9,ncol = 4)
df_comp  <- as.data.frame(table_mt)
colnames(df_comp) <- c(1,2,3,4)
table_comp$`1500.05`

for(i in 1:9){
  ind = colnames(df_comp)%in% names(table_comp[[i]])
  pos = which(ind)
  df_comp[i,pos] <- table_comp[[i]]
  
  
}



dados1 <- read.csv(file = "dados_sim_cen1//dados_sim_1_3.csv",header = TRUE)
dados2 <- read.csv(file = "dados_sim_cen1/dados_sim_4_5.csv",header = TRUE)
dados3 <- read.csv(file = "dados_sim_cen1/dados_sim_6_7.csv",header = TRUE)
dados4 <- read.csv(file = "dados_sim_cen1/dados_sim_8_9.csv",header = TRUE)

dados_final <- rbind(dados1,dados2,dados3,dados4)
dados_final$Tipo <- paste0(dados_final$n,dados_final$var)

dados_final$EQM_MIN <- apply(X = dados_final[,1:4],MARGIN = 1,FUN = which.min) 
write.csv(x = dados_final,file = "dados_sim_final_cen1.csv",row.names = FALSE)


table_comp <- tapply(X = dados_final$EQM_MIN,INDEX = dados_final$Tipo,FUN = table)

table_mt <- matrix(data = 0,nrow = 9,ncol = 4)
df_comp  <- as.data.frame(table_mt)
colnames(df_comp) <- c(1,2,3,4)
table_comp$`1500.05`

for(i in 1:9){
  ind = colnames(df_comp)%in% names(table_comp[[i]])
  pos = which(ind)
  df_comp[i,pos] <- table_comp[[i]]
  
  
}



################################




dados1 <- read.csv(file = "dados_sim_cen2//dados_sim1.csv",header = TRUE,sep = ";")
dados2 <- read.csv(file = "dados_sim_cen2//dados_sim2.csv",header = TRUE,sep = ";")
dados3 <- read.csv(file = "dados_sim_cen2//dados_sim3.csv",header = TRUE,sep = ";")
dados4 <- read.csv(file = "dados_sim_cen2//dados_sim4.csv",header = TRUE,sep = ";")
dados5 <- read.csv(file = "dados_sim_cen2//dados_sim5.csv",header = TRUE,sep = ";")
dados6 <- read.csv(file = "dados_sim_cen2//dados_sim6.csv",header = TRUE,sep = ";")
dados7 <- read.csv(file = "dados_sim_cen2//dados_sim7.csv",header = TRUE,sep = ";")
dados8 <- read.csv(file = "dados_sim_cen2//dados_sim8.csv",header = TRUE,sep = ";")
dados9 <- read.csv(file = "dados_sim_cen2//dados_sim9.csv",header = TRUE,sep = ";")

dados_final <- rbind(dados1,dados2,dados3,dados4,dados5,dados6,dados7,dados8,dados9)
dados_final$Tipo <- paste0(dados_final$n,dados_final$var)

dados_final$EQM_MIN <- apply(X = dados_final[,1:4],MARGIN = 1,FUN = which.min) 
write.csv(x = dados_final,file = "dados_sim_final_cen2.csv",row.names = FALSE)


table_comp <- tapply(X = dados_final$EQM_MIN,INDEX = dados_final$Tipo,FUN = table)

table_mt <- matrix(data = 0,nrow = 9,ncol = 4)
df_comp  <- as.data.frame(table_mt)
colnames(df_comp) <- c(1,2,3,4)
table_comp$`1500.05`

for(i in 1:9){
  ind = colnames(df_comp)%in% names(table_comp[[i]])
  pos = which(ind)
  df_comp[i,pos] <- table_comp[[i]]
  
  
}



data <- read.csv(file = "shotputt_powerclean.csv",header = TRUE)
plot.curves(x = data$V2,y = data$V1)


data <- read.table(textConnection(
  "
        .591E0         24.41E0  
       1.547E0         34.82E0  
       2.902E0         44.09E0  
       2.894E0         45.07E0  
       4.703E0         54.98E0  
       6.307E0         65.51E0  
       7.03E0          70.53E0  
       7.898E0         75.70E0  
       9.470E0         89.57E0  
       9.484E0         91.14E0  
      10.072E0         96.40E0  
      10.163E0         97.19E0  
      11.615E0        114.26E0  
      12.005E0        120.25E0  
      12.478E0        127.08E0  
      12.982E0        133.55E0  
      12.970E0        133.61E0  
      13.926E0        158.67E0  
      14.452E0        172.74E0  
      14.404E0        171.31E0  
      15.190E0        202.14E0  
      15.550E0        220.55E0  
      15.528E0        221.05E0  
      15.499E0        221.39E0  
      16.131E0        250.99E0  
      16.438E0        268.99E0  
      16.387E0        271.80E0  
      16.549E0        271.97E0  
      16.872E0        321.31E0  
      16.830E0        321.69E0  
      16.926E0        330.14E0  
      16.907E0        333.03E0  
      16.966E0        333.47E0  
      17.060E0        340.77E0  
      17.122E0        345.65E0  
      17.311E0        373.11E0  
      17.355E0        373.79E0  
      17.668E0        411.82E0  
      17.767E0        419.51E0  
      17.803E0        421.59E0  
      17.765E0        422.02E0  
      17.768E0        422.47E0  
      17.736E0        422.61E0  
      17.858E0        441.75E0  
      17.877E0        447.41E0  
      17.912E0        448.7E0   
      18.046E0        472.89E0  
      18.085E0        476.69E0  
      18.291E0        522.47E0  
      18.357E0        522.62E0  
      18.426E0        524.43E0  
      18.584E0        546.75E0  
      18.610E0        549.53E0  
      18.870E0        575.29E0  
      18.795E0        576.00E0  
      19.111E0        625.55E0  
        .367E0         20.15E0  
        .796E0         28.78E0  
       0.892E0         29.57E0  
       1.903E0         37.41E0  
       2.150E0         39.12E0  
       3.697E0         50.24E0  
       5.870E0         61.38E0  
       6.421E0         66.25E0  
       7.422E0         73.42E0  
       9.944E0         95.52E0  
      11.023E0        107.32E0  
      11.87E0         122.04E0  
      12.786E0        134.03E0  
      14.067E0        163.19E0  
      13.974E0        163.48E0  
      14.462E0        175.70E0  
      14.464E0        179.86E0  
      15.381E0        211.27E0  
      15.483E0        217.78E0  
      15.59E0         219.14E0  
      16.075E0        262.52E0  
      16.347E0        268.01E0  
      16.181E0        268.62E0  
      16.915E0        336.25E0  
      17.003E0        337.23E0  
      16.978E0        339.33E0  
      17.756E0        427.38E0  
      17.808E0        428.58E0  
      17.868E0        432.68E0  
      18.481E0        528.99E0  
      18.486E0        531.08E0  
      19.090E0        628.34E0  
      16.062E0        253.24E0  
      16.337E0        273.13E0  
      16.345E0        273.66E0  
      16.388E0        282.10E0  
      17.159E0        346.62E0  
      17.116E0        347.19E0  
      17.164E0        348.78E0  
      17.123E0        351.18E0  
      17.979E0        450.10E0  
      17.974E0        450.35E0  
      18.007E0        451.92E0  
      17.993E0        455.56E0  
      18.523E0        552.22E0  
      18.669E0        553.56E0  
      18.617E0        555.74E0  
      19.371E0        652.59E0  
      19.330E0        656.20E0  
       0.080E0         14.13E0  
       0.248E0         20.41E0  
       1.089E0         31.30E0  
       1.418E0         33.84E0  
       2.278E0         39.70E0  
       3.624E0         48.83E0  
       4.574E0         54.50E0  
       5.556E0         60.41E0  
       7.267E0         72.77E0  
       7.695E0         75.25E0  
       9.136E0         86.84E0  
       9.959E0         94.88E0  
       9.957E0         96.40E0  
      11.600E0        117.37E0  
      13.138E0        139.08E0  
      13.564E0        147.73E0  
      13.871E0        158.63E0  
      13.994E0        161.84E0  
      14.947E0        192.11E0  
      15.473E0        206.76E0  
      15.379E0        209.07E0  
      15.455E0        213.32E0  
      15.908E0        226.44E0  
      16.114E0        237.12E0  
      17.071E0        330.90E0  
      17.135E0        358.72E0  
      17.282E0        370.77E0  
      17.368E0        372.72E0  
      17.483E0        396.24E0  
      17.764E0        416.59E0  
      18.185E0        484.02E0  
      18.271E0        495.47E0  
      18.236E0        514.78E0  
      18.237E0        515.65E0  
      18.523E0        519.47E0  
      18.627E0        544.47E0  
      18.665E0        560.11E0  
      19.086E0        620.77E0  
       0.214E0         18.97E0  
       0.943E0         28.93E0  
       1.429E0         33.91E0  
       2.241E0         40.03E0  
       2.951E0         44.66E0  
       3.782E0         49.87E0  
       4.757E0         55.16E0  
       5.602E0         60.90E0  
       7.169E0         72.08E0  
       8.920E0         85.15E0  
      10.055E0         97.06E0  
      12.035E0        119.63E0  
      12.861E0        133.27E0  
      13.436E0        143.84E0  
      14.167E0        161.91E0  
      14.755E0        180.67E0  
      15.168E0        198.44E0  
      15.651E0        226.86E0  
      15.746E0        229.65E0  
      16.216E0        258.27E0  
      16.445E0        273.77E0  
      16.965E0        339.15E0  
      17.121E0        350.13E0  
      17.206E0        362.75E0  
      17.250E0        371.03E0  
      17.339E0        393.32E0  
      17.793E0        448.53E0  
      18.123E0        473.78E0  
      18.49E0         511.12E0  
      18.566E0        524.70E0  
      18.645E0        548.75E0  
      18.706E0        551.64E0  
      18.924E0        574.02E0  
      19.1E0          623.86E0  
       0.375E0         21.46E0  
       0.471E0         24.33E0  
       1.504E0         33.43E0  
       2.204E0         39.22E0  
       2.813E0         44.18E0  
       4.765E0         55.02E0  
       9.835E0         94.33E0  
      10.040E0         96.44E0  
      11.946E0        118.82E0  
      12.596E0        128.48E0  
      13.303E0        141.94E0  
      13.922E0        156.92E0  
      14.440E0        171.65E0  
      14.951E0        190.00E0  
      15.627E0        223.26E0  
      15.639E0        223.88E0  
      15.814E0        231.50E0  
      16.315E0        265.05E0  
      16.334E0        269.44E0  
      16.430E0        271.78E0  
      16.423E0        273.46E0  
      17.024E0        334.61E0  
      17.009E0        339.79E0  
      17.165E0        349.52E0  
      17.134E0        358.18E0  
      17.349E0        377.98E0  
      17.576E0        394.77E0  
      17.848E0        429.66E0  
      18.090E0        468.22E0  
      18.276E0        487.27E0  
      18.404E0        519.54E0  
      18.519E0        523.03E0  
      19.133E0        612.99E0  
      19.074E0        638.59E0  
      19.239E0        641.36E0  
      19.280E0        622.05E0  
      19.101E0        631.50E0  
      19.398E0        663.97E0  
      19.252E0        646.9E0   
      19.89E0         748.29E0  
      20.007E0        749.21E0  
      19.929E0        750.14E0  
      19.268E0        647.04E0  
      19.324E0        646.89E0  
      20.049E0        746.9E0   
      20.107E0        748.43E0  
      20.062E0        747.35E0  
      20.065E0        749.27E0  
      19.286E0        647.61E0  
      19.972E0        747.78E0  
      20.088E0        750.51E0  
      20.743E0        851.37E0  
      20.83E0         845.97E0  
      20.935E0        847.54E0  
      21.035E0        849.93E0  
      20.93E0         851.61E0  
      21.074E0        849.75E0  
      21.085E0        850.98E0  
      20.935E0        848.23E0
  
  
  
  
  "
  
  
  
))



# import numpy as np
# import matplotlib.pyplot as plt % matplotlib inline
# 
# x = np.arange(-5.0, 5.0, 0.1)
# 
# ## You can adjust the slope and intercept to verify the changes in the graph
# y = 1*(x**3) + 1*(x**2) + 1 * x + 3
# y_noise = 20 * np.random.normal(size = x.size)
# ydata = y + y_noise
# plt.plot(x, ydata,  'bo')
# plt.plot(x, y, 'r')
# plt.ylabel('Dependent Variable')
# plt.xlabel('Independent Variable')
# plt.show()

# registerDoParallel(cores = 4)
# 
# 
# 
# 
# 
# 
# #for( n in 1:9){
#   
#   mt <- c()
#   
#   mt <- foreach(k = 1:qtd.amostras,.combine = rbind)%dopar%{
#     library(tidyverse)
#     library(binsmooth)
#     library(knitr)
#     library(kableExtra)
#     library(additive.models)
#     library(locfit)
#     library(rms)
#     library(splines)
#     cat(n,"\\",k,"\n")
#     tr               = 1:length(dados[[1]][,k])
#     df               = cbind(x= dados_x[[1]], y = dados[[1]][,k]) %>% as.data.frame
#     dftr             = df[tr,]
#     number_of_bins   = seq(0.05,0.2,0.05)
#     cv.errors.loess  = matrix( NA,nrow =  nrow(df), ncol = length(number_of_bins) )
#     cv.errors.kernel = matrix( NA,nrow =  nrow(df), ncol = length(number_of_bins) )
#     
#     for( i in 1:length(number_of_bins)){ # for each number of knots to test
#       
#       for( j in tr ){ # for each fold
#         
#         
#         fit.loess              <- loess(y ~ x,degree=1, span = number_of_bins[i], data=df[tr!=j,])
#         loess.pred             = predict( fit.loess, newdata=data.frame(x = df$x[tr==j]) )
#         cv.errors.loess[j,i]   = mean( ( df$y[tr==j] - loess.pred )^2,na.rm = TRUE ) 
#         
#         fit.kernel             <- locfit(y ~ x,deg=1, alpha = number_of_bins[i],kern="gauss", data=df[tr!=j,])
#         kernel.pred            = predict( fit.kernel, newdata=data.frame(x = df$x[tr==j]) )
#         cv.errors.kernel[j,i]  = mean( ( df$y[tr==j] - kernel.pred )^2,na.rm = TRUE ) }
#       
#       
#     }
#     
#     
#     
#     
#     
#     number_of_bins_sp = seq(2,5)
#     cv.errors.sp1 = matrix( NA,nrow =  nrow(df), ncol = length(number_of_bins_sp) )
#     cv.errors.sp3 = matrix( NA,nrow =  nrow(df), ncol = length(number_of_bins_sp) )
#     
#     
#     for( i in number_of_bins_sp){ # for each number of knots to test
#       
#       for( j in tr ){ # for each fold
#         
#         
#         
#         p              <-  seq(1,i-1,1)/i
#         knots          <-  quantile(df$x  , p = p)
#         fit.spg1       <-  ols(as.formula(glue::glue("y ~ lsp(x, knots)")), data=df[tr!=j,])
#         spg1.pred    = predict( fit.spg1, newdata=data.frame(x = df$x[tr==j]) )
#         cv.errors.sp1[j,(i-1)] = mean( ( df$y[tr==j] - spg1.pred )^2,na.rm = TRUE ) 
#         
#         #p              <-  seq(1,i-1,1)/i
#         knots          <-  quantile(df$x  , p = p)
#         fit.spg3       <- lm(y ~ bs(x, knots = knots), data=df[tr!=j,] )
#         spg3.pred      = predict( fit.spg3, newdata=data.frame(x = df$x[tr==j]) )
#         cv.errors.sp3[j,(i-1)] = mean( ( df$y[tr==j] - spg3.pred )^2,na.rm = TRUE )
#         
#         
#       }
#       
#     }
#     
#     
#     
#     cv.errors.mean.sp1   = apply(cv.errors.sp1,2,mean,na.rm = TRUE)
#     min.cv.index.sp1     = which.min( cv.errors.mean.sp1 )
#     cv.min.sp1           = cv.errors.mean.sp1[min.cv.index.sp1]
#     par.sp1              = number_of_bins_sp[min.cv.index.sp1]-1
#     
#     cv.errors.mean.sp3   = apply(cv.errors.sp3,2,mean,na.rm = TRUE)
#     min.cv.index.sp3     = which.min( cv.errors.mean.sp3 )
#     cv.min.sp3           = cv.errors.mean.sp3[min.cv.index.sp3]
#     par.sp3              = number_of_bins_sp[min.cv.index.sp3]-1
#     
#     
#     cv.errors.mean.loess   = apply(cv.errors.loess,2,mean,na.rm = TRUE)
#     min.cv.index.loess     = which.min( cv.errors.mean.loess )
#     cv.min.loess           = cv.errors.mean.loess[min.cv.index.loess]
#     par.loess              = number_of_bins[min.cv.index.loess]
#     ### Kernel (Nadaraya-Watson) LOOCV
#     
#     cv.errors.mean.kernel  = apply(cv.errors.kernel,2,mean,na.rm = TRUE)
#     min.cv.index.kernel    = which.min( cv.errors.mean.kernel )
#     cv.min.kernel          = cv.errors.mean.kernel[min.cv.index.kernel]
#     par.kernel              = number_of_bins[min.cv.index.kernel]
#     
#     
#     
#     vec <- c(cv.min.kernel,cv.min.loess,cv.min.sp1,cv.min.sp3,par.loess,par.kernel,par.sp1,par.sp3)
#     vec
#   }
#   
#   dados_eqm[[ n ]] <- as.data.frame(mt)
#   cat("Fim",n)  
#   
# #}
# 
# 
# stopCluster()
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # dados       <- data.frame(x=x,y=mt[,1],variable = "Gamma(6,10)",value = dens_gamma)
# # df           <- as.tibble(dados) 
# # 
# # p1 <- plot.mult.curves(df = df,title = NULL)
# # 
# # 
# # dados       <- data.frame(x=x,y=mt[,2],variable = "Gamma(6,10)",value = dens_gamma)
# # df           <- as.tibble(dados) 
# # 
# # p2 <- plot.mult.curves(df = df,title = NULL)
# # 
# # dados       <- data.frame(x=x,y=mt[,3],variable = "Gamma(6,10)",value = dens_gamma)
# # df           <- as.tibble(dados) 
# # 
# # p3 <- plot.mult.curves(df = df,title = NULL)
# # 
# # 
# # partial_plots <- cowplot::plot_grid(p1,p2,p3,ncol=2,nrow = 2,labels=paste0(LETTERS[1:4]),vjust = 1,hjust = 0)
# # partial_plots==
