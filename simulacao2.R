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
var               <- c(0.5,1,2)
seed              <- 103159
dados             <- list()
dados_x           <- list()
i                 <- 1
for( n in 1:length(amostra)){
  mt <- matrix(data = NA,nrow = amostra[n],ncol = qtd.amostras)
  mt_x <- matrix(data = NA,nrow = amostra[n],ncol = qtd.amostras)
  for( v in 1:length(var)){
    
    for (j in 1:qtd.amostras){
      
      
      x           <- seq(0.1,50,length.out = amostra[n])
      set.seed(103159+j)
      norms       <- rnorm(length(x),0,var[v])
      d_sen       <- 10 + (5*sin(pi*x/24))
      y           <-  d_sen + norms
      
      mt[,j]      <- y
      
      #dados       <- data.frame(x=x,y=y,variable = "Gamma(6,10)",value = dens_gamma)
      
    }
    dados[[ i ]] <- as.data.frame(mt)
    dados_x[[ i ]] <- x
    i                    <- i + 1  
    
  }
  
}



dados_final <- read.csv(file = "dados_sim_final_cen1.csv",header = TRUE)


dados_n150 <- dados_final[1:3000,]
amostra_150 <- cbind(dados[[1]],dados[[2]],dados[[3]])
dados_x150  <-  cbind(dados_x[[1]],dados_x[[2]],dados_x[[3]])




cv.errors.loess  = matrix( NA,nrow =  3000, ncol = 1 )
# cv.errors.loess2  = matrix( NA,nrow =  3000, ncol = 1 )
# cv.errors.loess3  = matrix( NA,nrow =  3000, ncol = 1 )
cv.errors.kernel = matrix( NA,nrow =  3000, ncol = 1 )
cv.errors.sp1 = matrix( NA,nrow =  3000, ncol = 1 )
cv.errors.sp3 = matrix( NA,nrow =  3000, ncol = 1 ) 
cv.errors.poly  = matrix( NA,nrow =  3000, ncol = 1 ) 
mt1 <-c()


library(Metrics)
  
for( k in 1:3000){

  
  i <- as.integer( ceiling(k / 1000))
  
  df               = cbind(x= dados_x150[,1], y = amostra_150[,k]) %>% as.data.frame
  
    cat(i,"\\",k,"\n")
    
    
   
    
        
        fit.loess           <- loess(y ~ x,degree=1, span = dados_n150$V5[k], data=df)
        loess.pred     <- predict( fit.loess, newdata=data.frame(x = df$x) )
        cv.errors.loess[k,1]   = mean( ( df$y - loess.pred )^2,na.rm = TRUE ) 
        # cv.errors.loess2[k,1]   = rmse(  df$y , loess.pred )
        # cv.errors.loess3[k,1]   = mse(  df$y , loess.pred )
        
        fit.kernel             <- locfit(y ~ x,deg=1, alpha = dados_n150$V6[k],kern="gauss", data=df)
        kernel.pred    <- predict( fit.kernel, newdata=data.frame(x = df$x) )
        #cv.errors.kernel[k,1]  <- mean( ( df$y - kernel.pred )^2,na.rm = TRUE )
        cv.errors.kernel[k,1]  <- mse( df$y,kernel.pred )


        p              <-  seq(1,(dados_n150$V7[k]),1)/(dados_n150$V7[k]+1)
        knots          <-  quantile(df$x, p = p)
        fit.spg1       <-  ols(as.formula(glue::glue("y ~ lsp(x, knots)")), data=df)
        spg1.pred    <- fitted.values( fit.spg1 )
        cv.errors.sp1[k,1] <- mean( ( df$y - spg1.pred )^2,na.rm = TRUE )

        p              <-  seq(1,(dados_n150$V8[k]),1)/(dados_n150$V8[k]+1)
        knots          <-  quantile(df$x  , p = p)
        fit.spg3       <- lm(y ~ bs(x, knots = knots), data=df )
        spg3.pred      <- predict( fit.spg3, newdata=data.frame(x = df$x) )
        #cv.errors.sp3[k,1] <- mean( ( df$y - spg3.pred )^2,na.rm = TRUE )
        cv.errors.sp3[k,1]  <- mse( df$y,spg3.pred )
        
        
        fit.poly <- lm(y ~ poly(x = x,degree = 3),data = df)
        poly.pred      <- predict( fit.poly, newdata=data.frame(x = df$x) )
        cv.errors.poly[k,1] <- mean( ( df$y - poly.pred )^2,na.rm = TRUE )
        
      }
    
  
    
    mt1 <- cbind(cv.errors.kernel,cv.errors.loess,cv.errors.sp1,cv.errors.sp3,cv.errors.poly)
    

 
    dados_n150 <- dados_final[3001:6000,]
    amostra_150 <- cbind(dados[[4]],dados[[5]],dados[[6]])
    dados_x150  <-  cbind(dados_x[[4]],dados_x[[5]],dados_x[[6]])
    
    
    
    
    cv.errors.loess1  = matrix( NA,nrow =  3000, ncol = 1 )
    cv.errors.kernel1 = matrix( NA,nrow =  3000, ncol = 1 )
    cv.errors.sp11 = matrix( NA,nrow =  3000, ncol = 1 )
    cv.errors.sp31 = matrix( NA,nrow =  3000, ncol = 1 ) 
    cv.errors.poly1  = matrix( NA,nrow =  3000, ncol = 1 ) 
    mt2 <-c()
    
    
    library(Metrics)
    
    for( k in 1:3000){
      
      
      i <- as.integer( ceiling(k / 1000))
      
      df               = cbind(x= dados_x150[,1], y = amostra_150[,k]) %>% as.data.frame
      
      cat(i,"\\",k,"\n")
      
      
      
      
      
      fit.loess           <- loess(y ~ x,degree=1, span = dados_n150$V5[k], data=df)
      loess.pred     <- predict( fit.loess, newdata=data.frame(x = df$x) )
      cv.errors.loess1[k,1]   = mean( ( df$y - loess.pred )^2,na.rm = TRUE ) 
      # cv.errors.loess2[k,1]   = rmse(  df$y , loess.pred )
      # cv.errors.loess3[k,1]   = mse(  df$y , loess.pred )
      
      fit.kernel             <- locfit(y ~ x,deg=1, alpha = dados_n150$V6[k],kern="gauss", data=df)
      kernel.pred    <- predict( fit.kernel, newdata=data.frame(x = df$x) )
      #cv.errors.kernel[k,1]  <- mean( ( df$y - kernel.pred )^2,na.rm = TRUE )
      cv.errors.kernel1[k,1]  <- mse( df$y,kernel.pred )
      
      
      p              <-  seq(1,(dados_n150$V7[k]),1)/(dados_n150$V7[k]+1)
      knots          <-  quantile(df$x, p = p)
      fit.spg1       <-  ols(as.formula(glue::glue("y ~ lsp(x, knots)")), data=df)
      spg1.pred    <- fitted.values( fit.spg1 )
      cv.errors.sp11[k,1] <- mean( ( df$y - spg1.pred )^2,na.rm = TRUE )
      
      p              <-  seq(1,(dados_n150$V8[k]),1)/(dados_n150$V8[k]+1)
      knots          <-  quantile(df$x  , p = p)
      fit.spg3       <- lm(y ~ bs(x, knots = knots), data=df )
      spg3.pred      <- predict( fit.spg3, newdata=data.frame(x = df$x) )
      #cv.errors.sp3[k,1] <- mean( ( df$y - spg3.pred )^2,na.rm = TRUE )
      cv.errors.sp31[k,1]  <- mse( df$y,spg3.pred )
      
      
      fit.poly <- lm(y ~ poly(x = x,degree = 3),data = df)
      poly.pred      <- predict( fit.poly, newdata=data.frame(x = df$x) )
      cv.errors.poly1[k,1] <- mean( ( df$y - poly.pred )^2,na.rm = TRUE )
      
    }
    
    
    
    mt2 <- cbind(cv.errors.kernel1,cv.errors.loess1,cv.errors.sp11,cv.errors.sp31,cv.errors.poly1) 


    dados_n150 <- dados_final[6001:9000,]
    amostra_150 <- cbind(dados[[7]],dados[[8]],dados[[9]])
    dados_x150  <-  cbind(dados_x[[7]],dados_x[[8]],dados_x[[9]])
    
    #plot(amostra_150[,2500] ~ dados_x150[,1])
    
    
    cv.errors.loess2  = matrix( NA,nrow =  3000, ncol = 1 )
    # cv.errors.loess2  = matrix( NA,nrow =  3000, ncol = 1 )
    # cv.errors.loess3  = matrix( NA,nrow =  3000, ncol = 1 )
    cv.errors.kernel2 = matrix( NA,nrow =  3000, ncol = 1 )
    cv.errors.sp12 = matrix( NA,nrow =  3000, ncol = 1 )
    cv.errors.sp32 = matrix( NA,nrow =  3000, ncol = 1 ) 
    cv.errors.poly2  = matrix( NA,nrow =  3000, ncol = 1 ) 
    mt3 <-c()
    
    
    library(Metrics)
    
    for( k in 1:3000){
      
      
      i <- as.integer( ceiling(k / 1000))
      
      df               = cbind(x= dados_x150[,1], y = amostra_150[,k]) %>% as.data.frame
      
      cat(i,"\\",k,"\n")
      
      
      
      
      
      fit.loess           <- loess(y ~ x,degree=1, span = dados_n150$V5[k], data=df)
      loess.pred     <- predict( fit.loess, newdata=data.frame(x = df$x) )
      cv.errors.loess2[k,1]   = mean( ( df$y - loess.pred )^2,na.rm = TRUE ) 
      # cv.errors.loess2[k,1]   = rmse(  df$y , loess.pred )
      # cv.errors.loess3[k,1]   = mse(  df$y , loess.pred )
      
      fit.kernel             <- locfit(y ~ x,deg=1, alpha = dados_n150$V6[k],kern="gauss", data=df)
      kernel.pred    <- predict( fit.kernel, newdata=data.frame(x = df$x) )
      #cv.errors.kernel[k,1]  <- mean( ( df$y - kernel.pred )^2,na.rm = TRUE )
      cv.errors.kernel2[k,1]  <- mse( df$y,kernel.pred )
      
      
      p              <-  seq(1,(dados_n150$V7[k]),1)/(dados_n150$V7[k]+1)
      knots          <-  quantile(df$x, p = p)
      fit.spg1       <-  ols(as.formula(glue::glue("y ~ lsp(x, knots)")), data=df)
      spg1.pred    <- fitted.values( fit.spg1 )
      cv.errors.sp12[k,1] <- mean( ( df$y - spg1.pred )^2,na.rm = TRUE )
      
      p              <-  seq(1,(dados_n150$V8[k]),1)/(dados_n150$V8[k]+1)
      knots          <-  quantile(df$x  , p = p)
      fit.spg3       <- lm(y ~ bs(x, knots = knots), data=df )
      spg3.pred      <- predict( fit.spg3, newdata=data.frame(x = df$x) )
      #cv.errors.sp3[k,1] <- mean( ( df$y - spg3.pred )^2,na.rm = TRUE )
      cv.errors.sp32[k,1]  <- mse( df$y,spg3.pred )
      
      
      fit.poly <- lm(y ~ poly(x = x,degree = 3),data = df)
      poly.pred      <- predict( fit.poly, newdata=data.frame(x = df$x) )
      cv.errors.poly2[k,1] <- mean( ( df$y - poly.pred )^2,na.rm = TRUE )
      
    }
    
    
    
    mt3 <- cbind(cv.errors.kernel2,cv.errors.loess2,cv.errors.sp12,cv.errors.sp32,cv.errors.poly2)
    mt <- rbind(mt1,mt2,mt3)
    dados_final2 <- cbind(dados_final,mt)   
    dados_final2$EQM_MIN2 <- apply(X = dados_final2[,13:17],MARGIN = 1,FUN = which.min) 
    
    
    table_comp <- tapply(X = dados_final2$EQM_MIN2,INDEX = dados_final2$Tipo,FUN = table)
    
    table_mt <- matrix(data = 0,nrow = 9,ncol = 5)
    df_comp  <- as.data.frame(table_mt)
    colnames(df_comp) <- c(1,2,3,4,5)
    
    
    for(i in 1:9){
      ind = colnames(df_comp)%in% names(table_comp[[i]])
      pos = which(ind)
      df_comp[i,pos] <- table_comp[[i]]
      
      
    }
    
    
    
write.csv(x = dados_final2,file = "dados_sim_final_cen1_2.csv",row.names = FALSE)


############################################################

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




dados_final <- read.csv(file = "dados_sim_final_cen2.csv",header = TRUE)


dados_n150 <- dados_final[1:3000,]
amostra_150 <- cbind(dados[[1]],dados[[2]],dados[[3]])
dados_x150  <-  cbind(dados_x[[1]],dados_x[[2]],dados_x[[3]])




cv.errors.loess  = matrix( NA,nrow =  3000, ncol = 1 )
# cv.errors.loess2  = matrix( NA,nrow =  3000, ncol = 1 )
# cv.errors.loess3  = matrix( NA,nrow =  3000, ncol = 1 )
cv.errors.kernel = matrix( NA,nrow =  3000, ncol = 1 )
cv.errors.sp1 = matrix( NA,nrow =  3000, ncol = 1 )
cv.errors.sp3 = matrix( NA,nrow =  3000, ncol = 1 ) 
cv.errors.poly  = matrix( NA,nrow =  3000, ncol = 1 ) 
mt1 <-c()


library(Metrics)

for( k in 1:3000){
  
  
  i <- as.integer( ceiling(k / 1000))
  
  df               = cbind(x= dados_x150[,1], y = amostra_150[,k]) %>% as.data.frame
  
  cat(i,"\\",k,"\n")
  
  
  
  
  
  fit.loess           <- loess(y ~ x,degree=1, span = dados_n150$V5[k], data=df)
  loess.pred     <- predict( fit.loess, newdata=data.frame(x = df$x) )
  cv.errors.loess[k,1]   = mean( ( df$y - loess.pred )^2,na.rm = TRUE ) 
  # cv.errors.loess2[k,1]   = rmse(  df$y , loess.pred )
  # cv.errors.loess3[k,1]   = mse(  df$y , loess.pred )
  
  fit.kernel             <- locfit(y ~ x,deg=1, alpha = dados_n150$V6[k],kern="gauss", data=df)
  kernel.pred    <- predict( fit.kernel, newdata=data.frame(x = df$x) )
  #cv.errors.kernel[k,1]  <- mean( ( df$y - kernel.pred )^2,na.rm = TRUE )
  cv.errors.kernel[k,1]  <- mse( df$y,kernel.pred )
  
  
  p              <-  seq(1,(dados_n150$V7[k]),1)/(dados_n150$V7[k]+1)
  knots          <-  quantile(df$x, p = p)
  fit.spg1       <-  ols(as.formula(glue::glue("y ~ lsp(x, knots)")), data=df)
  spg1.pred    <- fitted.values( fit.spg1 )
  cv.errors.sp1[k,1] <- mean( ( df$y - spg1.pred )^2,na.rm = TRUE )
  
  p              <-  seq(1,(dados_n150$V8[k]),1)/(dados_n150$V8[k]+1)
  knots          <-  quantile(df$x  , p = p)
  fit.spg3       <- lm(y ~ bs(x, knots = knots), data=df )
  spg3.pred      <- predict( fit.spg3, newdata=data.frame(x = df$x) )
  #cv.errors.sp3[k,1] <- mean( ( df$y - spg3.pred )^2,na.rm = TRUE )
  cv.errors.sp3[k,1]  <- mse( df$y,spg3.pred )
  
  
  fit.poly <- lm(y ~ poly(x = x,degree = 3),data = df)
  poly.pred      <- predict( fit.poly, newdata=data.frame(x = df$x) )
  cv.errors.poly[k,1] <- mean( ( df$y - poly.pred )^2,na.rm = TRUE )
  
}



mt1 <- cbind(cv.errors.kernel,cv.errors.loess,cv.errors.sp1,cv.errors.sp3,cv.errors.poly)



dados_n150 <- dados_final[3001:6000,]
amostra_150 <- cbind(dados[[4]],dados[[5]],dados[[6]])
dados_x150  <-  cbind(dados_x[[4]],dados_x[[5]],dados_x[[6]])




cv.errors.loess1  = matrix( NA,nrow =  3000, ncol = 1 )
cv.errors.kernel1 = matrix( NA,nrow =  3000, ncol = 1 )
cv.errors.sp11 = matrix( NA,nrow =  3000, ncol = 1 )
cv.errors.sp31 = matrix( NA,nrow =  3000, ncol = 1 ) 
cv.errors.poly1  = matrix( NA,nrow =  3000, ncol = 1 ) 
mt2 <-c()


library(Metrics)

for( k in 1:3000){
  
  
  i <- as.integer( ceiling(k / 1000))
  
  df               = cbind(x= dados_x150[,1], y = amostra_150[,k]) %>% as.data.frame
  
  cat(i,"\\",k,"\n")
  
  
  
  
  
  fit.loess           <- loess(y ~ x,degree=1, span = dados_n150$V5[k], data=df)
  loess.pred     <- predict( fit.loess, newdata=data.frame(x = df$x) )
  cv.errors.loess1[k,1]   = mean( ( df$y - loess.pred )^2,na.rm = TRUE ) 
  # cv.errors.loess2[k,1]   = rmse(  df$y , loess.pred )
  # cv.errors.loess3[k,1]   = mse(  df$y , loess.pred )
  
  fit.kernel             <- locfit(y ~ x,deg=1, alpha = dados_n150$V6[k],kern="gauss", data=df)
  kernel.pred    <- predict( fit.kernel, newdata=data.frame(x = df$x) )
  #cv.errors.kernel[k,1]  <- mean( ( df$y - kernel.pred )^2,na.rm = TRUE )
  cv.errors.kernel1[k,1]  <- mse( df$y,kernel.pred )
  
  
  p              <-  seq(1,(dados_n150$V7[k]),1)/(dados_n150$V7[k]+1)
  knots          <-  quantile(df$x, p = p)
  fit.spg1       <-  ols(as.formula(glue::glue("y ~ lsp(x, knots)")), data=df)
  spg1.pred    <- fitted.values( fit.spg1 )
  cv.errors.sp11[k,1] <- mean( ( df$y - spg1.pred )^2,na.rm = TRUE )
  
  p              <-  seq(1,(dados_n150$V8[k]),1)/(dados_n150$V8[k]+1)
  knots          <-  quantile(df$x  , p = p)
  fit.spg3       <- lm(y ~ bs(x, knots = knots), data=df )
  spg3.pred      <- predict( fit.spg3, newdata=data.frame(x = df$x) )
  #cv.errors.sp3[k,1] <- mean( ( df$y - spg3.pred )^2,na.rm = TRUE )
  cv.errors.sp31[k,1]  <- mse( df$y,spg3.pred )
  
  
  fit.poly <- lm(y ~ poly(x = x,degree = 3),data = df)
  poly.pred      <- predict( fit.poly, newdata=data.frame(x = df$x) )
  cv.errors.poly1[k,1] <- mean( ( df$y - poly.pred )^2,na.rm = TRUE )
  
}



mt2 <- cbind(cv.errors.kernel1,cv.errors.loess1,cv.errors.sp11,cv.errors.sp31,cv.errors.poly1) 


dados_n150 <- dados_final[6001:9000,]
amostra_150 <- cbind(dados[[7]],dados[[8]],dados[[9]])
dados_x150  <-  cbind(dados_x[[7]],dados_x[[8]],dados_x[[9]])

#plot(amostra_150[,2500] ~ dados_x150[,1])


cv.errors.loess2  = matrix( NA,nrow =  3000, ncol = 1 )
# cv.errors.loess2  = matrix( NA,nrow =  3000, ncol = 1 )
# cv.errors.loess3  = matrix( NA,nrow =  3000, ncol = 1 )
cv.errors.kernel2 = matrix( NA,nrow =  3000, ncol = 1 )
cv.errors.sp12 = matrix( NA,nrow =  3000, ncol = 1 )
cv.errors.sp32 = matrix( NA,nrow =  3000, ncol = 1 ) 
cv.errors.poly2  = matrix( NA,nrow =  3000, ncol = 1 ) 
mt3 <-c()


library(Metrics)

for( k in 1:3000){
  
  
  i <- as.integer( ceiling(k / 1000))
  
  df               = cbind(x= dados_x150[,1], y = amostra_150[,k]) %>% as.data.frame
  
  cat(i,"\\",k,"\n")
  
  
  
  
  
  fit.loess           <- loess(y ~ x,degree=1, span = dados_n150$V5[k], data=df)
  loess.pred     <- predict( fit.loess, newdata=data.frame(x = df$x) )
  cv.errors.loess2[k,1]   = mean( ( df$y - loess.pred )^2,na.rm = TRUE ) 
  # cv.errors.loess2[k,1]   = rmse(  df$y , loess.pred )
  # cv.errors.loess3[k,1]   = mse(  df$y , loess.pred )
  
  fit.kernel             <- locfit(y ~ x,deg=1, alpha = dados_n150$V6[k],kern="gauss", data=df)
  kernel.pred    <- predict( fit.kernel, newdata=data.frame(x = df$x) )
  #cv.errors.kernel[k,1]  <- mean( ( df$y - kernel.pred )^2,na.rm = TRUE )
  cv.errors.kernel2[k,1]  <- mse( df$y,kernel.pred )
  
  
  p              <-  seq(1,(dados_n150$V7[k]),1)/(dados_n150$V7[k]+1)
  knots          <-  quantile(df$x, p = p)
  fit.spg1       <-  ols(as.formula(glue::glue("y ~ lsp(x, knots)")), data=df)
  spg1.pred    <- fitted.values( fit.spg1 )
  cv.errors.sp12[k,1] <- mean( ( df$y - spg1.pred )^2,na.rm = TRUE )
  
  p              <-  seq(1,(dados_n150$V8[k]),1)/(dados_n150$V8[k]+1)
  knots          <-  quantile(df$x  , p = p)
  fit.spg3       <- lm(y ~ bs(x, knots = knots), data=df )
  spg3.pred      <- predict( fit.spg3, newdata=data.frame(x = df$x) )
  #cv.errors.sp3[k,1] <- mean( ( df$y - spg3.pred )^2,na.rm = TRUE )
  cv.errors.sp32[k,1]  <- mse( df$y,spg3.pred )
  
  
  fit.poly <- lm(y ~ poly(x = x,degree = 3),data = df)
  poly.pred      <- predict( fit.poly, newdata=data.frame(x = df$x) )
  cv.errors.poly2[k,1] <- mean( ( df$y - poly.pred )^2,na.rm = TRUE )
  
}



mt3 <- cbind(cv.errors.kernel2,cv.errors.loess2,cv.errors.sp12,cv.errors.sp32,cv.errors.poly2)
mt <- rbind(mt1,mt2,mt3)
dados_final2 <- cbind(dados_final,mt)   
dados_final2$EQM_MIN2 <- apply(X = dados_final2[,13:17],MARGIN = 1,FUN = which.min) 


table_comp <- tapply(X = dados_final2$EQM_MIN2,INDEX = dados_final2$Tipo,FUN = table)

table_mt <- matrix(data = 0,nrow = 9,ncol = 5)
df_comp  <- as.data.frame(table_mt)
colnames(df_comp) <- c(1,2,3,4,5)


for(i in 1:9){
  ind = colnames(df_comp)%in% names(table_comp[[i]])
  pos = which(ind)
  df_comp[i,pos] <- table_comp[[i]]
  
  
}



write.csv(x = dados_final2,file = "dados_sim_final_cen2_2.csv")