








\subsubsection{Avaliando melhores parâmetros para os suavizadores}

\subsubsection{Leave One Out Cross Validation}





```{r,warning=FALSE}

library(splines)
library(ISLR)
set.seed(103159)
df <- data.frame(x = x,y = y)


#### Bin Smoother LOOCV

# cv.error.min <- c(7L)
# tr= 1:nrow(df)
# dftr= df[tr,]
# number_of_bins = seq(2,30)
# cv.errors = matrix( NA,nrow =  nrow(df), ncol = length(number_of_bins) )
# 
# for( i in number_of_bins){ # for each number of knots to test
#   
#   for( j in tr ){ # for each fold
#     x1 <- with(dados, cut(x,i))
#     df_cut <- data.frame(y = df$y, x = x1)
#     fit.bin <- glm(y ~ x, data=df_cut[tr!=j,])
#     bin.pred = predict( fit.bin, newdata=data.frame(x = df_cut[tr==j,c("x")]) )
#     cv.errors[j,(i-1)] = mean( ( df_cut$y[tr==j] - bin.pred )^2 ) }
#   
# }
# 
# cv.errors.mean = apply(cv.errors,2,mean)
# cv.errors.stderr = apply(cv.errors,2,sd)/sqrt(nrow(dftr)-1)
# 
# min.cv.index = which.min( cv.errors.mean )+1
# cv1 = cv.errors.mean[min.cv.index-1]
# one_se_up_value = ( cv.errors.mean+cv.errors.stderr )[min.cv.index] 
# 
# # Set up the x-y limits for plotting:
# min_lim=min( one_se_up_value, cv.errors.mean, cv.errors.mean-cv.errors.stderr, cv.errors.mean+cv.errors.stderr ) * 0.9
# max_lim=max( one_se_up_value, cv.errors.mean, cv.errors.mean-cv.errors.stderr, cv.errors.mean+cv.errors.stderr ) * 1.1
# 
# df <- data.frame(x = x, y = y)
# df1 <- data.frame(x = number_of_bins, y = cv.errors.mean,EQM = cv.errors.mean)
# 
# 
# 
# colnames(df1) <- c("x","y",
#                     paste("EQM"))
# 
# df1 <- as.tibble(df1) %>%
#   gather(key = "variable", value = "value",-x,-y)
# 
# p1.cv <-ggplot(df1,aes(x = x,y=y))+
#   geom_line()+
#   geom_point()+
#   labs(x = "Binwidth",y = "EQM")+
#   geom_vline(xintercept = min.cv.index,color ="red")+
#   annotate("text",x = min.cv.index-8,y = max(cv.errors.mean)*0.8,label=paste0("Largura:",min.cv.index,"\nEQM:",round(cv1,4))) +
#   axis.theme()
# par1 = min.cv.index




```


```{r,echo=FALSE}

### Loess LOOCV

tr= 1:nrow(df)
dftr= df[tr,]
number_of_bins = seq(0.1,0.35,0.01)
cv.errors = matrix( NA,nrow =  nrow(df), ncol = length(number_of_bins) )


for( i in 1:length(number_of_bins)){ # for each number of knots to test
  
  for( j in tr ){ # for each fold
    
    
    fit.loess <- loess(y ~ x,degree=1, span = number_of_bins[i], data=df[tr!=j,])
    loess.pred = predict( fit.loess, newdata=data.frame(x = df$x[tr==j]) )
    cv.errors[j,i] = mean( ( df$y[tr==j] - loess.pred )^2,na.rm = TRUE ) }
  
}

cv.errors.mean = apply(cv.errors,2,mean,na.rm = TRUE)
cv.errors.stderr = apply(cv.errors,2,sd)/sqrt(nrow(dftr)-1)

min.cv.index = which.min( cv.errors.mean )
cv4 = cv.errors.mean[min.cv.index]
one_se_up_value = ( cv.errors.mean+cv.errors.stderr )[min.cv.index] 

# Set up the x-y limits for plotting:
min_lim=min( one_se_up_value, cv.errors.mean, cv.errors.mean-cv.errors.stderr, cv.errors.mean+cv.errors.stderr ) * 0.9
max_lim=max( one_se_up_value, cv.errors.mean, cv.errors.mean-cv.errors.stderr, cv.errors.mean+cv.errors.stderr ) * 1.1

df <- data.frame(x = x, y = y)
df1 <- data.frame(x = number_of_bins, y = cv.errors.mean,EQM = cv.errors.mean)



colnames(df1) <- c("x","y",
                   paste("EQM"))

df1 <- as.tibble(df1) %>%
  gather(key = "variable", value = "value",-x,-y)

p4.cv <-ggplot(df1,aes(x = x,y=y))+
  geom_line()+
  geom_point()+
  labs(x = "Binwidth",y = "EQM") +
  geom_vline(xintercept = number_of_bins[min.cv.index],color ="red")+
  annotate("text",x = number_of_bins[min.cv.index]-0.05,y = max(cv.errors.mean)*0.99,label=paste0("Span:",number_of_bins[min.cv.index],"\nEQM:",round(cv4,4))) +
  axis.theme()

par4 = number_of_bins[min.cv.index]


```


```{r,echo=FALSE}

### Kernel (Nadaraya-Watson) LOOCV
library(locfit)
tr= 1:nrow(df)
dftr= df[tr,]
number_of_bins = seq(0.1,0.35,0.01)
cv.errors = matrix( NA,nrow =  nrow(df), ncol = length(number_of_bins) )


for( i in 1:length(number_of_bins)){ # for each number of knots to test
  
  for( j in tr ){ # for each fold
    
    fit.kernel <- locfit(y ~ x,deg=1, alpha = number_of_bins[i],kern="gauss", data=df[tr!=j,])
    kernel.pred = predict( fit.kernel, newdata=data.frame(x = df$x[tr==j]) )
    cv.errors[j,i] = mean( ( df$y[tr==j] - kernel.pred )^2,na.rm = TRUE ) }
  
}

cv.errors.mean = apply(cv.errors,2,mean,na.rm = TRUE)
cv.errors.stderr = apply(cv.errors,2,sd)/sqrt(nrow(dftr)-1)

min.cv.index = which.min( cv.errors.mean )
cv5 = cv.errors.mean[min.cv.index]
one_se_up_value = ( cv.errors.mean+cv.errors.stderr )[min.cv.index] 

# Set up the x-y limits for plotting:
min_lim=min( one_se_up_value, cv.errors.mean, cv.errors.mean-cv.errors.stderr, cv.errors.mean+cv.errors.stderr ) * 0.9
max_lim=max( one_se_up_value, cv.errors.mean, cv.errors.mean-cv.errors.stderr, cv.errors.mean+cv.errors.stderr ) * 1.1

df <- data.frame(x = x, y = y)
df1 <- data.frame(x = number_of_bins, y = cv.errors.mean,EQM = cv.errors.mean)



colnames(df1) <- c("x","y",
                   paste("EQM"))

df1 <- as.tibble(df1) %>%
  gather(key = "variable", value = "value",-x,-y)

p5.cv <-ggplot(df1,aes(x = x,y=y))+
  geom_line()+
  geom_point()+
  labs(x = "Binwidth",y = "EQM") +
  geom_vline(xintercept = number_of_bins[min.cv.index],color ="red")+
  annotate("text",x = number_of_bins[min.cv.index]-0.05,y = max(cv.errors.mean)*.99,label=paste0("Span:",number_of_bins[min.cv.index],"\nEQM:",round(cv5,4))) +
  axis.theme()

par5 = number_of_bins[min.cv.index]

```


```{r,echo=FALSE}
library(rms)
### Kernel (Nadaraya-Watson) LOOCV
library(locfit)
tr= 1:nrow(df)
dftr= df[tr,]
number_of_bins = seq(2,20)
cv.errors = matrix( NA,nrow =  nrow(df), ncol = length(number_of_bins) )


for( i in number_of_bins){ # for each number of knots to test
  
  for( j in tr ){ # for each fold
    
    
    
    p              <-  seq(1,i-1,1)/i
    knots          <-  quantile(dados$x  , p = p)
    fit.spg1       <-  ols(as.formula(glue::glue("y ~ lsp(x, knots)")), data=df[tr!=j,])
    spg1.pred    = predict( fit.spg1, newdata=data.frame(x = df$x[tr==j]) )
    cv.errors[j,(i-1)] = mean( ( df$y[tr==j] - spg1.pred )^2,na.rm = TRUE ) }
  
}

cv.errors.mean = apply(cv.errors,2,mean,na.rm = TRUE)
cv.errors.stderr = apply(cv.errors,2,sd)/sqrt(nrow(dftr)-1)

min.cv.index = which.min( cv.errors.mean )-1
cv6 = cv.errors.mean[min.cv.index]
one_se_up_value = ( cv.errors.mean+cv.errors.stderr )[min.cv.index] 

# Set up the x-y limits for plotting:
min_lim=min( one_se_up_value, cv.errors.mean, cv.errors.mean-cv.errors.stderr, cv.errors.mean+cv.errors.stderr ) * 0.9
max_lim=max( one_se_up_value, cv.errors.mean, cv.errors.mean-cv.errors.stderr, cv.errors.mean+cv.errors.stderr ) * 1.1

df <- data.frame(x = x, y = y)
df1 <- data.frame(x = number_of_bins, y = cv.errors.mean,EQM = cv.errors.mean)



colnames(df1) <- c("x","y",
                   paste("EQM"))

df1 <- as.tibble(df1) %>%
  gather(key = "variable", value = "value",-x,-y)

p6.cv <-ggplot(df1,aes(x = x,y=y))+
  geom_line()+
  geom_point()+
  labs(x = "Binwidth",y = "EQM") +
  geom_vline(xintercept = number_of_bins[min.cv.index],color ="red")+
  annotate("text",x = number_of_bins[min.cv.index]+5,y = max(cv.errors.mean)*0.8,label=paste0("Span:",number_of_bins[min.cv.index]-1,"\nEQM:",round(cv6,4))) +
  axis.theme()

par6 = number_of_bins[min.cv.index]-1

```

```{r,echo=FALSE}
library(rms)
library(locfit)
### Kernel (Nadaraya-Watson) LOOCV

tr= 1:nrow(df)
dftr= df[tr,]
number_of_bins = seq(2,20)
cv.errors = matrix( NA,nrow =  nrow(df), ncol = length(number_of_bins) )


for( i in number_of_bins){ # for each number of knots to test
  
  for( j in tr ){ # for each fold
    
    
    p              <-  seq(1,i-1,1)/i
    knots          <-  quantile(dados$x  , p = p)
    fit.spg3       <- lm(y ~ bs(x, knots = knots), data=df[tr!=j,] )
    spg3.pred      = predict( fit.spg3, newdata=data.frame(x = df$x[tr==j]) )
    cv.errors[j,(i-1)] = mean( ( df$y[tr==j] - spg3.pred )^2,na.rm = TRUE ) }
  
}

cv.errors.mean = apply(cv.errors,2,mean,na.rm = TRUE)
cv.errors.stderr = apply(cv.errors,2,sd)/sqrt(nrow(dftr)-1)

min.cv.index = which.min( cv.errors.mean ) 
cv7 = cv.errors.mean[min.cv.index]
one_se_up_value = ( cv.errors.mean+cv.errors.stderr )[min.cv.index] 

# Set up the x-y limits for plotting:
min_lim=min( one_se_up_value, cv.errors.mean, cv.errors.mean-cv.errors.stderr, cv.errors.mean+cv.errors.stderr ) * 0.9
max_lim=max( one_se_up_value, cv.errors.mean, cv.errors.mean-cv.errors.stderr, cv.errors.mean+cv.errors.stderr ) * 1.1

df <- data.frame(x = x, y = y)
df1 <- data.frame(x = number_of_bins, y = cv.errors.mean,EQM = cv.errors.mean)



colnames(df1) <- c("x","y",
                   paste("EQM"))

df1 <- as.tibble(df1) %>%
  gather(key = "variable", value = "value",-x,-y)

p7.cv <-ggplot(df1,aes(x = x,y=y))+
  geom_line()+
  geom_point()+
  labs(x = "Binwidth",y = "EQM") +
  geom_vline(xintercept = number_of_bins[min.cv.index],color ="red")+
  annotate("text",x = number_of_bins[min.cv.index]+5,y = max(cv.errors.mean)*0.99,label=paste0("Span:",number_of_bins[min.cv.index]-1,"\nEQM:",round(cv7,4))) +
  axis.theme()

par7 = number_of_bins[min.cv.index]-1
```





```{r,fig.cap="\\label{smoothers_fit}Comparação entre diferentes ajustes para valores de parâmetros distintos.", fig.height=4.5}

partial_plots <- cowplot::plot_grid(p4.cv,p5.cv,p6.cv,p7.cv,ncol=2,nrow = 2,labels=paste0(LETTERS[1:4],c("-Loess","-Kernel","-Splines:Grau1","-Splines:Grau3")),vjust = -0.25,hjust = 0)
partial_plots

```



```{r,echo=FALSE,warning=FALSE}

library(Metrics)
df_metrics <- data.frame(Smoother = c("Loess", "Kernel","Splines Grau 1","Splines Grau 3"),
                         EQM      =  c(
                           cv4,
                           cv5,
                           cv6,
                           cv7))


kable_data(data = df_metrics,cap = "\\label{tab:tab_eqm_cenario1}Erro Quadrático Médio para os suavizadores Loess, Kernel e Spline Cúbico",foot = NULL)

```


```{r,  echo=F}
# bin  -  span = 11
# fit1  <- glm(y ~ cut(x, par1), data=dados)
# 
# 
# 
# # RM  -  span = 10
# fit2 <- ma(x = dados$y,order = 10)
# 
# # RL  -  span = 0.1
# fit3 <- with(dados, running.line(x = x, y = y, f = 0.1))
# 

# loess -  span = 0.2
fit4 <- loess(y ~ x, degree=1, span = par4, data=dados)



# kernel - span = 6
fit5 <- it.kernel <- locfit(y ~ x,deg=1, alpha = par5,kern="gauss", data=dados)


#Spline grau 1


k = par6

p         <-  seq(1,k-1,1)/k
knots     <-  quantile(dados$x  , p = p)
fit6       <-  ols(as.formula(glue::glue("y ~ lsp(x, knots)")), data = df)

# spline cubico
require(splines)
k = par7

p         <-  seq(1,k-1,1)/k
knots     <-  quantile(dados$x  , p = p)
fit7 <- lm(y ~ bs(x, knots = knots),data = dados )


spans = c(par4,par5,par6,par7)
df = cbind(x= dados$x, y = dados$y, Loess = fit4$fitted,Kernel = fitted.values(fit5),`Spline Grau 1`= predict(fit6),`Spline Cubico`= fit7$fitted)
df1 = df %>% as.data.frame
colnames(df) = c("x", "y",
                 paste0("A4\n", "Loess:", par4),
                 paste0("A5\n", "Kernel:", par5),
                 paste0("A6\n", "Spline Grau 1:", par6),
                 paste0("A7\n", "Spline Grau 3:", par7))

df = as_tibble(df) %>%
  gather(key = "variable", value = "value", -x, -y )

```

```{r echo=F,fig.height=3.5, fig.cap="\\label{fig:smoothers_fit_bestloocv}Comparação dos ajustes entre os métodos"}
plot.mult.curves(df = dados,df_fit = df,title = NULL,labelx = "X",labely = "Y",line.s = 1.05,alpha.o = .99)
```







\clearpage

\subsubsection{Cenário 2}

```{r , include=FALSE}
rm(list=ls())
source(file = "funcoes.R",encoding = "UTF-8")
library(tidyverse)
library(binsmooth)
library(knitr)
library(kableExtra)
library(additive.models)

knitr::opts_chunk$set(echo = FALSE,warning= FALSE, message= FALSE,
                      out.width = "100%",fig.align = "center",size ="large",fig.height = 3)

library(additive.models)

n <- 1e3
set.seed(103159)
n           <- 50
x           <- seq(0.1,2,length.out = 250)
norms       <- rnorm(length(x),0,0.1)
dens_gamma  <- dgamma(x = x,shape = 6,rate = 10)
y           <- dens_gamma + norms

dados       <- data.frame(x=x,y=y,variable = "Gamma(6,10)",value = dens_gamma)




```


\hspace{1.25cm} Para este cenário, foram geradas 201 observações, sendo x uma sequência de 0 à 2 com intervalos de 0.01. Ainda, temos que $y = f(x) + e$, com $f(x) \sim Gamma(6,10)$ e $e ~ N(0,0.10)$. O gráfico de dispersão para estes dados pode ser verificado na Figura 7, onde podemos observar o seu comportamento.


```{r,fig.cap="Gráfico de dispersão dos dados gerados para o estudo de simulação", fig.height=2.5}

####### Aqui
#plot.curves(data = dados,x =  x, y = y,title = NULL)
df           <- as.tibble(dados) 
#gather(key = "variable", value = "value",-x,-y)
plot.mult.curves(df = df,title = NULL)
```






```{r , echo = FALSE,fig.cap="\\label{fig:cenario2_loess}Alguns ajustes utilizando a ténica LOESS"}


spans = c(0.05,.15, 0.3,.6,.8)


df        <- cbind(dados$x,dados$y)
for(s in spans){
  
  fit     <-  loess(y ~ x, degree=1, span = s, data=dados)$fitted
  df      <-  cbind(df,fit)
}

colnames(df) <- c("x","y",
                  paste0("A1\n", "S1:",spans[1]),
                  paste0("A2\n", "S2:",spans[2]),
                  paste0("A3\n", "S3:",spans[3]),
                  paste0("A4\n", "S4:",spans[4]),
                  paste0("A5\n", "S5:",spans[5]))

df <- as.tibble(df) %>%
  gather(key = "variable", value = "value",-x,-y)


p4 <- plot.mult.curves(df = dados,df_fit = df,title = NULL,labelx = "X",labely = "Y",line.s = 1.05,alpha.o = .99)
```





```{r , echo=F, fig.cap="\\label{fig:cenario2_kernel}Alguns ajustes utilizando os Suavizadores com Kernel"}


spans = c(0.05,.15, 0.3,.6,.8)


df        <- cbind(dados$x,dados$y)
for(s in spans){
  fit <- locfit(y ~ x,deg=1, alpha = s,kern="gauss", data=dados)
  #fit     <-  ksmooth(x, y, kernel = "normal", bandwidth = s)$y
  df      <-  cbind(df,fitted(fit))
}

colnames(df) <- c("x","y",
                  paste0("A1\n", "L1:",spans[1]),
                  paste0("A2\n", "L2:",spans[2]),
                  paste0("A3\n", "L3:",spans[3]),
                  paste0("A4\n", "L4:",spans[4]),
                  paste0("A5\n", "L5:",spans[5]))

df <- as.tibble(df) %>%
  gather(key = "variable", value = "value",-x,-y)
p5 <- plot.mult.curves(df = dados,df_fit = df,title = NULL,labelx = "X",labely = "Y",line.s = 1.05,alpha.o = .99)

```




```{r , echo=F,fig.cap="\\label{fig:cenario2_sp1}Alguns ajustes utilizando a técnica Spliness de Regressão de Grau 1.",warning=FALSE}
no <- function(x, no) ifelse(x < no, 0, x-no)  # funcao truncada, so pega os positivos

library(segmented)
#p_load(kirkegaard, rms)
library(rms)
library(pacman)
library(igraph)
library(zoo)
library(splines)


df        <- as.data.frame(cbind(dados$x,dados$y))
nos       = c(2,3,5,11,21) 
colnames(df)<-c("x","y")
for(k in nos){
  
  p         <-  seq(1,k-1,1)/k
  knots     <-  quantile(df$x  , p = p)
  
  
  fit       <-  ols(as.formula(glue::glue("y ~ lsp(x, knots)")), data = df)
  df        <-  cbind(df,predict(fit))
}

colnames(df) <- c("x","y",
                  paste("A1\n", "Nós:",nos[1]-1),
                  paste("A2\n", "Nós:",nos[2]-1),
                  paste("A3\n", "Nós:",nos[3]-1),
                  paste("A4\n", "Nós:",nos[4]-1),
                  paste("A5\n", "Nós:",nos[5]-1)
)

df <- as.tibble(df) %>%
  gather(key = "variable", value = "value",-x,-y)

p6 <- plot.mult.curves(df = dados,df_fit = df,title = NULL,labelx = "X",labely = "Y",line.s = 1.05,alpha.o = .99)

```




```{r, echo=F,fig.cap="\\label{fig:cenario2_sp3}Alguns ajustes utilizando a técnica Spliness de Regressão"}
library(igraph)
library(zoo)
library(splines)


df        <- as.data.frame(cbind(dados$x,dados$y))
nos       = c(2,3,5,11,21) 


for(k in nos){
  
  
  p         <-  seq(1,k-1,1)/k
  knots     <-  quantile(dados$x, p = p)
  fit       <-  additive.spline.cubic(x = dados$x, y = dados$y, k = k,knots = knots)$fitted.values
  df        <-  cbind(df,fit)
}

colnames(df) <- c("x","y",
                  paste("A1\n", "Nós:",nos[1]-1),
                  paste("A2\n", "Nós:",nos[2]-1),
                  paste("A3\n", "Nós:",nos[3]-1),
                  paste("A4\n", "Nós:",nos[4]-1),
                  paste("A5\n", "Nós:",nos[5]-1)
)

df <- as.tibble(df) %>%
  gather(key = "variable", value = "value",-x,-y)





p7 <- plot.mult.curves(df = dados,df_fit = df,title = NULL,labelx = "X",labely = "Y",line.s = 1.05,alpha.o = .99)

```







```{r,fig.cap="\\label{smoothers_fit}Comparação entre diferentes ajustes para valores de parâmetros distintos.", fig.height=4.5}

partial_plots <- cowplot::plot_grid(p4,p5,p6,p7,ncol=2,nrow = 2,labels=paste0(LETTERS[1:4],c("-Loess","-Kernel","-Splines:Grau1","-Splines:Grau3")),vjust = 1,hjust = 0)
partial_plots

```


Figura \ref{fig:smoothers_fit_comp1}



Tabela \ref{tab:tab_eqm_cenario2}



