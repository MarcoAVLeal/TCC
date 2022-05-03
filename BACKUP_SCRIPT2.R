\clearpage

\subsection{Aplicação 1}


```{r}
source(file = "funcoes.R",encoding = "UTF-8")
library(locfit)
library(splines)
library(segmented)
library(rms)
library(pacman)
library(igraph)
library(zoo)

dados        <-    read.table(file = "aquiferborehole.csv",
                              header = T,sep = ",",encoding = "UTF-8")

x1            <- dados[,4]
x2            <- dados[,5]
x3            <- dados[,6]
y             <- dados[,7]
dados         <- data.frame(x = x2,  y = y)
x = dados$x
y = dados$y
```




```{r,fig.cap="\\label{fig:dispersao_aplicacao1}Gráfico de dispersão dos dados gerados para o estudo de simulação", fig.height=2.25}

#plot.curves(x = x1,y = y)
plot.curves(x = dados$x,y = dados$y)
#plot.curves(x = x3,y = y)
```


```{r , echo=F}

#fit <- locfit(y ~ x,deg=1, alpha = 0.1,kern="gauss", data=dados)
spans = c(0.1,.2, 0.3,.6,.8)


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
p4 <- plot.mult.curves(df = dados,df_fit = df,title = NULL,labelx = "X",labely = "Y",line.s = 1.05,alpha.o = .99)

```


```{r , echo = FALSE,fig.cap="Alguns ajustes utilizando a ténica LOESS"}


spans = c(0.05,.15, 0.3,.6,.8)


df        <- as.data.frame(cbind(dados$x,dados$y))
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


p5 <- plot.mult.curves(df = dados,df_fit = df,title = NULL,labelx = "X",labely = "Y",line.s = 1.05,alpha.o = .99)
```









```{r , echo=F,fig.cap="Alguns ajustes utilizando a técnica Spliness de Regressão de Grau 1.",warning=FALSE}
no <- function(x, no) ifelse(x < no, 0, x-no)  # funcao truncada, so pega os positivos



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




```{r, echo=F,fig.cap="Alguns ajustes utilizando a técnica Spliness de Regressão"}



df        <- as.data.frame(cbind(dados$x,dados$y))
nos       = c(2,3,5,11,21) 


for(k in nos){
  
  
  p         <-  seq(1,k-1,1)/k
  knots     <-  quantile(dados$x, p = p)
  fit       <-  lm(y ~ bs(x, knots = knots),data = dados )$fitted
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







Pela Figura \ref{fig:smoothers_fit_aplicacao_1}


```{r,fig.cap="\\label{fig:smoothers_fit_aplicacao_1}Comparação entre diferentes ajustes para valores de parâmetros distintos, considerando os suavizadores (A) Kernel, (B) Loess, (C) Splines de Regressão Grau 1 e (D) Splines de Regressão Grau 3.", fig.height=4.0}

partial_plots <- cowplot::plot_grid(p4,p5,p6,p7,ncol=2,nrow = 2,labels=LETTERS[1:4],vjust = 1,hjust = 0)
partial_plots

```

\subsubsection{Leave One Out Cross Validation}


```{r,echo=FALSE}

### Loess LOOCV
df        <- cbind(x= dados$x,y = dados$y) %>% as.data.frame
tr= 1:nrow(df)
dftr= df[tr,]
number_of_bins   = seq(0.1,0.35,0.01)
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
par.sp1              = number_of_bins_sp[min.cv.index.sp1]  

cv.errors.mean.sp3   = apply(cv.errors.sp3,2,mean,na.rm = TRUE)
min.cv.index.sp3     = which.min( cv.errors.mean.sp3 )
cv.min.sp3           = cv.errors.mean.sp3[min.cv.index.sp3]
par.sp3              = number_of_bins_sp[min.cv.index.sp3]


cv.errors.mean.loess   = apply(cv.errors.loess,2,mean,na.rm = TRUE)
min.cv.index.loess     = which.min( cv.errors.mean.loess )
cv.min.loess           = cv.errors.mean.loess[min.cv.index.loess]
par.loess              = number_of_bins[min.cv.index.loess]
### Kernel (Nadaraya-Watson) LOOCV

cv.errors.mean.kernel  = apply(cv.errors.kernel,2,mean,na.rm = TRUE)
min.cv.index.kernel    = which.min( cv.errors.mean.kernel )
cv.min.kernel          = cv.errors.mean.kernel[min.cv.index.kernel]
par.kernel              = number_of_bins[min.cv.index.kernel]


mt <- c()
mt <- rbind(mt,c(cv.min.kernel,cv.min.loess,cv.min.sp1,cv.min.sp3,par.loess,par.kernel,par.sp1,par.sp3))





```




```{r,echo=FALSE}



df <- data.frame(x = x, y = y)
df1 <- data.frame(x = number_of_bins, y = cv.errors.mean.kernel,EQM = cv.errors.mean.kernel)



colnames(df1) <- c("x","y",
                   paste("EQM"))

df1 <- as.tibble(df1) %>%
  gather(key = "variable", value = "value",-x,-y)

p4.cv <-ggplot(df1,aes(x = x,y=y))+
  geom_line()+
  geom_point()+
  labs(x = "Binwidth",y = "EQM") +
  geom_vline(xintercept = par.kernel,color ="red")+
  annotate("text",x = number_of_bins[min.cv.index.kernel]-0.05,y = max(cv.errors.mean.kernel)*0.97,label=paste0("Span:",number_of_bins[min.cv.index.kernel],"\nEQM:",round(cv.min.kernel,4))) +
  axis.theme()

par4 = number_of_bins[min.cv.index.kernel]

```


```{r,echo=FALSE}

### Kernel (Nadaraya-Watson) LOOCV
df <- data.frame(x = x, y = y)
df1 <- data.frame(x = number_of_bins, y = cv.errors.mean.loess,EQM = cv.errors.mean.loess)



colnames(df1) <- c("x","y",
                   paste("EQM"))

df1 <- as.tibble(df1) %>%
  gather(key = "variable", value = "value",-x,-y)

p5.cv <-ggplot(df1,aes(x = x,y=y))+
  geom_line()+
  geom_point()+
  labs(x = "Binwidth",y = "EQM") +
  geom_vline(xintercept = par.loess,color ="red")+
  annotate("text",x = number_of_bins[min.cv.index.loess]-0.05,y = max(cv.errors.mean.loess)*0.97,label=paste0("Span:",number_of_bins[min.cv.index.loess],"\nEQM:",round(cv.min.loess,4))) +
  axis.theme()

par5 = number_of_bins[min.cv.index.loess]
```


```{r,echo=FALSE}
df <- data.frame(x = dados$x, y = dados$y)
df1 <- data.frame(x = number_of_bins_sp, y = cv.errors.mean.sp1,EQM = cv.errors.mean.sp1)



colnames(df1) <- c("x","y",
                   paste("EQM"))

df1 <- as.tibble(df1) %>%
  gather(key = "variable", value = "value",-x,-y)

p6.cv <-ggplot(df1,aes(x = x,y=y))+
  geom_line()+
  geom_point()+
  labs(x = "Binwidth",y = "EQM") +
  geom_vline(xintercept = number_of_bins_sp[min.cv.index.sp1],color ="red")+
  annotate("text",x = number_of_bins_sp[min.cv.index.sp1]-3,y = max(cv.errors.mean.sp1)*0.90,label=paste0("Nº Nós:",par.sp1,"\nEQM:",round(cv.min.sp1,4))) +
  axis.theme()

par6 = number_of_bins[min.cv.index.sp1]
```

```{r,echo=FALSE}
df <- data.frame(x = x, y = y)
df1 <- data.frame(x = number_of_bins_sp, y = cv.errors.mean.sp3,EQM = cv.errors.mean.sp3)



colnames(df1) <- c("x","y",
                   paste("EQM"))

df1 <- as.tibble(df1) %>%
  gather(key = "variable", value = "value",-x,-y)

p7.cv <-ggplot(df1,aes(x = x,y=y))+
  geom_line()+
  geom_point()+
  labs(x = "Binwidth",y = "EQM") +
  geom_vline(xintercept = number_of_bins_sp[min.cv.index.sp3],color ="red")+
  annotate("text",x = number_of_bins_sp[min.cv.index.sp3]+5,y = max(cv.errors.mean.sp3)*0.90,label=paste0("Nº Nós:",par.sp3,"\nEQM:",round(cv.min.sp3,4))) +
  axis.theme()

par7 = number_of_bins[min.cv.index.sp3]

```


```{r,fig.cap="\\label{fig:smoothers_cv_aplicacao_1}Comparação entre diferentes ajustes para valores de parâmetros distintos, considerando os suavizadores (A) Kernel, (B) Loess, (C) Splines de Regressão Grau 1 e (D) Splines de Regressão Grau 3.", fig.height=4.0}

partial_plots <- cowplot::plot_grid(p4.cv,p5.cv,p6.cv,p7.cv,ncol=2,nrow = 2,labels=LETTERS[1:4])
partial_plots

```




```{r,echo=FALSE,warning=FALSE}

library(Metrics)
df_metrics <- data.frame(Smoother = c("Kernel", "Loess","Splines Grau 1","Splines Grau 3"),
                         EQM      =  c(
                           cv.min.kernel,
                           cv.min.loess,
                           cv.min.sp1,
                           cv.min.sp3))


kable_data(data = df_metrics,cap = "\\label{tab:tab_eqm_aplicao1}Erro Quadrático Médio para os suavizadores Loess, Kernel e Spline Cúbico",foot = NULL)

```


```{r,  echo=F}

fit4 <- loess(y ~ x, degree=1, span = par.loess, data=dados)



# kernel - span = 6
fit5 <-  locfit(y ~ x,deg=1, alpha = par.kernel,kern="gauss", data=dados)


#Spline grau 1


k = par.sp1

p         <-  seq(1,k-1,1)/(k)
knots     <-  quantile(dados$x  , p = p)
fit6       <-  ols(as.formula(glue::glue("y ~ lsp(x, knots)")), data = dados)

# spline cubico
require(splines)

k = par.sp3 
p         <-  seq(1,k-1,1)/(k)
knots     <-  quantile(dados$x  , p = p)
fit7 <- lm(y ~ bs(x, knots = knots),data = dados )


spans = c(par4,par5,par6,par7)
df = cbind(x= dados$x, y = dados$y, Loess = fit4$fitted,Kernel = fitted.values(fit5),`Spline Grau 1`= predict(fit6),`Spline Cubico`= fit7$fitted)
df1 = df %>% as.data.frame
colnames(df) = c("x", "y",
                 paste0("A4\n", "Loess|Span: ", par.loess),
                 paste0("A5\n", "Kernel|Width: ", par.kernel),
                 paste0("A6\n", "Spline Grau 1|Nº Nós: ", par.sp1),
                 paste0("A7\n", "Spline Grau 3|Nº Nós: ", par.sp3))

df = as_tibble(df) %>%
  gather(key = "variable", value = "value", -x, -y )

```

```{r echo=F,fig.height=3.0, fig.cap="\\label{fig:smoothers_fit_aplicacao_1_bestloocv}Comparação dos ajustes entre os métodos"}
plot.mult.curves(df = dados,df_fit = df,title = NULL,labelx = "X",labely = "Y",line.s = 1.05,alpha.o = .99)
```


\clearpage

\subsection{Aplicação 2}


```{r}
source(file = "funcoes.R",encoding = "UTF-8")
library(locfit)
library(splines)
library(segmented)
library(rms)
library(pacman)
library(igraph)
library(zoo)

dados        <-    read.table(file = "aquiferborehole.csv",
                              header = T,sep = ",",encoding = "UTF-8")

x1            <- dados[,4]
x2            <- dados[,5]
x3            <- dados[,6]
y             <- dados[,7]
dados         <- data.frame(x = x1,  y = y)
x = dados$x
y = dados$y
```




```{r,fig.cap="\\label{fig:dispersao_aplicacao2}Gráfico de dispersão dos dados gerados para o estudo de simulação", fig.height=2.25}

plot.curves(x = x1,y = y)
#plot.curves(x = dados$x,y = dados$y)
#plot.curves(x = x3,y = y)
```


```{r , echo=F}

#fit <- locfit(y ~ x,deg=1, alpha = 0.1,kern="gauss", data=dados)
spans = c(0.1,.2, 0.3,.6,.8)


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
p4 <- plot.mult.curves(df = dados,df_fit = df,title = NULL,labelx = "X",labely = "Y",line.s = 1.05,alpha.o = .99)

```


```{r , echo = FALSE,fig.cap="Alguns ajustes utilizando a ténica LOESS"}


spans = c(0.05,.15, 0.3,.6,.8)


df        <- as.data.frame(cbind(dados$x,dados$y))
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


p5 <- plot.mult.curves(df = dados,df_fit = df,title = NULL,labelx = "X",labely = "Y",line.s = 1.05,alpha.o = .99)
```









```{r , echo=F,fig.cap="Alguns ajustes utilizando a técnica Spliness de Regressão de Grau 1.",warning=FALSE}
no <- function(x, no) ifelse(x < no, 0, x-no)  # funcao truncada, so pega os positivos



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




```{r, echo=F,fig.cap="Alguns ajustes utilizando a técnica Spliness de Regressão"}



df        <- as.data.frame(cbind(dados$x,dados$y))
nos       = c(2,3,5,11,21) 


for(k in nos){
  
  
  p         <-  seq(1,k-1,1)/k
  knots     <-  quantile(dados$x, p = p)
  fit       <-  lm(y ~ bs(x, knots = knots),data = dados )$fitted
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







Pela Figura \ref{fig:smoothers_fit_aplicacao_2}


```{r,fig.cap="\\label{fig:smoothers_fit_aplicacao_2}Comparação entre diferentes ajustes para valores de parâmetros distintos, considerando os suavizadores (A) Kernel, (B) Loess, (C) Splines de Regressão Grau 1 e (D) Splines de Regressão Grau 3.", fig.height=4.0}

partial_plots <- cowplot::plot_grid(p4,p5,p6,p7,ncol=2,nrow = 2,labels=LETTERS[1:4],vjust = 1,hjust = 0)
partial_plots

```

\subsubsection{Leave One Out Cross Validation}


```{r,echo=FALSE}

### Loess LOOCV
df        <- cbind(x= dados$x,y = dados$y) %>% as.data.frame
tr= 1:nrow(df)
dftr= df[tr,]
number_of_bins   = seq(0.1,0.25,0.01)
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
par.sp1              = number_of_bins_sp[min.cv.index.sp1]  

cv.errors.mean.sp3   = apply(cv.errors.sp3,2,mean,na.rm = TRUE)
min.cv.index.sp3     = which.min( cv.errors.mean.sp3 )
cv.min.sp3           = cv.errors.mean.sp3[min.cv.index.sp3]
par.sp3              = number_of_bins_sp[min.cv.index.sp3]


cv.errors.mean.loess   = apply(cv.errors.loess,2,mean,na.rm = TRUE)
min.cv.index.loess     = which.min( cv.errors.mean.loess )
cv.min.loess           = cv.errors.mean.loess[min.cv.index.loess]
par.loess              = number_of_bins[min.cv.index.loess]
### Kernel (Nadaraya-Watson) LOOCV

cv.errors.mean.kernel  = apply(cv.errors.kernel,2,mean,na.rm = TRUE)
min.cv.index.kernel    = which.min( cv.errors.mean.kernel )
cv.min.kernel          = cv.errors.mean.kernel[min.cv.index.kernel]
par.kernel              = number_of_bins[min.cv.index.kernel]


mt <- c()
mt <- rbind(mt,c(cv.min.kernel,cv.min.loess,cv.min.sp1,cv.min.sp3,par.loess,par.kernel,par.sp1,par.sp3))





```




```{r,echo=FALSE}



df <- data.frame(x = x, y = y)
df1 <- data.frame(x = number_of_bins, y = cv.errors.mean.kernel,EQM = cv.errors.mean.kernel)



colnames(df1) <- c("x","y",
                   paste("EQM"))

df1 <- as.tibble(df1) %>%
  gather(key = "variable", value = "value",-x,-y)

p4.cv <-ggplot(df1,aes(x = x,y=y))+
  geom_line()+
  geom_point()+
  labs(x = "Binwidth",y = "EQM") +
  geom_vline(xintercept = par.kernel,color ="red")+
  annotate("text",x = number_of_bins[min.cv.index.kernel]-0.05,y = max(cv.errors.mean.kernel)*0.97,label=paste0("Span:",number_of_bins[min.cv.index.kernel],"\nEQM:",round(cv.min.kernel,4))) +
  axis.theme()

par4 = number_of_bins[min.cv.index.kernel]

```


```{r,echo=FALSE}

### Kernel (Nadaraya-Watson) LOOCV
df <- data.frame(x = x, y = y)
df1 <- data.frame(x = number_of_bins, y = cv.errors.mean.loess,EQM = cv.errors.mean.loess)



colnames(df1) <- c("x","y",
                   paste("EQM"))

df1 <- as.tibble(df1) %>%
  gather(key = "variable", value = "value",-x,-y)

p5.cv <-ggplot(df1,aes(x = x,y=y))+
  geom_line()+
  geom_point()+
  labs(x = "Binwidth",y = "EQM") +
  geom_vline(xintercept = par.loess,color ="red")+
  annotate("text",x = number_of_bins[min.cv.index.loess]-0.05,y = max(cv.errors.mean.loess)*0.97,label=paste0("Span:",number_of_bins[min.cv.index.loess],"\nEQM:",round(cv.min.loess,4))) +
  axis.theme()

par5 = number_of_bins[min.cv.index.loess]
```


```{r,echo=FALSE}
df <- data.frame(x = dados$x, y = dados$y)
df1 <- data.frame(x = number_of_bins_sp, y = cv.errors.mean.sp1,EQM = cv.errors.mean.sp1)



colnames(df1) <- c("x","y",
                   paste("EQM"))

df1 <- as.tibble(df1) %>%
  gather(key = "variable", value = "value",-x,-y)

p6.cv <-ggplot(df1,aes(x = x,y=y))+
  geom_line()+
  geom_point()+
  labs(x = "Binwidth",y = "EQM") +
  geom_vline(xintercept = number_of_bins_sp[min.cv.index.sp1],color ="red")+
  annotate("text",x = number_of_bins_sp[min.cv.index.sp1]-3,y = max(cv.errors.mean.sp1)*0.90,label=paste0("Nº Nós:",par.sp1,"\nEQM:",round(cv.min.sp1,4))) +
  axis.theme()

par6 = number_of_bins[min.cv.index.sp1]
```

```{r,echo=FALSE}
df <- data.frame(x = x, y = y)
df1 <- data.frame(x = number_of_bins_sp, y = cv.errors.mean.sp3,EQM = cv.errors.mean.sp3)



colnames(df1) <- c("x","y",
                   paste("EQM"))

df1 <- as.tibble(df1) %>%
  gather(key = "variable", value = "value",-x,-y)

p7.cv <-ggplot(df1,aes(x = x,y=y))+
  geom_line()+
  geom_point()+
  labs(x = "Binwidth",y = "EQM") +
  geom_vline(xintercept = number_of_bins_sp[min.cv.index.sp3],color ="red")+
  annotate("text",x = number_of_bins_sp[min.cv.index.sp3]+5,y = max(cv.errors.mean.sp3)*0.90,label=paste0("Nº Nós:",par.sp3,"\nEQM:",round(cv.min.sp3,4))) +
  axis.theme()

par7 = number_of_bins[min.cv.index.sp3]

```


```{r,fig.cap="\\label{fig:smoothers_cv_aplicacao_2}Comparação entre diferentes ajustes para valores de parâmetros distintos, considerando os suavizadores (A) Kernel, (B) Loess, (C) Splines de Regressão Grau 1 e (D) Splines de Regressão Grau 3.", fig.height=4.0}

partial_plots <- cowplot::plot_grid(p4.cv,p5.cv,p6.cv,p7.cv,ncol=2,nrow = 2,labels=LETTERS[1:4])
partial_plots

```




```{r,echo=FALSE,warning=FALSE}

library(Metrics)
df_metrics <- data.frame(Smoother = c("Kernel", "Loess","Splines Grau 1","Splines Grau 3"),
                         EQM      =  c(
                           cv.min.kernel,
                           cv.min.loess,
                           cv.min.sp1,
                           cv.min.sp3))


kable_data(data = df_metrics,cap = "\\label{tab:tab_eqm_aplicacao2} Erro Quadrático Médio para os suavizadores Loess, Kernel e Spline Cúbico",foot = NULL)

```


```{r,  echo=F}

fit4 <- loess(y ~ x, degree=1, span = par.loess, data=dados)



# kernel - span = 6
fit5 <-  locfit(y ~ x,deg=1, alpha = par.kernel,kern="gauss", data=dados)


#Spline grau 1


k = par.sp1

p         <-  seq(1,k,1)/(k+1)
knots     <-  quantile(dados$x  , p = p)
fit6       <-  ols(as.formula(glue::glue("y ~ lsp(x, knots)")), data = dados)

# spline cubico
require(splines)

k = par.sp3 
p         <-  seq(1,k,1)/(k+1)
knots     <-  quantile(dados$x  , p = p)
fit7 <- lm(y ~ bs(x, knots = knots),data = dados )


spans = c(par4,par5,par6,par7)
df = cbind(x= dados$x, y = dados$y, Loess = fit4$fitted,Kernel = fitted.values(fit5),`Spline Grau 1`= predict(fit6),`Spline Cubico`= fit7$fitted)
df1 = df %>% as.data.frame
colnames(df) = c("x", "y",
                 paste0("A4\n", "Loess|Span: ", par.loess),
                 paste0("A5\n", "Kernel|Width: ", par.kernel),
                 paste0("A6\n", "Spline Grau 1|Nº Nós: ", par.sp1),
                 paste0("A7\n", "Spline Grau 3|Nº Nós: ", par.sp3))

df = as_tibble(df) %>%
  gather(key = "variable", value = "value", -x, -y )

```

```{r echo=F,fig.height=3.0, fig.cap="\\label{fig:smoothers_fit_aplicacao_1_bestloocv}Comparação dos ajustes entre os métodos"}
plot.mult.curves(df = dados,df_fit = df,title = NULL,labelx = "X",labely = "Y",line.s = 1.05,alpha.o = .99)
```
