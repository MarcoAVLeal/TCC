## Aplicação 2



\begin{columns}
\begin{column}{0.5\textwidth}

* Os dados para esta aplicação são referentes ao preço mediano de casas de Boston. Este conjunto de dados possui 506 observações e 14 variáveis. 

* A variável resposta $y$, a coluna *medv* que representa o valor mediano das casas ocupadas pelos proprietários em $1000s. 

* Para a variável preditora $x$, será considerado a coluna *lstat* que representa o percentual de status baixo da população. Tendo em conta o modelo aditivo da seguinte forma,

\begin{equation*}
y = \alpha + f(x) + \varepsilon,
\end{equation*}

onde os erros $\varepsilon$ são independentes, com $E(\varepsilon) = 0$ e $var(\varepsilon) = \sigma^2$. A $f(x)$ é uma função univariada arbitrária, que será suavizada pelos métodos vistos até o momento.

\end{column}

\begin{column}{0.5\textwidth}

```{r,echo=FALSE}
source(file = "funcoes.R",encoding = "UTF-8")
library(locfit)
library(splines)
library(segmented)
library(rms)
library(pacman)
library(igraph)
library(zoo)
data("Boston")


x1            <- Boston$lstat
y             <- Boston$medv
dados         <- data.frame(x = x1,  y = y)
x = dados$x
y = dados$y
```




```{r,fig.cap="\\label{fig:dispersao_aplicacao2}Gráfico de dispersão referente ao estudo de espansão térmica do cobre versus a temperatura em graus Kelvin.", fig.height=4}

#plot.curves(x = x1,y = y)
plot.curves(x = dados$x,y = dados$y,labelx = "% de Status Baixo da População",labely = "Preço Mediano de Casas",title = NULL)+
  axis.theme(lengend_text_size = 16,lengend_title_size = 16,textsize = 18,leg = FALSE,pos_leg = "none")
#plot.curves(x = x3,y = y)
```



\end{column}
\end{columns}

## Avaliando os $EQM_{loocv}$ e os $EQM_c$


\begin{columns}
\begin{column}{0.5\textwidth}


```{r,echo=FALSE}

###  LOOCV
df        <- cbind(x= dados$x,y = dados$y) %>% as.data.frame
tr= 1:nrow(df)
dftr= df[tr,]
number_of_bins   = seq(0.1,0.29,0.01)
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
  labs(x = "Bandwidth",y = "EQM") +
  geom_vline(xintercept = par.kernel,color ="red")+
  annotate("text",x = number_of_bins[min.cv.index.kernel]-0.06,y = max(cv.errors.mean.kernel)*0.985,label=paste0("Span:",number_of_bins[min.cv.index.kernel],"\nEQM:",round(cv.min.kernel,2))) +
  axis.theme(lengend_text_size = 16,lengend_title_size = 16,textsize = 18,leg = FALSE,pos_leg = "none")

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
  labs(x = "Span",y = "EQM") +
  geom_vline(xintercept = par.loess,color ="red")+
  annotate("text",x = number_of_bins[min.cv.index.loess]-0.06,y = max(cv.errors.mean.loess)*0.985,label=paste0("Span:",number_of_bins[min.cv.index.loess],"\nEQM:",round(cv.min.loess,2))) +
  axis.theme(lengend_text_size = 16,lengend_title_size = 16,textsize = 18,leg = FALSE,pos_leg = "none")

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
  labs(x = "Nº de Nós",y = "EQM") +
  geom_vline(xintercept = number_of_bins_sp[min.cv.index.sp1],color ="red")+
  annotate("text",x = number_of_bins_sp[min.cv.index.sp1]+6,y = max(cv.errors.mean.sp1)*0.90,label=paste0("Nº Nós:",par.sp1,"\nEQM:",round(cv.min.sp1,2))) +
  axis.theme(lengend_text_size = 16,lengend_title_size = 16,textsize = 18,leg = FALSE,pos_leg = "none")

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
  labs(x = "Nº de Nós",y = "EQM") +
  geom_vline(xintercept = number_of_bins_sp[min.cv.index.sp3],color ="red")+
  annotate("text",x = number_of_bins_sp[min.cv.index.sp3]+6,y = max(cv.errors.mean.sp3)*0.97,label=paste0("Nº Nós:",par.sp3,"\nEQM:",round(cv.min.sp3,2))) +
  axis.theme(lengend_text_size = 16,lengend_title_size = 16,textsize = 18,leg = FALSE,pos_leg = "none")

par7 = number_of_bins[min.cv.index.sp3]

```




```{r,fig.cap="\\label{fig:smoothers_best_par_aplicacao2} Erro quadrático médio versus parâmetro de suavização (Cenário 1) pós aplicação do Leave One Out Cross-Validation. (A) Kernel, (B) Loess, (C) Splines de Regressão Linear e (D) Splines de Regressão Cúbico.", fig.height=2.5,include=FALSE,echo=FALSE}

partial_plots <- cowplot::plot_grid(p4.cv,p5.cv,p6.cv,p7.cv,ncol=2,nrow = 2,labels=LETTERS[1:4])
partial_plots

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


require(splines)

k = par.sp3 
p         <-  seq(1,k-1,1)/(k)
knots     <-  quantile(dados$x  , p = p)
fit7 <- lm(y ~ bs(x, knots = knots),data = dados )

fit8 <- lm(y ~ poly(x = x,degree = 3),data = dados)


spans = c(par4,par5,par6,par7)
df = cbind(x= dados$x, y = dados$y,Kernel = fitted.values(fit5), Loess = fit4$fitted,`Spline Grau 1`= predict(fit6),`Spline Cubico`= fit7$fitted,`Pol. Cúbico` = fit8$fitted.values)
df1 = df %>% as.data.frame
colnames(df) = c("x", "y",
                 paste0("A1 ", "Kernel|Width: ", par.kernel),
                 paste0("A2 ", "Loess|Span: ", par.loess),
                 paste0("A3 ", "Spline Grau 1|Nº Nós: ", par.sp1),
                 paste0("A4 ", "Spline Grau 3|Nº Nós: ", par.sp3),
                 paste0("A5 ", "Polinómio Cúbico"))

df = as_tibble(df) %>%
  gather(key = "variable", value = "value", -x, -y )

```


```{r,fig.height=8, echo=FALSE, message=FALSE, warning=FALSE, include=TRUE, fig.cap="\\label{fig:smoothers_fit_bestloocv_aplicacao2} Comparação entre os ajustes, considerando parâmetros de suavização obtidos por meio da validação cruzada. (A) Kernel - Parâm. 0,21. (B) Loess - Parâm. 0,18. (C) Splines de Regressão Linear - Parâm. 5 (D) Splines de Regressão Cúbico - Parâm. 5  (E) Regressão Polinômial Cúbica."}


dados = data.frame(x,y)
true = fit4$fitted

df = cbind(dados$x, dados$y, true) %>% as.data.frame
colnames(df) = c("x", "y","Curva real")

df = as_tibble(df) %>%
  gather(key = "variable", value = "value", -x, -y )
p1 <- plot.mult.curves(df = dados,df_fit = df,title = NULL,labelx = "X",labely = "Y",legend.pos = "none")+
  axis.theme(lengend_text_size = 16,lengend_title_size = 16,textsize = 18,leg = FALSE,pos_leg = "none")

true = fitted(fit5)

df = cbind(dados$x, dados$y, true) %>% as.data.frame
colnames(df) = c("x", "y","Curva real")

df = as_tibble(df) %>%
  gather(key = "variable", value = "value", -x, -y )
p2 <- plot.mult.curves(df = dados,df_fit = df,title = NULL,labelx = "X",labely = "Y",legend.pos = "none")+
  axis.theme(lengend_text_size = 16,lengend_title_size = 16,textsize = 18,leg = FALSE,pos_leg = "none")



true = fit6$fitted.values

df = cbind(dados$x, dados$y, true) %>% as.data.frame
colnames(df) = c("x", "y","Curva real")


df = as_tibble(df) %>%
  gather(key = "variable", value = "value", -x, -y )
p3 <- plot.mult.curves(df = dados,df_fit = df,title = NULL,labelx = "X",labely = "Y",legend.pos = "none")+
  axis.theme(lengend_text_size = 16,lengend_title_size = 16,textsize = 18,leg = FALSE,pos_leg = "none")

true = fit7$fitted.values

df = cbind(dados$x, dados$y, true) %>% as.data.frame
colnames(df) = c("x", "y","Curva real")

df = as_tibble(df) %>%
  gather(key = "variable", value = "value", -x, -y )
p4 <- plot.mult.curves(df = dados,df_fit = df,title = NULL,labelx = "X",labely = "Y",legend.pos = "none")+
  axis.theme(lengend_text_size = 16,lengend_title_size = 16,textsize = 18,leg = FALSE,pos_leg = "none")


true = fit8$fitted.values

df = cbind(dados$x, dados$y, true) %>% as.data.frame
colnames(df) = c("x", "y","Curva real")

df = as_tibble(df) %>%
  gather(key = "variable", value = "value", -x, -y )
p5 <- plot.mult.curves(df = dados,df_fit = df,title = NULL,labelx = "X",labely = "Y",legend.pos = "none")+
  axis.theme(lengend_text_size = 16,lengend_title_size = 16,textsize = 18,leg = FALSE,pos_leg = "none")

pcol    <- cowplot::plot_grid(p1,p2,p3,p4,p5, align = "hv",ncol = 2, nrow = 3,labels = LETTERS[1:6])

pcol

```





\end{column}
\begin{column}{0.5\textwidth}

\scriptsize

```{r,echo=FALSE,warning=FALSE}

library(Metrics)
df_metrics <- data.frame(Smoother = c("Kernel", "Loess","Sp. Reg. Linear","Sp. Reg. Cúbico"),
                         `Parâm. Suavizador` = c(par.kernel,par.loess,as.integer(round(par.sp1,0)),as.integer(round(par.sp3,0))),
                         EQM      =  c( round(cv.min.kernel,4),round(cv.min.loess,4),round(cv.min.sp1,4),round(cv.min.sp3,4))
)

colnames(df_metrics) <- c("Suavizador","Parâm. Suavizador","EQM")
kable_data(data = df_metrics,cap = "\\label{tab:tab_eqm_aplicacao1} Erro Quadrático Médio ($EQM_{loocv}$) para os suavizadores Loess, Kernel e Splines de Regressão Linear e Cúbico (Aplicação 1). ",foot = NULL,c_names = c("Suavizador","Parâm. Suavizador","EQM"))


```

\scriptsize

```{r,echo=FALSE,warning=FALSE}

library(Metrics)
df_metrics <- data.frame(Smoother = c("Kernel", "Loess","Sp. Reg. Linear","Sp. Reg. Cúbico","Polinômio Cúbico"),
                         `Parâm. Suavizador` = c(par.kernel,par.loess,par.sp1,par.sp3,""),
                         EQM      =  c( round(rse(df$y,fitted.values(fit4)),6),
                                        round(rse(df$y,fitted.values(fit5)),6),
                                        round(rse(df$y,fitted.values(fit6)),6),
                                        round(rse(df$y,fitted.values(fit7)),6),
                                        round(rse(df$y,fitted.values(fit8)),6))
)

colnames(df_metrics) <- c("Suavizador","Parâm. Suavizador","EQM")
kable_data(data = df_metrics,cap = "\\label{tab:tab_eqm_aplicacao2_geral} Erro Quadrático Médio ($EQM_c$) para os suavizadores Loess, Kernel e Splines de Regressão Linear e Cúbico.",foot = NULL,c_names = c("Suavizador","Parâm. Suavizador","EQM"))

```




\end{column}
\end{columns}



