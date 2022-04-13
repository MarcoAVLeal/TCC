library(fitdistrplus)
library(ggplot2)
library(gridExtra)
library(knitr)
library(tidyverse)
library(kableExtra)
library(wesanderson)
library(RColorBrewer)



running.mean <- function(x = NULL, y = NULL, span = 5, align="full"){
  data <- data.frame(x = x,y = y)
  data <- data %>% arrange(x)
  estimates <- c()
  
  k=1
  j=span
  x_r <- c()
  for(i in 1:length(data[,1])){
    
    l <- length(data[,1])
    if(align=="left"){
      
      if((l-(i + (span))) > 0  ){
        estimates[i] <- mean(data[(i):(i+(span-1)),2])
        x_r[i] <- data$x[i]
        k = k + 1
      }
      
      
    }else if (align == "right"){
      if(i >= span){
        estimates[i] <- mean(data[(i):(i+(span-1)),2])
        x_r[i] <- data$x[i+(span-1)]
        k = k + 1
      }
    }else{
      if((i - span) < 0){
        
      }else if((l-(i + (span))) < 0  ){
        
      }else{
        
        estimates[k] <- mean(data[(i - (span - 1)):(i + (span - 1)),2])
        x_r[k] <- data$x[i]
        k = k + 1
      }
    }
  }
  d = data.frame(x=x_r,y=estimates)
  return(list(name = "Running Mean",fitted.values = estimates,span = span,data=d))
}


bin.smoother2  <- function(y,k){
  # y=dados$y
  # k=3
  library(reshape2)
  v   = split(y, ceiling(seq_along(y)/k))
  fit = sapply(v, function(item){rep(mean(item), length(item))})
  
  fit =  melt(data = fit,id.vars=1:k) %>% select(value)
  
  return(fit)
  
}

locv1 <- function(x1, y1, nd, span, ntrial)
{
  locvgcv <- function(sp, x1, y1)
  {
    nd <- length(x1)
    
    assign("data1", data.frame(xx1 = x1, yy1 = y1))
    fit.lo <- loess(yy1 ~ xx1, data = data1, span = sp, degree = 1, surface = "direct")
    res <- residuals(fit.lo)
    
    dhat2 <- function(x1, sp)
    {
      nd2 <- length(x1)
      diag1 <- diag(nd2)
      dhat <- rep(0, length = nd2)
      
      for(jj in 1:nd2){
        y2 <- diag1[, jj]
        assign("data1", data.frame(xx1 = x1, yy1 = y2))
        fit.lo <- loess(yy1 ~ xx1, data = data1, span = sp, family = "gaussian", degree = 2, surface = "direct")
        ey <- fitted.values(fit.lo)
        dhat[jj] <- ey[jj]
      }
      return(dhat)
    }
    
    dhat <- dhat2(x1, sp)
    trhat <- sum(dhat)
    sse <- sum(res^2)
    
    cv <- sum((res/(1 - dhat))^2)/nd
    gcv <- sse/(nd * (1 - (trhat/nd))^2)
    
    return(gcv)
  }
  
  gcv <- sapply(as.list(span1), locvgcv, x1 = x1, y1 = y1)
  #cvgcv <- unlist(cvgcv)
  #cv <- cvgcv[attr(cvgcv, "names") == "cv"]
  #gcv <- cvgcv[attr(cvgcv, "names") == "gcv"]
  
  return(gcv)
}



library(np)
cv_bws_npreg <- function(x,y,bandwidths=(1:50)/50,
                         num.folds=10) {
  require(np)
  n <- length(x)
  stopifnot(n> 1, length(y) == n)
  stopifnot(length(bandwidths) > 1)
  stopifnot(num.folds > 0, num.folds==trunc(num.folds))
  fold_MSEs <- matrix(0,nrow=num.folds,
                      ncol=length(bandwidths))
  colnames(fold_MSEs) = bandwidths
  case.folds <- rep(1:num.folds,length.out=n)
  case.folds <- sample(case.folds)
  for (fold in 1:num.folds) {
    train.rows = which(case.folds==fold)
    x.train = x[train.rows]
    y.train = y[train.rows]
    x.test = x[-train.rows]
    y.test = y[-train.rows]
    for (bw in bandwidths) {
      fit <- npreg(txdat=x.train,tydat=y.train,
                   exdat=x.test,eydat=y.test,bws=bw)
      fold_MSEs[fold,paste(bw)] <- fit$MSE
    }
  }
  CV_MSEs = colMeans(fold_MSEs)
  best.bw = bandwidths[which.min(CV_MSEs)]
  return(list(best.bw=best.bw,
              CV_MSEs=CV_MSEs,
              fold_MSEs=fold_MSEs))
}




textsize <<- 10
textsize2 <<- 18
point.size  <<- 3
point.alpha <<- .25
point.color <<- "black"
line.size  <<- 1.2
line.alpha <<- 1
line.color <<- "red"
cores               <<- c()
col.brew            <<- brewer.pal(n = 9, name = "Set1")
cores              <<- col.brew

mytheme <- theme(
                 axis.title = element_text(size = textsize),
                 axis.text = element_text(size = textsize),
                 legend.title = element_text(size = textsize),
                 legend.text = element_text(size = textsize))

axis_theme <- theme_bw()  +
  theme(
    axis.text.x = element_text(angle = 0,face = "bold",size = textsize),
    axis.text.y = element_text(angle = 0,face = "bold",size = textsize),
    legend.background = element_rect(fill = "transparent", colour = NA,size = 2),
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.title.x = element_text(colour = "black",size = textsize,face = "bold"),
    axis.title.y = element_text(colour = "black",size = textsize,face = "bold"),
    legend.title = element_text(colour = "black",size = 10),
    legend.position = "top",
    legend.text = element_text(colour = "black",size = 8,face = "bold"),
    panel.grid = element_line(linetype="dashed"),
    panel.grid.major = element_line(colour = "gray"),
    title =element_text(size=textsize, face='bold',hjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(color="#000000", face="bold", size=textsize,lineheight = 2))


axis.theme <- function(x.angle = 0,vjust=0,hjust=0.5,pos_leg="top",textsize = 10,lengend_title_size = 10,lengend_text_size = 8,title_size = 16,axis_x=F,axis_y=F){
  
  if(isTRUE(axis_x)){
    theme_bw()  +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(angle = 0,face = "bold",size = textsize),
        legend.background = element_rect(fill = "transparent", colour = NA,size = 2),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA),
        axis.title.x = element_text(colour = "black",size = textsize,face = "bold"),
        axis.title.y = element_text(colour = "black",size = textsize,face = "bold"),
        legend.title = element_text(colour = "black",size = lengend_title_size),
        legend.position = pos_leg,
        legend.text = element_text(colour = "black",size = lengend_text_size,face = "bold"),
        panel.grid = element_line(linetype="dashed"),
        panel.grid.major = element_line(colour = "gray"),
        title =element_text(size=title_size, face='bold',hjust = 0.5),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(color="#000000", face="bold", size=textsize,lineheight = 2))
    
    
  }else if(isTRUE(axis_y)){
    
    theme_bw()  +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle = x.angle,face = "bold",size = textsize,hjust=hjust, vjust=vjust),
        legend.background = element_rect(fill = "transparent", colour = NA,size = 2),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA),
        axis.title.x = element_text(colour = "black",size = textsize,face = "bold"),
        axis.title.y = element_text(colour = "black",size = textsize,face = "bold"),
        legend.title = element_text(colour = "black",size = lengend_title_size),
        legend.position = pos_leg,
        legend.text = element_text(colour = "black",size = lengend_text_size,face = "bold"),
        panel.grid = element_line(linetype="dashed"),
        panel.grid.major = element_line(colour = "gray"),
        title =element_text(size=title_size, face='bold',hjust = 0.5),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(color="#000000", face="bold", size=textsize,lineheight = 2))
    
    
  }else{
    theme_bw()  +
      theme(
        axis.text.x = element_text(angle = x.angle,face = "bold",size = textsize,hjust=hjust, vjust=vjust),
        axis.text.y = element_text(angle = 0,face = "bold",size = textsize),
        legend.background = element_rect(fill = "transparent", colour = NA,size = 2),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA),
        axis.title.x = element_text(colour = "black",size = textsize,face = "bold"),
        axis.title.y = element_text(colour = "black",size = textsize,face = "bold"),
        legend.title = element_text(colour = "black",size = lengend_title_size),
        legend.position = pos_leg,
        legend.text = element_text(colour = "black",size = lengend_text_size,face = "bold"),
        panel.grid = element_line(linetype="dashed"),
        panel.grid.major = element_line(colour = "gray"),
        title =element_text(size=title_size, face='bold',hjust = 0.5),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(color="#000000", face="bold", size=textsize,lineheight = 2))
    
    
    
  }
  
  
}


ci.compute<-function(means,mse,n,a,c.lev=95,txt,diff=FALSE){

statistic.val=qt(p=.5+c.lev/200, df = (n*a)-a )

if(isTRUE(diff)){
confidence.interval=round(means+(c(-1,1)*statistic.val*sqrt((2*mse)/n)),4)
}else{confidence.interval=round(means+(c(-1,1)*statistic.val*sqrt(mse/n)),4)}

return(confidence.interval)



}


plot_design <- function(data){
  
  
  library(ggplot2)
  ggplot(data = data, aes(x=data[,1], y=data[,2]))+
    geom_point(alpha=.5)+labs(x = "Tratamentos", y = "Resposta")+
    ggtitle("Gráfico contendo o comportamento do Experimento") +
    theme(plot.title = element_text(hjust = 0.5))+
    stat_summary(fun.y=mean, geom="line",aes(y = data[,2],group=1),
                 colour="red",lwd=0.8,alpha=.6)+
    stat_summary(fun.y=mean, geom="point",colour="blue",size=3,alpha=.8)
    
  
}


confidence.interval <- function(modelo){
  library(agricolae)
  library(ggplot2)
  aov      <- aov(modelo)
  anova    <- anova(modelo)
  hsd.test <- HSD.test(y = aov, trt = "V1",console =  F,group = T)
  mse      <- anova$`Mean Sq`[2]
  means    <- sort(hsd.test$means[,1],decreasing = T)
  a        <- length(unique(modelo$model[,(attr(modelo$terms,which = "dataClasses")=="factor")]))
  n        <- length(modelo$model[,(attr(modelo$terms,which = "dataClasses")=="factor")])/a
   # unique(modelo$model[,(attr(modelo$terms,which = "dataClasses")=="factor")])
  
  intervalos <- sapply(X = means,ci.compute,mse,n,a)
  intervalos <- as.data.frame(t(rbind(means,intervalos)),rownames=T)
  label   <- rownames(HSD.test(y = aov, trt = "V1",console =  F,group = T)$groups)
  colnames(intervalos) <- c("means","LI","LS")
  row.names(intervalos) <- label
  require(ggplot2)

  x_axis  <- sort(as.numeric(label))
  plot <- ggplot() +
    geom_point(size = 1.8,aes(x = x_axis,y = means )) +
    geom_errorbar(aes(x = x_axis,
                      ymax = intervalos$LS,
                      ymin = intervalos$LI))+
    labs(x = "Tratamentos", y = "Resposta")+
    ggtitle("Gráfico Comparativo dos Intervalos de Confianças dos Tratamentos") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(face="bold", 
                                     size=10, angle=45),
          axis.text.y = element_text(face="bold", 
                                     size=10, angle=45))+ 
          scale_x_discrete(limits=x_axis, labels=label)
  #plot
  
  
  return(list(intervalos,plot))
  
}

qq_gplot <- function(data){

fit    <- fitdist(data[["values"]], distr = "norm")
empir  <- ecdf(data[["values"]])
Fempir <- empir(data[["values"]])
Fteori <- pnorm(data[["values"]], fit$estimate[1],fit$estimate[2])

Qteori <-sapply(Fteori, function(x) qnorm(x, fit$estimate[1],fit$estimate[2]))
Qempir <-sapply(Fempir, function(x) qnorm(x, fit$estimate[1],fit$estimate[2]))

ggplot()+
geom_point(aes(x=Qteori, y=Qempir),size=2)+labs(x = "Distribuição Teórica", y = "Distribuição Empírica")+
xlim(range(data[["values"]])+c(-0,+0))+
ylim(range(data[["values"]])+c(-0,+0))+
ggtitle(paste0("Gráfico de Dispersão QQplot ",data[["label"]]))+
geom_abline(slope = 1,intercept = 0,col="red",size=1.20)+
theme(plot.title = element_text(hjust = 0.5))
}

fit_gplot <- function(data,axis=c(1,2),grau){
  ggplot(data = data, aes(x = as.numeric(data[,axis[1]]),y = data[,axis[2]]))+
    geom_point(alpha=.5)+
    labs(x = toupper(colnames(data)[1]), y = toupper(colnames(data)[2]))+
    ggtitle(paste0("Regressão Polinomial de ordem ", grau))+
    stat_smooth(method = "lm", formula = y ~ poly(x,grau),col="red",se = F)+
    theme(plot.title = element_text(hjust = 0.5))
}



summary_descritive <-function(dados){
  library(tidyverse)
  library(knitr)
  library(kableExtra)
  
  min    <- min(dados)
  q1     <- quantile(dados)[2]
  median <- quantile(dados)[3]
  q3     <- quantile(dados)[4]
  mean   <- mean(dados)
  max    <- max(dados)
  
  data <- data.frame(Min.=min,"1st Qu."=q1,Median=median,Mean=mean,"3rd Qu."=q3, Max=max,row.names = "Value")
  
  data  %>%
    kable(booktabs=T,caption = "Tabela contendo medidas descritivas",digits = 4,col.names = c("Min.","1st Qu.","Median","Mean","3rd Qu","Max."),align = "c")%>%
    kable_styling(full_width = F, latex_options = "hold_position") %>%
    row_spec(0, align = "c",bold=T ) %>%
    column_spec(1, bold = T)
}


boxggplot <- function(data,title=" dos Resíduos",xtitle="Resíduos",ytitle="Valores"){
  library(ggplot2)
  outliers<-0L
  #outliers <- as.numeric(boxplot.stats(data)$out)
  outliers <- data[IsOutlier(data = data)]
  if(length(outliers)>0){
    pos <- which(data %in% outliers)
    out <- data.frame(Pos=pos,outliers=outliers)
  }else{
    out=data
    pos <- NA
    outliers <-NA
    label <- NA
  }
  
 plot <- ggplot(data = as.data.frame(data),aes(x=xtitle, y=data))+
    geom_boxplot(outlier.colour="black",alpha=0.8, outlier.shape=16,outlier.size=2, notch=FALSE)+
   
    labs(x = NULL,y = ytitle) +
    ggtitle(paste0(title))+
    axis_theme
   return(list(plot=plot,outliers=out)) 
  
}

r_standard <- function(modelo){
  
 resid(modelo)/summary(modelo)$sigma
}


r_student <- function(modelo){
  X=model.matrix(modelo)  
  h=hat(x = X,intercept = T)
 resid(modelo)/(summary(modelo)$sigma*(1-h)^.5)
}

r_student_external <- function(modelo){
X=model.matrix(modelo)  
p=length(coef(modelo))
n=length(modelo$residuals)
rstudent(modelo)*((n-p-1)/(n-p-rstudent(modelo)))^.5

}



envelope=function(modelo,title=NULL){
  dados=na.omit(modelo$data)
  nsim=100
  n=modelo$df.null+1
  r1=sort(rstandard(modelo,type='deviance'))
  m1=matrix(0,nrow=n,ncol=nsim)
  a2=simulate(modelo,nsim=nsim)
  
  for (i in 1:nsim){
    dados$y=a2[,i]
    aj=update(modelo,y~.,data=dados)
    m1[,i]=sort(rstandard(aj,type='deviance'))}
  
  li=apply(m1,1,quantile,0.025)
  m=apply(m1,1,quantile,0.5)
  ls=apply(m1,1,quantile,0.975)
  
  quantis=qnorm((1:n-0.5)/n)
  
  
  plot <- ggplot()+
    geom_point(aes(x=quantis, y=r1), alpha=.5)+
    labs(x = 'Percentil da N(0,1)', y = 'Resíduos')+
    ggtitle(title)+
    theme(plot.title = element_text(hjust = 0.5))+
    geom_line(aes(x = quantis,y=li,type="l"),size=0.5)+
    geom_line(aes(x = quantis,y=m,type="l"),size=.5,colour="red")+
    geom_line(aes(x = quantis,y=ls,type="l"),size=.5)
  
  
  plot
}


residual_qqnorm <- function(data,title= title){
  
  ggplot(size=.75)+
    stat_qq(aes(sample=data), alpha=.5)+
    stat_qq_line(aes(sample=data), size=.75,colour="red")+
    labs(x = "Quantís Teóricos", y = "Quantís Empíricos")+
    ggtitle(paste0(title))+
    axis_theme
}



plot_resvfit <- function(x,y,labelx="Eixo X",labely="Eixo Y", title="Gráfico de Dispersão"){
 
  ggplot()+
    geom_point(aes(x=x, y=y),alpha=.5)+labs(x = labelx, y = labely)+
    xlim(range(x)+c(-0,+0))+
    ylim(range(y)+c(-0,+0))+
    geom_abline(slope = 0,intercept = 0,col="red",size=.75, alpha=.5)+
    ggtitle(paste0(title))+
    axis_theme
  
  
}





# kable_data <- function(data,cap,foot=" ",align="c"){
#   library(kableExtra)
#   
#   data %>%
#     kable(booktabs=T,caption = cap,align = align) %>%
#     add_footnote(foot) %>%
#     kable_styling(full_width = F, latex_options = "hold_position") %>%
#     row_spec(0, align = align,bold=T ) %>%
#     column_spec(1, bold = T)
#   
# }

kable_data <- function(data,cap,foot=" ",align="c",label=NULL){
  library(kableExtra)
  
  data %>%
    kable(booktabs=T,caption = cap,align = align,label = label) %>%
    add_footnote(foot) %>%
    kable_styling(full_width = F, latex_options = "hold_position") %>%
    row_spec(0, align = align ) %>%
    column_spec(1)
  
}



kable_cov <- function(data,cap){
  library(kableExtra)
  
  
  cov(data) %>%
    kable(booktabs=T,caption = cap,digits = 4)%>%
    kable_styling(full_width = F, latex_options = "hold_position") %>%
    row_spec(0, align = "c",bold=T ) %>%
    column_spec(1, bold = T)
  
  
}


kable_cor <- function(data,cap){
  
  
  cor(data) %>%
    kable(booktabs=T,caption = cap,digits = 4)%>%
    kable_styling(full_width = F, latex_options = "hold_position") %>%
    row_spec(0, align = "c",bold=T ) %>%
    column_spec(1, bold = T)
  
  
}


kable_estimates <- function(model,cap,label= NULL){
  library(kableExtra)
  if(any(apply(modelo$model,2,typeof)=="character")){

    vif=c(" ")



  } else if (length(modelo$coefficients) > 2){
    vif=round(car::vif(modelo),4)

  }
  else{vif=c(" ")}


  names <- names(model$coefficients)
  i=1
  col_label<-c(length(names))
  for(item in names){

    paste = paste0("$\\mathbf{",item,"}$")
    col_label[i]=paste
    i=i+1
  }
  #col_labelf


  summary <- summary(model)
  anova   <- anova(model)
  indices <- c(summary$r.squared,summary$adj.r.squared,summary$fstatistic[1],
                        pf(q = summary$fstatistic[1],df1 = summary$fstatistic[2],summary$fstatistic[3],lower.tail = F ))
  res <- cbind(round(summary$coefficients,4),VIF=c("NA",vif))
  row.names(res) <- col_label
  res <- rbind(res,c("$\\mathbf{R^2}$","$\\mathbf{R^2_{Adj}}$","$\\mathbf{F_0}$","$\\mathbf{Pr(>|F|)}$","NA"),c(round(indices,4), "NA"))
  n <- length(model$coefficients)



     res %>%
    kable(booktabs=T,caption = cap,digits = 4,escape = F,label=label)%>%
    kable_styling(full_width = F, latex_options = "hold_position") %>%
    row_spec(0, align = "l" ) %>%
       row_spec(n+1, align = "c" ) %>%
       row_spec(n+2, align = "c" ) %>%
    column_spec(1)



}

# kable_estimates <- function(model,cap){
#   library(kableExtra)
#   if(any(apply(modelo$model,2,typeof)=="character")){
#     
#     vif=c(" ")
#     
#     
#     
#   } else if (length(modelo$coefficients) > 2){
#     vif=round(car::vif(modelo),4)
#     
#   }
#   else{vif=c(" ")}
#   
#   
#   names <- names(model$coefficients) 
#   i=1
#   col_label<-c(length(names))
#   for(item in names){
#     
#     paste = paste0("$\\mathbf{",item,"}$")
#     col_label[i]=paste
#     i=i+1
#   }
#   #col_label
#   
#   
#   summary <- summary(model)
#   anova   <- anova(model)
#   indices <- c(summary$r.squared,summary$adj.r.squared,summary$fstatistic[1],
#                pf(q = summary$fstatistic[1],df1 = summary$fstatistic[2],summary$fstatistic[3],lower.tail = F ))
#   res <- cbind(round(summary$coefficients,4),VIF=c("NA",vif))
#   row.names(res) <- col_label
#   res <- rbind(res,c("$R^2$","$R^2_{Adj}$","$F_0$","$Pr(>|F|)$", "NA"),c(round(indices,4), "NA"))
#   n <- length(model$coefficients)
#   
#   
#   
#   res %>%
#     kable(caption = cap,digits = 4,escape = F)%>%
#     kable_styling(full_width = F, latex_options = "hold_position") %>%
#     row_spec(0, align = "l" ) %>%
#     row_spec(n+1, align = "c" ) %>%
#     row_spec(n+2, align = "c" ) %>%
#     column_spec(1)
#   
#   
#   
# }


kable_fit <- function(model){
  summary <- summary(model)
  indices <- data.frame(R2=summary$r.squared,R2Adj=summary$adj.r.squared,Fstatistic=summary$fstatistic[1],pvalor=pf(q = summary$fstatistic[1],df1 = summary$fstatistic[2],summary$fstatistic[3],lower.tail = F ))
  kable(indices,booktabs=T,caption = "Coeficiente de Determinação e Estatística para ANOVA",digits = 4)%>%
    kable_styling(full_width = F, latex_options = "hold_position") %>%
    row_spec(0, align = "c",bold=T ) %>%
    column_spec(1, bold = T)
  
}
################################################## ########################################
#################################### FUNÇÕES LATEX ########################################
###########################################################################################

print_ajuste <- function(coeficientes){
  coeficientes <- round(x = coeficientes,digits = 4)
  #coeficientes <- as.character(coeficientes)
  expressao <- c()
  i=1
  sinal <- NULL
  for(item in coeficientes){
    
    if(item > 0 && i==1){
      
      expressao[i] <- paste0(" + ",item)
      
    }else if (item < 0 && i==1){
      
      expressao[i] <- paste0(" - ",abs(item))
      
    }else if(item > 0 && i>1){
      
      expressao[i] <- paste0(" + ",item,"X_",i-1) 
      
    }else {
      
      expressao[i] <- paste0(" - ", abs(item),"X_",i-1)
      
    }
    i=i+1
  }
  exp <- "$$ \\hat{Y}="
  for ( item in expressao){
    
    exp <- str_c(exp,item)
    
  }
  exp <- str_c(exp," $$")
  exp %>% knitr::raw_latex()
  
  
  
  
}


write_matex <- function(x) {
  begin <- "$$\\begin{bmatrix}"
  end <- "\\end{bmatrix}$$"
  X <-
    apply(x, 1, function(x) {
      paste(
        paste(x, collapse = "&"),
        "\\\\"
      )
    })
  writeLines(c(begin, X, end))
}


write_matex2 <- function(x) {
  begin <- "\\begin{bmatrix}"
  end <- "\\end{bmatrix}"
  X <-
    apply(x, 1, function(x) {
      paste(
        paste(x, collapse = "&"),
        "\\\\"
      )
    })
  paste(c(begin, X, end), collapse = "")
}

################################################## ########################################
################################ FIM FUNÇÕES LATEX ########################################
###########################################################################################
IsOutlier <- function(data) {
  lowerq = quantile(data, na.rm = TRUE)[2]
  upperq = quantile(data, na.rm = TRUE)[4]
  iqr = upperq - lowerq 
  threshold_upper = (iqr * 1.5) + upperq
  threshold_lower = lowerq - (iqr * 1.5)
  data > threshold_upper | data <  threshold_lower 
}





################################ MODELOS ADITIVOS ##########################################


plot.curves <- function(data  = NULL, x,y,labelx="Eixo X",labely="Eixo Y",
                        title = "Gráfico de Dispersão",
                        type  = NULL,
                        se    = FALSE,
                        span  = 3/4,
                        mult  = NULL,
                        fit   = NULL){
  
  data <- data.frame(x=x,y=y)
  data <- data %>% arrange(x)
  plot <- ggplot(data)+
    geom_point(aes(x=x, y=y),alpha=point.alpha,size=point.size, color = "black")+
    labs(x = labelx, y = labely)+
    xlim(range(x)+c(-0,+0))+
    ylim(range(y)+c(-0,+0))+
    ggtitle(paste0(title))+
    axis.theme()
    
    
    if(is.null(type)){
      
      plot
      
      }else if(type == "res"){
        plot +
          geom_abline(slope = 0,intercept = 0,col=line.color,size=line.size,alpha=line.alpha)+
          ggtitle(paste0(title," \nResíduos vs Valores Ajustados"))
        
        }else if(type == "lm"){
          
          plot +
            geom_smooth(method = lm,se = se, col=line.color,size=line.size, alpha=line.alpha, aes(x = x,y = y))+
            ggtitle(paste0(title," \nCom ajuste Linear"))
          
          }else if(type == "loess"){
            
            plot + 
              stat_smooth(method = "loess",se = se, col=line.color,size=line.size, alpha=line.alpha,span = span, aes(x = x,y = y))+
              ggtitle(paste0(title," \nCom ajuste Loess","\n span = ",round(span,2)))
            
            }else if(type == "mult"){
              library(wesanderson)
              library(RColorBrewer)
              
              colors   <- as.list(brewer.pal(7,"Dark2"))
              n_colors <- c("C1","C2","C3","C4","C5","C6","C7")
              names(colors) <- n_colors
              l <- length(mult)
              n_colors1 <-c()
              for(item in mult){
                n_colors1[i] <- paste(n_colors[i], "\nSpan : ", item$span,"\nType : ",item$type)
                names(colors) <- n_colors1
                plot <- ggplot()+
                  geom_point(aes(x=x, y=y),alpha=point.alpha,size=point.size,color=point.color)+
                  labs(x = labelx, 
                       y = labely,
                       color = "Legend")+
                  xlim(range(x)+c(-0,+0))+
                  ylim(range(y)+c(-0,+0))+
                  ggtitle(paste0(title," \nCom ajuste Loess"))+
                  
                  axis.theme()
                size = 1
                alpha = 1
                if(l <= 1){
                  
                  plot <- plot + 
                    stat_smooth(method = mult[[1]]$type,se = se,size=size, alpha=alpha,span = mult[[1]]$span, aes(x = x,y = y,color = n_colors1[1]))
                  }else if(l <= 2){
                    
                    plot <- plot + 
                      stat_smooth(method = mult[[1]]$type,se = se,size=size, alpha=alpha,span = mult[[1]]$span, aes(x = x,y = y,color = n_colors1[1]))+ 
                      stat_smooth(method = mult[[2]]$type,se = se,size=size, alpha=alpha,span = mult[[2]]$span, aes(x = x,y = y,color = n_colors1[2]))
      
                  }else if(l <= 3){
                    plot <- plot + 
                      stat_smooth(method = mult[[1]]$type,se = se,size=size, alpha=alpha,span = mult[[1]]$span, aes(x = x,y = y,color = n_colors1[1]))+ 
                      stat_smooth(method = mult[[2]]$type,se = se,size=size, alpha=alpha,span = mult[[2]]$span, aes(x = x,y = y,color = n_colors1[2]))+ 
                      stat_smooth(method = mult[[3]]$type,se = se,size=size, alpha=alpha,span = mult[[3]]$span, aes(x = x,y = y,color = n_colors1[3]))
                    
                  }else if(l <= 4){
                    plot <- plot + 
                      stat_smooth(method = mult[[1]]$type,se = se,size=size, alpha=alpha,span = mult[[1]]$span, aes(x = x,y = y,color = n_colors1[1]))+ 
                      stat_smooth(method = mult[[2]]$type,se = se,size=size, alpha=alpha,span = mult[[2]]$span, aes(x = x,y = y,color = n_colors1[2]))+ 
                      stat_smooth(method = mult[[3]]$type,se = se,size=size, alpha=alpha,span = mult[[3]]$span, aes(x = x,y = y,color = n_colors1[3]))+ 
                      stat_smooth(method = mult[[4]]$type,se = se,size=size, alpha=alpha,span = mult[[4]]$span, aes(x = x,y = y,color = n_colors1[4]))
                    
                  }else if(l <= 5){
                    plot <- plot + 
                      stat_smooth(method = mult[[1]]$type,se = se,size=size, alpha=alpha,span = mult[[1]]$span, aes(x = x,y = y,color = n_colors1[1]))+ 
                      stat_smooth(method = mult[[2]]$type,se = se,size=size, alpha=alpha,span = mult[[2]]$span, aes(x = x,y = y,color = n_colors1[2]))+ 
                      stat_smooth(method = mult[[3]]$type,se = se,size=size, alpha=alpha,span = mult[[3]]$span, aes(x = x,y = y,color = n_colors1[3]))+ 
                      stat_smooth(method = mult[[4]]$type,se = se,size=size, alpha=alpha,span = mult[[4]]$span, aes(x = x,y = y,color = n_colors1[4]))+ 
                      stat_smooth(method = mult[[5]]$type,se = se,size=size, alpha=alpha,span = mult[[5]]$span, aes(x = x,y = y,color = n_colors1[5]))
                    
                  }else if(l <= 6){
                    plot <- plot + 
                      stat_smooth(method = mult[[1]]$type,se = se,size=size, alpha=alpha,span = mult[[1]]$span, aes(x = x,y = y,color = n_colors1[1]))+ 
                      stat_smooth(method = mult[[2]]$type,se = se,size=size, alpha=alpha,span = mult[[2]]$span, aes(x = x,y = y,color = n_colors1[2]))+ 
                      stat_smooth(method = mult[[3]]$type,se = se,size=size, alpha=alpha,span = mult[[3]]$span, aes(x = x,y = y,color = n_colors1[3]))+ 
                      stat_smooth(method = mult[[4]]$type,se = se,size=size, alpha=alpha,span = mult[[4]]$span, aes(x = x,y = y,color = n_colors1[4]))+ 
                      stat_smooth(method = mult[[5]]$type,se = se,size=size, alpha=alpha,span = mult[[5]]$span, aes(x = x,y = y,color = n_colors1[5]))+ 
                      stat_smooth(method = mult[[6]]$type,se = se,size=size, alpha=alpha,span = mult[[6]]$span, aes(x = x,y = y,color = n_colors1[6]))
                    
                  }else if(l <= 7){
                    plot <- plot + 
                      stat_smooth(method = mult[[1]]$type,se = se,size=size, alpha=alpha,span = mult[[1]]$span, aes(x = x,y = y,color = n_colors1[1]))+ 
                      stat_smooth(method = mult[[2]]$type,se = se,size=size, alpha=alpha,span = mult[[2]]$span, aes(x = x,y = y,color = n_colors1[2]))+ 
                      stat_smooth(method = mult[[3]]$type,se = se,size=size, alpha=alpha,span = mult[[3]]$span, aes(x = x,y = y,color = n_colors1[3]))+ 
                      stat_smooth(method = mult[[4]]$type,se = se,size=size, alpha=alpha,span = mult[[4]]$span, aes(x = x,y = y,color = n_colors1[4]))+ 
                      stat_smooth(method = mult[[5]]$type,se = se,size=size, alpha=alpha,span = mult[[5]]$span, aes(x = x,y = y,color = n_colors1[5]))+ 
                      stat_smooth(method = mult[[6]]$type,se = se,size=size, alpha=alpha,span = mult[[6]]$span, aes(x = x,y = y,color = n_colors1[6]))+ 
                      stat_smooth(method = mult[[7]]$type,se = se,size=size, alpha=alpha,span = mult[[7]]$span, aes(x = x,y = y,color = n_colors1[7]))
                    
                  }
              }
                plot <- plot + scale_color_manual(values = colors)
                
                }else if (type == 'rm'){
                  
                  plot + 
                    geom_line(data = fit$data,aes(x = fit$data$x ,y = fit$data$y), col=line.color,size=line.size,alpha=line.alpha) +
                    ggtitle(paste0(title," \n",fit$name,"\n span = ",round(fit$span,2)))
                }else if (type == 'rl'){
                  
                  plot + 
                    geom_line(aes(x, fit$fitted.values), col=line.color,size=line.size,alpha=line.alpha) +
                    ggtitle(paste0(title," \n",fit$name,"\n span = ",round(fit$span,2)))
                }else if (type == 'g'){
                  
                  plot + 
                    geom_line(aes(x, fit), col=line.color,size=line.size,alpha=line.alpha)
                  # +
                  #   ggtitle(paste0(title," \n",fit$name,"\n span = ",round(fit$span,2)))
                }
  
  
  
}




plot.mult.curves <- function(df,df_fit = NULL,labelx="Eixo X",labely="Eixo Y",
                        title = "Gráfico de Dispersão",legend.pos = "right",line.s = 1.2,alpha.o = 1,point.alpha = .25,
                        point.color = "black"
                        ){
  
 
  plot <- ggplot(data = df,aes(x=x,y=y))+
    geom_point(alpha=point.alpha,size=point.size, color = point.color)+
    labs(x = labelx, 
         y = labely,
         color = "Legenda")+
    xlim(range(x)+c(-0,+0))+
    ylim(range(y)+c(-0,+0))+
    ggtitle(paste0(title))+
    geom_line(data = df_fit,aes(x=x,y = value,color=variable),size=line.s,alpha=alpha.o) +
    scale_color_manual(values = cores) + axis.theme(pos_leg =legend.pos)
  
  
  plot
  
  
  
}




backfitting <- function(x,y,tol   = 0.001*sd(y)/sqrt(length(y)),span  = c(3/4,3/4),bdw  = c(10,10), k = c(5,5), knots = NULL,method = "loess"){
  n     <- length(y)
  fx_11 <- numeric(n)
  fx_21 <- numeric(n)
  
  
  alpha     <- mean(y)
  converged <- FALSE
  
  while( !converged ){
    
    fx_10           <- fx_11
    fx_20           <- fx_21
    
    if(method=="loess"){
      y_x2            <- y - alpha - fx_11
      fx_21[x2_ord]   <- lowess(x = x[,2],y = y_x2,f = span[2])$y
      y_x1            <-  y - alpha - fx_21
      fx_11[x1_ord]   <- lowess(x = x[,1],y_x1,f = span[1])$y
      
    }else if( method == "kernel"){
      y_x2            <- y - alpha - fx_11
      fx_21[x2_ord]   <- ksmooth(x = x[,2],y = y_x2,kernel = "normal",bandwidth =  bdw[2])$y
      y_x1            <-  y - alpha - fx_21
      fx_11[x1_ord]   <- ksmooth(x = x[,1],y = y_x1,kernel = "normal",bandwidth = bdw[1])$y
      
    }else if( method == "sr.linear"){
      y_x2            <- y - alpha - fx_11
      fx_21[x2_ord]   <- additive.spline.linear(x = x[,2],y = y_x2,k = k[2])$fitted.values
      y_x1            <-  y - alpha - fx_21
      fx_11[x1_ord]   <- additive.spline.linear(x = x[,1],y = y_x1,k = k[1])$fitted.values
      
    }else if( method == "sr.cubic"){
      y_x2            <- y - alpha - fx_11
      fx_21[x2_ord]   <- additive.spline.cubic(x = x[,2],y = y_x2,k = k[2],knots = knots$v2)$fitted.values
      y_x1            <-  y - alpha - fx_21
      fx_11[x1_ord]   <- additive.spline.cubic(x = x[,1],y = y_x1,k = k[1],knots = knots$v1)$fitted.values
      
    }
    
    Norma_f1_0 <- sqrt(t(fx_10)%*%(fx_10))
    Norma_f2_0 <- sqrt(t(fx_20)%*%(fx_20))
    
    Norma_dif_f1 <- sqrt(t(fx_11 - fx_10)%*%(fx_11 - fx_10))
    Norma_dif_f2 <- sqrt(t(fx_21 - fx_20)%*%(fx_21 - fx_20))
    
    Dif_rel <- (Norma_dif_f1 + Norma_dif_f2)/(Norma_f1_0 + Norma_f2_0) # Diferen?a Relativa
    if(Dif_rel > tol){
      
      converged = TRUE
      
    }
    
    return(list(v1=fx_11,v2=fx_21,v3=y_x1,v4=y_x2))
  }
}



backfitting2 <- function(x,y,tol   = 0.001*sd(y)/sqrt(length(y)),span  = c(3/4,3/4),bdw  = c(10,10), k = c(5,5), knots = list(v1=NULL,v2=NULL),method = "loess"){
  n     <- length(y)
  fx_11 <- numeric(n)
  fx_21 <- numeric(n)
  
  
  alpha     <- mean(y)
  converged <- FALSE
  
  while( !converged ){
    
    fx_10           <- fx_11
    fx_20           <- fx_21
    
    if(method=="loess"){
      y_x2            <- y - alpha - fx_11
      fx_21[x2_ord]   <- lowess(x = x[,2],y = y_x2,f = span[2])$y
      y_x1            <-  y - alpha - fx_21
      fx_11[x1_ord]   <- lowess(x = x[,1],y_x1,f = span[1])$y
      
    }else if( method == "kernel"){
      y_x2            <- y - alpha - fx_11
      fx_21[x2_ord]   <- ksmooth(x = x[,2],y = y_x2,kernel = "normal",bandwidth =  bdw[2])$y
      y_x1            <-  y - alpha - fx_21
      fx_11[x1_ord]   <- ksmooth(x = x[,1],y = y_x1,kernel = "normal",bandwidth = bdw[1])$y
      
    }else if( method == "sr.linear"){
      y_x2            <- y - alpha - fx_11
      fx_21[x2_ord]   <- additive.spline.linear(x = x[,2],y = y_x2,k = k[2])$fitted.values
      y_x1            <-  y - alpha - fx_21
      fx_11[x1_ord]   <- additive.spline.linear(x = x[,1],y = y_x1,k = k[1])$fitted.values
      
    }else if( method == "sr.cubic"){
      y_x2            <- y - alpha - fx_11
      fx_21[x2_ord]   <- additive.spline.cubic(x = x[,2],y = y_x2,k = k[2],knots = knots$v2)$fitted.values
      y_x1            <-  y - alpha - fx_21
      fx_11[x1_ord]   <- additive.spline.cubic(x = x[,1],y = y_x1,k = k[1],knots = knots$v1)$fitted.values
      
    }
    
    Norma_f1_0 <- sqrt(t(fx_10)%*%(fx_10))
    Norma_f2_0 <- sqrt(t(fx_20)%*%(fx_20))
    
    Norma_dif_f1 <- sqrt(t(fx_11 - fx_10)%*%(fx_11 - fx_10))
    Norma_dif_f2 <- sqrt(t(fx_21 - fx_20)%*%(fx_21 - fx_20))
    
    Dif_rel <- (Norma_dif_f1 + Norma_dif_f2)/(Norma_f1_0 + Norma_f2_0) # Diferen?a Relativa
    if(Dif_rel > tol){
      
      converged = TRUE
      
    }
    
    return(list(y.x1=(y - alpha - (fx_21 - mean(fx_21))) - mean((y - alpha - (fx_21 - mean(fx_21)))),
                y.x2=(y - alpha - (fx_11 - mean(fx_11))) - mean((y - alpha - (fx_11 - mean(fx_11)))),
                fx_11 = fx_11 - mean(fx_11),
                fx_21 = fx_21 - mean(fx_21)))
  }
}
