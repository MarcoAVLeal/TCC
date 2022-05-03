###  LOOCV
df        <- cbind(x= dados$x,y = dados$y) %>% as.data.frame
tr= 1:nrow(df)
dftr= df[tr,]
number_of_bins   = seq(0.05,0.95,0.01)            
number_of_bins_sp= 1:length(number_of_bins) 
cv.errors.loess  = matrix( NA,nrow =  nrow(df), ncol = length(number_of_bins) )
cv.errors.kernel = matrix( NA,nrow =  nrow(df), ncol = length(number_of_bins) )
cv.errors.sp1 = matrix( NA,nrow =  nrow(df), ncol = length(number_of_bins_sp) )
cv.errors.sp3 = matrix( NA,nrow =  nrow(df), ncol = length(number_of_bins_sp) ) 


for( i in 1:length(number_of_bins)){ # Laço de repetição que irá percorrer todos os parâmetros possíveis.
  
  for( j in tr ){ # Laço de repetição que irá percorrer todas as observações.
    
    
    fit.loess              <- loess(y ~ x, span = number_of_bins[i], data=df[tr!=j,]) # Ajuste para i-ésimo parâmetro, retirando  a j-ésima observação
    loess.pred             = predict( fit.loess, newdata=data.frame(x = df$x[tr==j]) ) # PRedicação para j-ésima observação
    # Percorre todas as observações armazenando os EQM's na matrix abaixo, para o i-ésimo parâmetro (coluna) na j-ésima observação(linha).
    # Ao final do laço de repetição cada coluna (parâmetro) será um vetor contendo "j" valores para os EQM's.
    cv.errors.loess[j,i]   = mean( ( df$y[tr==j] - loess.pred )^2,na.rm = TRUE ) # Marix que contém todos os EQM's possíveis (por observação)
    
    
  }}

cv.errors.mean.sp1   = apply(cv.errors.sp1,2,mean,na.rm = TRUE) # Aplica a média por coluna para obter o ERRO qudrático médio final,ter-se-a um EQM para cada parâmetro de suavização.
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