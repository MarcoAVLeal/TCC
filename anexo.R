library(igraph)
library(zoo)
widths    <- c(0.05,0.15,0.25,0.35,.45)   ## Diferentes tamanha de span em um vetor
df        <- cbind(dados$x,dados$y)  # cria um data frame contendo uma coluna com os valores de x e outra para y
for(size in widths){        ## loop for, para cada span, faz um ajuste e concantena no dataframe atraves do cbind
  
  fit     <-  lowess(x = dados$x,y = dados$y,f = size)$y  ## valores ajustados para a técnica lowess, e paga somente os valores y
  df      <-  cbind(df,fit)  ## concatena por coluna dentro do data frame para cada ajuste com span diferentes
}


### Renomeias as colunas do dataframe df, as duas primeira sao os valores originais dos dados simulados. para cada span temos uma coluna diferente para valores ajustados no dataframe df, logo temos atribuimos um nome especifico para cada coluna.
colnames(df) <- c("x","y",
                  paste("Ajuste1 - ", "Largura : ",widths[1]),
                  paste("Ajuste2 - ", "Largura : ",widths[2]),
                  paste("Ajuste3 - ", "Largura : ",widths[3]),
                  paste("Ajuste4 - ", "Largura : ",widths[4]),
                  paste("Ajuste5 - ", "Largura : ",widths[5]))

## procedimento para colocar os valores ajustados empilhados em um unica coluna e criasse uma coluna com os rotulos para os valores ajustados
df <- as.tibble(df) %>%
  gather(key = "variable", value = "value",-x,-y)


## funcao do arquivo funcoes.R que na pasa Relatorio final
plot.mult.curves(df = df,title = "Gráfico de dispersão com ajustes Lowess")


