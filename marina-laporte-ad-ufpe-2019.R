## Script do trabalho final - Análise de Dados 
## Marina Laporte 

## Definir diretório 

setwd('C:/Users/marin/OneDrive/Área de Trabalho/Trabalho Marina')

## Ler arquivo em .csv 

library(readxl)

bd <- read_excel("dados_indice_2.xlsx")

## Tranformar variáveis em Z score (Como recomendado por Hair)

bd$delegacia_padron <- scale(bd$delegacia, center = TRUE, scale = TRUE)

bd$centros_especializados_padron <- scale(bd$centros_especializados, center = TRUE, scale = TRUE)

bd$Abrigamento_padron <- scale(bd$Abrigamento, center = TRUE, scale = TRUE) 

bd$varas_juizados_padron <- scale(bd$varas_juizados, center = TRUE, scale = TRUE)

bd$promotoria_padron <- scale(bd$Promotoria, center = TRUE, scale = TRUE)

bd$Defensoria_padron <- scale(bd$Defensoria, center = TRUE, scale = TRUE) 

## Criar um data frame com as variáveis padronizadas 

indice <- data.frame(bd$delegacia_padron,bd$centros_especializados_padron, bd$Abrigamento_padron,bd$varas_juizados_padron,
                     bd$promotoria_padron, bd$Defensoria_padron)

## Adequabilidade dos dados para à ACP

### Teste de Kaiser-Meyer-Olkin

KMO(fit) ###  0.71

### Teste de Alfa de Cronbach 

alpha(indice)  ###  0.848

### Teste de esferecidade 

bart<-function(dat){ #dat is your raw data
  R<-cor(dat)
  p<-ncol(dat)
  n<-nrow(dat)
  chi2<- -((n-1)-((2*p)+5)/6 ) * log(det(R)) #this is the formula
  df<-(p*(p-1)/2)
  crit<-qchisq(.95,df) #critical value
  p<-pchisq(chi2,df,lower.tail=F) #pvalue
  cat("Bartlett's test of sphericity: X2(",
      df,")=",chi2,", p=",
      round(p,3),sep="" )   
}

bart(indice) ###  p-value < 0.005

## Criar o índice 

## Análise de Compotentes Principais 

fit <- princomp(indice, cor=TRUE) ### inserindo dados brutos e extraindo CPs da matriz de correlação

summary(fit) ### Resumo 

loadings(fit) ### Loadings dos CPs

plot(fit,type="lines") ### scree plot

fit$scores ## Componentes principais 

## Extrair e rotacionar fatores 

## Varimax 

install.packages('psych')
library(psych)

fit <- principal(indice, nfactors=6, rotate="varimax") ### Gerando fatores 

## Extrair 1 fator --> analisando o scree plot 

fit <- principal(indice, nfactors=1, rotate="varimax") ### Gerando 1 fator 

fit ## print : cargas fatoriais

fit$scores ### scores 

## Comunalidades

fit$communality

## Criar nova base de dados

indice$scores <- fit$scores ### Adicionando uma coluna dos scores em "indice"

indice$NM_MUN <- bd$Codmun6 ### Adiconando uma coluna com o nome dos municípios 

indice$Codmun6 <- bd$Codmun6

indice$pop <- bd$População

## Normalizar os scores da ACP

normalized <- (fit$scores-min(fit$scores))/(max(fit$scores)-min(fit$scores))

indice$Indice_de_proteção <- normalized 

indice[1:6] <- NULL ## excluir variáveis da coluna 1 a 6

## Construindo banco de dados para a análise

library(readxl)

## dados do índice
bd_indice <- indice

## dados de violência contra a mulher (2014)
bd_vio <- read_excel('./dados_hom_2014.xlsx')

## dados atlas
bd_atlas <- read_excel('./banco_atlas.xlsx')

## Dados PIBm 
bd_PIBm <- read_excel('./PIB_dados_mun.xlsx')


# selecionar municípios de interesse (mais de 100 mil hab)
## dos dados do atlas
bd_1 <- subset(bd_vio, bd_vio$Codmun6 %in% bd_indice$Codmun6)
colnames(bd_1)[1]<- 'Codmun6'
## dos dados de violência
bd_2 <- subset(bd_atlas, bd_atlas$Codmun6 %in% bd_indice$Codmun6)
## dos dados do PIBm 
bd_3 <- subset(bd_PIBm, bd_PIBm$Codmun6 %in% bd_indice$Codmun6)

## criar dataframe final com dados para o modelo
df_1 <- merge(bd_1, bd_2)
df_2 <- merge(df_1, bd_3)
df_final <- merge(df_2, bd_indice)

## Método bayesiano empírico 
## Calcular as estimativas

Total_morte_fem <- sum(df_final$Fem)
total_pop <- sum(df_final$pop)
n_barra <- total_pop/284

m <- (Total_morte_fem/total_pop)
df_final$r <- (df_final$Fem/df_final$pop)

s <- sum((df_final$pop*(df_final$r-m)^2)/total_pop)

df_final$c <- (s-(m/n_barra))/((s-(m/n_barra)+m/df_final$pop))

df_final$bay <- df_final$c*df_final$r + (1 + df_final$c)*m 

## Inserir efeitos fixos de Estados 

bd_estados <- read_excel('./dados_estados.xlsx')

df_final_2 <- merge(df_final, bd_estados)

## Modelo de regressão Multivariado com o método bayesiano empírico 

a <- lm(bay ~ Indice_de_proteção + GINI + log(Densidade_pop) + log(PIB_per_capita) + 
          IDHM_2010 + log(PIB) + UF, data = df_final_2)


summary(a)

## Remover coeficiente da variável 'UF' 

a$coefficients <- a$coefficients [1:7]

## Tabela de regressão 

install.packages('stargazer')
library(stargazer)


stargazer(a, type="text", dep.var.labels=c("Taxa bayesiana de Homicídio de mulheres em residência "),
  covariate.labels=c("Índice de Proteção","Índice de Gini","Densidade populacional",
                     "PIBm per capita","IDHm", "PIB"), out="models.txt")



## Gráfico do intervalo de confiança das estimativas 

if(require(dotwhisker) == F) install.packages('dotwhisker'); require(dotwhisker)
if(require(broom) == F) install.packages('broom'); require(broom)

dwplot(a, vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2))


## Análise exploratória dos dados 

## Cálcular o RMSE do modelo 

rmse <- function(x){
  sqrt(mean(x^2))  
}

rmse(a$residuals)

## Análise descritiva 

describe(df_final_2[c("bay", "Indice_de_proteção", "GINI", "Densidade_pop", "PIB_per_capita", 
                      "IDHM_2010", "PIB")])


## Correlação entre as variáveis 

## Criar data frame com as variáveis usadas no modelo

taxa_bayesiana <- df_final_2$bay

Indice_de_Proteção <- df_final_2$Indice_de_proteção

Indice_Gini <- df_final_2$GINI

IDHm <- df_final_2$IDHM_2010

PIBm <- df_final_2$PIB

PIBm_percapita <- df_final_2$PIB_per_capita

Densidade_pop <- df_final_2$Densidade_pop

## Data frame

matriz <- data.frame(taxa_bayesiana, Indice_de_Proteção, Indice_Gini, IDHm, PIBm, PIBm_percapita, Densidade_pop)

## Mudar os nomes das variáveis 

names(matriz)[names(matriz) == "PC1"] <- "Índ Proteção"

names(matriz)[names(matriz) == "taxa_bayesiana"] <- " Feminicídio"

names(matriz)[names(matriz) == "Indice_Gini"] <- "Índ. Gini"

names(matriz)[names(matriz) == "PIBm_percapita"] <- "PIBm pc"

names(matriz)[names(matriz) == "Densidade_pop"] <- "Densidade"

C <- matriz[-c(31,45,74,136,207,263), ]

## Gráfico de correlação 

install.packages("corrplot")
library(corrplot)

M <- cor(C)

head(round(M,2))

corrplot(M, method = 'circle') ## Gráfico


## Índice por Estado (municípios com mais de 100 mil habitantes)

## Box plot 

ggplot (a , aes ( y = Indice_de_proteção, x = UF ))  + 
  geom_boxplot ()


## Histograma 

library(ggplot2)

qplot(df_final_2$UF) ## Municípios contidos em cada estado 

qplot(df_final_2$Indice_de_proteção) ## Os resultados sugerem que a maioria dos municípios apresentam índice 
                                     ## igual a 1 ou próximo de 1    

## Verificar pressupostos e ajuste do modelo 

## Verificar linearidade 

## Gráfico de dispersão 

## Indice_de_proteção 

library(ggplot2)

ggplot(data = df_final_2, aes(y = bay, x = Indice_de_proteção)) + geom_point(color = "blue") + 
  theme_classic() + geom_smooth(method = "lm", color = "black", se = FALSE) + 
  labs(y = "Taxa de Feminicídio", x = "Índice de Proteção") + geom_smooth(method=lm, color='#2C3E50')

## GINI 

ggplot(data = df_final, aes(y = df_final$bay, x = df_final$GINI)) + geom_point(color = "blue") + 
  theme_classic() + geom_smooth(method = "lm", color = "black", se = FALSE) 

## Densidade_pop 

## Função de logarítimo natural 

df_final$densidade_pop <- log(df_final$Densidade_pop) 

ggplot(data = df_final, aes(y = df_final$bay, x = df_final$densidade_pop)) + geom_point(color = "blue") + 
  theme_classic() + geom_smooth(method = "lm", color = "black", se = FALSE) 

## PIB_per_capita

## Função de logarítimo natural 

df_final$pib_percapita <- log(df_final$PIB_per_capita) 

ggplot(data = df_final, aes(y = df_final$bay, x = df_final$pib_percapita)) + geom_point(color = "blue") + 
  theme_classic() + geom_smooth(method = "lm", color = "black", se = FALSE) 

## PIBm 

## Função de logarítimo natural 

df_final$PIBm <- log(df_final$PIB) 

ggplot(data = df_final, aes(y = df_final$bay, x = df_final$PIBm)) + geom_point(color = "blue") + 
  theme_classic() + geom_smooth(method = "lm", color = "black", se = FALSE) 

## IDHm_2010 

ggplot(data = df_final, aes(y = df_final$bay, x = df_final$IDHM_2010)) + geom_point(color = "blue") + 
  theme_classic() + geom_smooth(method = "lm", color = "black", se = FALSE)

## Verificar homocedasticidade 

## Gráfico dos resíduos versus valores ajustados 

ggplot(lm(a)) + geom_point(aes(x = .fitted, y = .resid)) + geom_abline(slope = 0)

## Verificar normalidade 

ggplot(data = a, aes(sample = a$residuals)) + 
  stat_qq(color = 'blue') + stat_qq_line(color = 'red', lty = 2)

## Verificar a média dos erros 

mean(a$residuals)

## Verificar multicolinariedade 

library(car)

vif(a) ## Se for maior que 10 a multicolineridade é fortemente sugerida 

## Uma diretriz geral sugere que um VIF maior que 5 ou 10 é grande, indicando que o modelo 
## tem probelemas para estimar os coeficientes. No entanto isso não prejudica a qualidade das previsões

## Salvar base de dados final em XLSX

install.packages('xlsx')
library(xlsx)

write.xlsx(df_final_2, 'df_final_2.xlsx')
