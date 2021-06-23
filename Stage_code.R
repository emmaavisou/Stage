# Ouvrir la base de données pour les corrélations
library(readr)
corr <- read.delim("Desktop/test_corr.csv", sep = ";")

corr$GMT_Vaccin <- as.numeric(as.character(corr$GMT_Vaccin))
head(corr)
dat$Nom_vaccin <- as.factor(dat$Nom_vaccin)
head(corr)

#####################################################################
# Analyses de corrélations
####################################################################

# Corrélation de rang de spearman entre efficacité et GMT neutralisant
library("ggpubr")
p <- ggscatter(corr, x = "GMT_Vaccin", y = "Efficacite_.", repel = TRUE,
          color = "black", size = "Cas", label = "Nom_vaccin", label.rectangle = TRUE,
          add = "reg.line",  # Add regressin line
          add.params = list(color = "orange", linetype="dashed"), 
          conf.int = FALSE, # Add confidence interval
          cor.coef = FALSE, # Add correlation coefficient. see ?stat_cor
          cor.method = "spearman"
)
p+ scale_y_continuous(breaks=c(20,50,60,70,80,90,95,99)) + coord_trans(x="log2") +scale_x_continuous(breaks=c(10, 50, 100, 500,1000,4000))
ggpar(p, xscale= "log2", xlim=c(10,5000), ylim=c(40,99), xlab = "SARS-Cov-2 neutralization (GMT)", ylab="Efficacité des vaccins" )+ scale_y_continuous(breaks=c(0,20,50,60,70,80,90,95,99)) + scale_x_log10(breaks=c(10,50,100,500,1000,5000))+ stat_cor(method = "spearman", label.x = 1, label.y = 99)


# Corrélation de rang de spearman entre efficacité et GMT IgG ELISA
p3 <- ggscatter(corr, x = "GMT_Elisa", y = "Efficacite_.",repel = TRUE,
               color = "black", size = "Cas", label = "Nom_vaccin", label.rectangle = TRUE,
               add = "reg.line",  # Add regressin line
               add.params = list(color = "red", linetype="dashed"), # Customize reg. line
               conf.int = FALSE, # Add confidence interval
               cor.coef = FALSE, # Add correlation coefficient. see ?stat_cor
               cor.method = "spearman"
)
p3+ scale_y_continuous(breaks=c(20,50,60,70,80,90,95,99)) + coord_trans(x="log2") +scale_x_continuous(breaks=c(10, 50, 100, 500,1000,4000))
library(scales)
ggpar(p3, ylim=c(40,99), xlab = "SARS-Cov-2 Spike IgG ELISA", ylab="Efficacité des vaccins" )+ scale_y_continuous(breaks=c(0,20,50,60,70,80,90,95,99)) + stat_cor(method = "spearman", label.x = 3, label.y = 99) + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))


# Corrélation de rang de spearman entre efficacité et ratio GMT neutralisant/Convalescents
p1 <- ggscatter(corr, x = "Vaccin.HCS", y = "Efficacite_.",
               color = "black", size = "Cas", label = "Nom_vaccin", repel = TRUE,label.rectangle = TRUE,# Points color, shape and size
               add = "reg.line",  # Add regressin line
               add.params = list(color = "orange", linetype="dashed"), # Customize reg. line
               conf.int = FALSE, # Add confidence interval
               cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
               cor.method = "spearman"
)
p1
ggpar(p1, xlim=c(0.1,5), ylim=c(40,99), xlab = "SARS-Cov-2 neutralisation (ratio Vaccinés/Convalescents)", ylab="Efficacité des vaccins (%)" )+ scale_y_log10(breaks=c(0,20,50,60,70,80,90,95,99)) + scale_x_log10(breaks=c(0.1,0.2,0.5,0.8,1,1.5,2,3,5))


# Corrélation de rang de spearman entre efficacité et ratio GMT IgG ELISA/Convalescents
p2 <- ggscatter(corr, x = "Vaccin.HCSElisa", y = "Efficacite_.",
                color = "black", size = "Cas", label = "Nom_vaccin", repel = TRUE,label.rectangle = TRUE,# Points color, shape and size
                add = "reg.line",  # Add regressin line
                add.params = list(color = "red", linetype="dashed"), # Customize reg. line
                conf.int = FALSE, # Add confidence interval
                cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                cor.method = "spearman"
)
p2
ggpar(p2, xlim=c(0.1,20), ylim=c(40,99), xlab = "SARS-Cov-2 neutralisation (ratio Vaccinés/Convalescents)", ylab="Efficacité des vaccins (%)" )+ scale_y_continuous(breaks=c(0,20,50,60,70,80,90,95,99)) + scale_x_log10(breaks=c(0.1,0.2,0.5,0.8,1,1.5,2,3,5, 15, 20)) + stat_cor(method = "spearman", label.x=3, label.y = 97) 

# Etude des corrélations de Spearman et variance expliqué (carré du coefficient de pearson) sans graphique
# Ratio neutralisant/Convalesents
res <- cor.test(corr$Vaccin.HCS, corr$Efficacite_., 
                method = "spearman")
res
res4 <- cor.test(corr$Vaccin.HCS, corr$Efficacite_., 
                method = "pearson")
((res4$estimate)^2)*100 # Variance expliquée en pourcentage

# Ratio IgG Elisa/Convalesents
cor.test(corr$Vaccin.HCSElisa, corr$Efficacite_., 
                method = "spearman")
res6 <- cor.test(corr$Vaccin.HCSElisa, corr$Efficacite_., 
                 method = "pearson")
((res6$estimate)^2)*100 # Variance expliquée en pourcentage

# Efficacité et GMT ELISA
cor.test( ~ GMT_Elisa + Efficacite_., 
          data=corr,
          method = "spearman",
          conf.level = 0.95)
res5 <- cor.test(corr$GMT_Elisa, corr$Efficacite_., 
                method = "pearson")
((res5$estimate)^2)*100

# Efficacité et GMT neutralisant
cor.test( ~ GMT_Vaccin + Efficacite_., 
          data=corr,
          method = "spearman",
          conf.level = 0.95)
res7 <- cor.test(corr$GMT_Vaccin, corr$Efficacite_., 
                 method = "pearson")
((res7$estimate)^2)*100


############################################
# Meta-analyse avec le package métafor
###########################################
library(metafor)
dat  <-  read.delim("Desktop/Metadonnees.csv", sep = ";")
dat1 <-  escalc ( measure = "RR" , ai = ai, bi=bi, ci = ci, di=di, data = dat)
dat1
print(dat1, row.names = FALSE)

res <- rma(ai = ai, bi = bi, ci = ci, di = di, data = dat1, measure = "RR")
res
# Resultat du RR moyen entre les études 0.18
predict(res, transf=exp)


?forest

# Creation du forest plot 
forest(res, slab = paste(dat1$Nom_vaccin), xlim = c(-16, 6), at = log(c(0.05, 0.25, 1, 4)), atransf = exp, ilab = cbind(dat1$ai, dat1$bi, dat1$ci, dat1$di),ilab.xpos = c(-9.5, -8, -6, -4.5), cex = 0.75, header="Vaccins")
op <- par(cex = 0.5, font = 2) 
text(c(-9.5, -8, -6, -4.5), 15, c("Covid+", "Covid-", "Covid+", "Covid-")) 
text(c(-8.75, -5.25), 16, c("Vaccinés", "Convalescents")) 
text(-16, 15, "Vaccins", pos = 4) 
text(6, 15, "Risque Relatif [95% CI]", pos = 2) 
par(op)

### Ajout du texte Q-value, dfs, p-value, and I^2 statistique au forest plot
text(-16, -1, pos=4, cex=0.75, bquote(paste("RE Model (Q = ",
                                            .(formatC(res$QE, digits=2, format="f")), ", df = ", .(res$k - res$p),
                                            ", p = ", .(formatC(res$QEp, digits=2, format="f")), "; ", I^2, " = ",
                                            .(formatC(res$I2, digits=1, format="f")), "%)")))



# Modèle à effet fixe
res1 <- rma(yi,vi, data=dat1, method = "FE")
res1
predict(res1, transf=exp)

#Modèle à effet aléatoire
res10 <- rma(yi,vi, data=dat1, method = "REML")
res10
predict(res10, transf=exp)

#Modèle à effet mixte en rajoutant la variable du ratio comme variable modératrice.
res2 <-  rma (yi, vi, mods = ~ Vaccin.HCS, data = dat1, method = "ML" ) 
res2

# Prédiction basé sur modèle à effet mixte
predict(res2, newmods = cbind(seq(from = 0.7, to = 16, by = 1.5)),transf = exp, addx = TRUE)
# Estimation pour vaccin Curevac : Ratio = 0.70 donc efficacité = 65% (1-0.35)*100
preds <- predict(res2, newmods = cbind(0.2:10), transf = exp) 
wi <- 1/sqrt(dat1$vi)
size <- 0.5 + 3 * (wi - min(wi))/(max(wi) - min(wi))
?plot
plot(dat1$Vaccin.HCS, exp(dat1$yi), pch = 19, xlab = "Ratio GMT Vaccinés/Convalescents", ylab = "Risque relatif", size=dat1$Cas,las = 1, bty = "l", log = "y") 
lines(0.2:10, preds$pred)
lines(0.2:10, preds$ci.lb, lty = "dashed")
lines(0.2:10, preds$ci.ub, lty = "dashed")
abline(h = 1, lty = "dotted")


#######################################################################
# Nombre de cas nécessaire pour efficacité
######################################################################

#Calcul du risque relatif
RR <- (dat$ai/(dat$ai+dat$bi))/(dat$ci/(dat$ci+dat$di))
RR

# Intervalle de confiance du RR
SE <- sqrt(((dat$bi/dat$ai)/(dat$ai+dat$bi))+((dat$di/dat$ci)/(dat$ci+dat$di)))
ICsup <- exp(log(RR)+1.96*SE)
ICsup
ICinf <- exp(log(RR)-1.96*SE)
ICinf

# Relative Risk reduction ou efficacy
RRR <- 1-RR
RRR

#Intervalle de confiance RRR
RRRICsup <- 1-ICsup
RRRICsup

RRRICinf <- 1-ICinf
RRRICinf
# Absolute Risk reduction
ARR <- (dat$ci/(dat$ci+dat$di)-dat$ai/(dat$ai+dat$bi))
ARR

# Intervalle de confiance ARR
SE2 <- sqrt(((dat$ai/(dat$ai+dat$bi))*(1-(dat$ai/(dat$ai+dat$bi))))/(dat$ai+dat$bi)+((dat$ci/(dat$ci+dat$di))*(1-(dat$ci/(dat$ci+dat$di))))/(dat$ci+dat$di))

ARRICsup <- ARR+1.96*SE2
ARRICsup

ARRICinf <- ARR-1.96*SE2
ARRICinf

# Nombre neccessaire a vacciné 
NNV <- 1/ARR
NNV

#Intervalle confiance NNV
NNVICsup <- 1/(ARRICsup)
NNVICsup

NNVICinf <- 1/(ARRICinf)
NNVICinf

library(ggplot2)
d=data.frame(Vaccins=c("Pfizer/BioNTech","Moderna","Gamaleya","Oxford/AstraZeneca","Sinovac","Novarax","Janssen"), mean=c(NNV), lower=c(NNVICinf), upper=c(NNVICsup))
d
# Test 1 graphe
ggplot() + geom_pointrange(data=d, mapping=aes(x=Vaccins, y=mean, ymin=upper, ymax=lower), size=1, color="red", fill="white", shape=22) 

# Graphe final
library(ggplot2)
library(dplyr)
library(patchwork) 
library(hrbrthemes)
coeff=0.2
Color1 <- "#0099CC"
Color2 <- "#F6B113"

p <- ggplot(df, aes(x=Vaccins),label=labelNNV) +
  
  geom_pointrange( aes(y=NNV,ymin = upperNNV, ymax = lowerNNV), size=0.5, color=Color1) + 
  geom_pointrange( aes(y=RRR/coeff, ymin = lowerRRR/coeff, ymax = upperRRR/coeff), size=0.5, color=Color2)+
  
  scale_y_continuous(
    
    # Features of the first axis
    name = "NNV (nombre de personnes)",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name="Efficacité ou RRR (%)")
  ) + 
  theme_classic() +
  
  theme(
    axis.title.y = element_text(color = Color1, size=13),
    axis.title.y.right = element_text(color = Color2, size=13)
  )
p

#######################################################################
# Severe Cases
#######################################################################

library(readr)
sev <- read.delim("Desktop/severe.csv", sep = ";")
sev1 <-  escalc ( measure = "RR" , ai = ai, bi=bi, ci = ci, di=di, data = sev)
resu <- rma(ai = ai, bi = bi, ci = ci, di = di, data = sev1, measure = "RR")
resu
predict(resu, transf=exp)


forest(resu, slab = paste(sev1$Nom_vaccin), xlim = c(-16, 6), at = log(c(0.1, 0.5, 1, 4)), atransf = exp, ilab = cbind(sev1$ai, sev1$bi, sev1$ci, sev1$di),ilab.xpos = c(-9.5, -8, -6, -4.5), cex = 0.75, header="Vaccins")
op <- par(cex = 0.5, font = 2) 
text(c(-9.5, -8, -6, -4.5), 15, c("Covid+", "Covid-", "Covid+", "Covid-")) 
text(c(-8.75, -5.25), 16, c("Vaccinés", "Convalescents")) 
text(-16, 15, "Vaccins", pos = 4) 
text(6, 15, "Risque Relatif [95% CI]", pos = 2) 
par(op)

resul <-  rma (yi, vi, mods = ~ Vaccin.HCS , data = sev1, method = "REML" ) 
resul

preds <- predict(resul, newmods = cbind(0.2:10), transf = exp) 
wi <- 1/sqrt(dat1$vi)
size <- 0.5 + 3 * (wi - min(wi))/(max(wi) - min(wi))
plot(sev1$Vaccin.HCS, exp(sev1$yi), pch = 19, cex = size, xlab = "Ratio GMT Vaccinés/Convalescents", ylab = "Risque relatif", las = 1, bty = "l", log = "y") 
lines(0.2:10, preds$pred)
lines(0.2:10, preds$ci.lb, lty = "dashed")
lines(0.2:10, preds$ci.ub, lty = "dashed")
abline(h = 1, lty = "dotted")


# Graphe corrélation ELISA
p2sev <- ggscatter(sev, x = "GMT_Elisa", y = "Efficacite",
                color = "black", label = "Nom_vaccin", repel = TRUE,label.rectangle = TRUE,# Points color, shape and size
                add = "reg.line",  # Add regressin line
                add.params = list(color = "red", linetype="dashed"), # Customize reg. line
                conf.int = FALSE, # Add confidence interval
                cor.coef = FALSE, # Add correlation coefficient. see ?stat_cor
                cor.method = "spearman"
)

p2sev

ggpar(p2sev, ylim=c(80,100), xlab = "SARS-Cov-2 Spike IgG ELISA", ylab="Efficacité des vaccins" )+ scale_y_continuous(breaks=c(80,90,95,100)) + stat_cor(method = "spearman", label.x = 3, label.y = 99) + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))

res <- cor.test(sev$Vaccin.HCS, sev$Efficacite, 
                method = "spearman")
res

# Graphe corrélation GMT?HCSElisa

p2 <- ggscatter(sev, x = "Vaccin.HCSElisa", y = "Efficacite",
                color = "black", label = "Nom_vaccin", repel = TRUE,label.rectangle = TRUE,# Points color, shape and size
                add = "reg.line",  # Add regressin line
                add.params = list(color = "red", linetype="dashed"), # Customize reg. line
                conf.int = FALSE, # Add confidence interval
                cor.coef = FALSE, # Add correlation coefficient. see ?stat_cor
                cor.method = "spearman"
)

p2
ggpar(p2, xlim=c(0.5,20), ylim=c(75,100), xlab = "SARS-Cov-2 neutralisation (ratio Vaccinés/Convalescents)", ylab="Efficacité des vaccins (%)" )+ scale_y_continuous(breaks=c(75,80,90,95,100)) + scale_x_log10(breaks=c(0.1,0.2,0.5,0.8,1,1.5,2,3,5, 15, 20)) + stat_cor(method = "spearman", label.x=0.001, label.y = 97) 


# Graphe corrélation GMT_Neutralizing
p <- ggscatter(sev, x = "GMT_Vaccin", y = "Efficacite", repel = TRUE,
               color = "black", label = "Nom_vaccin", label.rectangle = TRUE,# Points color, shape and size
               add = "reg.line",  # Add regressin line
               add.params = list(color = "orange", linetype="dashed"), # Customize reg. line
               conf.int = FALSE, # Add confidence interval
               cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
               cor.method = "spearman",
               cor.coef.coord = c(1, 98)
)
p

# Graphe corrélation GMT_Neutralizing HCS
p <- ggscatter(sev, x = "Vaccin.HCS", y = "Efficacite", repel = TRUE,
               color = "black", label = "Nom_vaccin", label.rectangle = TRUE,# Points color, shape and size
               add = "reg.line",  # Add regressin line
               add.params = list(color = "orange", linetype="dashed"), # Customize reg. line
               conf.int = FALSE, # Add confidence interval
               cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
               cor.method = "spearman",
               cor.coef.coord = c(1, 98)
)
p

###############################
# Corrélation NNV / Efficacité
###############################
p5 <- ggscatter(df, x = "VaccinHCS", y = "NNV",
                color = "black", label = "Vaccins", repel = TRUE,label.rectangle = TRUE,# Points color, shape and size
                add = "reg.line",  # Add regressin line
                add.params = list(color = "#E55E5B", linetype="dashed"), # Customize reg. line
                conf.int = FALSE, # Add confidence interval
                cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                cor.method = "spearman"
)

p5
ggpar(p5, xlim=c(0.1,5),ylim=c(50,350), ylab = "SARS-Cov-2 NNV", xlab="Ratio GMT Vaccinés/Convalescents" )+ scale_y_log10(breaks=c(0,50,100,150,200,250,300,350)) + scale_x_log10(breaks=c(0.1,0.2,0.5,0.8,1,1.5,2,3,5))

res23 <- cor.test(df$NNV, df$VaccinHCS, 
                method = "spearman")
res23

