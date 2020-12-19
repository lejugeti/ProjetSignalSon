library("ggplot2")
library("rstatix")

df = read.csv("../data_eeg.csv")


e1_DFA = df[(df$electrode=="e1" & df$method=="DFA"),]
e2_DFA = df[(df$electrode=="e2" & df$method=="DFA"),]
e1_DMA = df[(df$electrode=="e1" & df$method=="DMA"),]
e2_DMA = df[(df$electrode=="e2" & df$method=="DMA"),]

ggplot(df, aes(x=etat, y=alpha, fill=method)) + 
  geom_boxplot()

ggplot(e1, aes(x=etat, y=alpha, fill=method)) + 
  geom_boxplot() +
  ggtitle("Boxplot Electrode 1")

ggplot(e2, aes(x=etat, y=alpha, fill=method)) + 
  geom_boxplot() +
  ggtitle("Boxplot Electrode 2")   

res1 = anova_test(e1, alpha~method * etat)
res2 = anova_test(e2, alpha~method * etat)

res1
res2


# Analyses attendues par le prof ------------------------------------------

# on fait des anova pour les alpha sur une électrode pour une méthode donnée 
# et en fonction du type d'état

aov_e1_DFA = anova_test(e1_DFA, alpha~etat)
aov_e2_DFA = anova_test(e2_DFA, alpha~etat)
aov_e1_DMA = anova_test(e1_DMA, alpha~etat)
aov_e2_DMA = anova_test(e2_DMA, alpha~etat)

aov_e1_DFA
aov_e2_DFA
aov_e1_DMA
aov_e2_DMA
