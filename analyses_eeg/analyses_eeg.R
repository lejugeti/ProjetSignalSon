library("ggplot2")
library("rstatix")

df = read.csv("../data_eeg.csv")


e1 = df[df$electrode=="A",]
e2 = df[df$electrode=="B",]

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
