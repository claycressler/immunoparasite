x = read.csv("Dimensionless_results_across_immune_parasite_initials_and_total_immune.csv", header=FALSE)

y = as.data.frame(matrix(unlist(x), ncol=4, byrow=TRUE))
colnames(y) <- c("totalI","initialP","initialTh1bias","equilP")
unique(round(y$equilP,3))
y$equilP <- round(y$equil,3)


mybreaks = c(sort(unique(y$equilP)),1)
library(patchwork)

a1 = ggplot(subset(y,totalI==1), aes(initialTh1bias,initialP,z=equilP)) + 
  geom_contour_filled(breaks=mybreaks) + 
  scale_fill_discrete(labels=c("0.0","0.02","0.87","0.88"))

a2 = ggplot(subset(y,totalI==0.5), aes(initialTh1bias,initialP,z=equilP)) + 
  geom_contour_filled(breaks=mybreaks) + 
  scale_fill_discrete(labels=c("0.0","0.02","0.87","0.88"))

a3 = ggplot(subset(y,totalI>1.49 & totalI < 1.51), aes(initialTh1bias,initialP,z=equilP)) + 
  geom_contour_filled(breaks=mybreaks) + 
  scale_fill_discrete(labels=c("0.0","0.02","0.87","0.88"))

a4 = ggplot(subset(y,totalI==2), aes(initialTh1bias,initialP,z=equilP)) + 
  geom_contour_filled(breaks=mybreaks) + 
  scale_fill_discrete(labels=c("0.0","0.02","0.87","0.88"))

ggplot(y, aes(initialTh1bias,initialP,z=equilP)) + 
  geom_contour_filled(breaks=mybreaks) + 
  scale_fill_discrete(labels=c("0.0","0.02","0.87","0.88")) + 
  facet_wrap(~totalI) + theme_bw()


