############################ Aurelia coerulea comparison ############################
#                       R Code - Author: Jonathan W. Lawley                         #
#                           Used in: Lawley et al., 2021                            #
#   The importance of molecular characters when morphological variability hinders   # 
#                   diagnosability: systematics of the moon jellyfish               #
#                           genus Aurelia (Cnidaria: Scyphozoa)                     #
#                       PeerJ https://doi.org/10.7717/peerj.11954                   #
#                                                                                   #
#                        Used for Welch's t-test in Fig. 8                          #
#                                                                                   #
#               Contact: jonathan.wanderleylawley@griffithuni.edu.au                #

coerulea=read.csv("coerulea.csv", sep=",", as.is=T)
compare=t(cbind(coerulea$mean.med, coerulea$mean.lab))
stdev=t(cbind(coerulea$sd.med, coerulea$sd.lab))
sampling=t(cbind(coerulea$n.med, coerulea$n.lab))
par(mar=c(5,5,2,2), cex=1.5)
bp=barplot(compare, beside=T, xpd=FALSE, ylim = c(0,70), ylab="% of bell diameter (f1)", xlab = "Morphological features", legend=c("Mediterranean", "North Sea (cultured in lab)"))
box()
arrows(x0=bp, y0=compare-stdev, y1=compare+stdev, code=3, angle=90, length=.07)
axis(1, at=c(2,5,8,11,14,17,20,23), labels=coerulea$feature)

#Welch's two samples t-test, independent samples, equal or unequal variances
result=matrix(NA, 8, 1)
for (i in 1:8) {
  temp1=rnorm(sampling[2,i], compare[2,i], stdev[2,i])
  temp2=rnorm(sampling[1,i], compare[1,i], stdev[1,i])
  value=t.test(temp1, temp2)
  result[i,1]=value$p.value
}
result



