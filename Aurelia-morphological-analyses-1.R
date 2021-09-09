######################### Aurelia morphological analyses 1 ##########################
#                       R Code - Author: Jonathan W. Lawley                         #
#                           Used in: Lawley et al., 2021                            #
#   The importance of molecular characters when morphological variability hinders   # 
#                   diagnosability: systematics of the moon jellyfish               #
#                           genus Aurelia (Cnidaria: Scyphozoa)                     #
#                       PeerJ https://doi.org/10.7717/peerj.11954                   #
#                                                                                   #
#     Used for MDS analysis in Fig. 3 and Mantel test, including all specimens      #
#                           WITH estimation of missing data                         #
#                                                                                   #
#               Contact: jonathan.wanderleylawley@griffithuni.edu.au                #

library(vegan)
library(ade4)

#Data prep, log transform, slope calculations and isometry test
aurelia.raw=read.csv("aurelia-all.csv", sep=",", as.is=T) #raw data inserted
str(aurelia.raw)

aurelia.raw1=aurelia.raw[,-c(8,12,16)] #removes categorical features f8, f12 and f21

reference=aurelia.raw1$f1 #setting f1 as reference variable
response=aurelia.raw1[,3:36] #setting the rest as response variables
slope=rep(NA, ncol(response)) #empty vector for slope values
result=rep(NA, ncol(response)) #empty vector for isometry test results (Allometry, Isometry or lm not significant)

#Loop to calculate slopes of log-transformed variables and isometry test
for (i in 1:ncol(response)){
  response.lm=lm(log(response[,i]+1)~log(reference+1))
  f=(summary(response.lm))$fstatistic #creates the object "f" with the F-statistic values of the linear model.
  p=pf(f[1],f[2],f[3],lower.tail=F) #creates the object "p" with the p-value from the F-statistic values above.
  
  if (p > 0.05) #verifies if "p" is bigger than "alpha".
  {
    slope[i]=NA #if yes, attribute "lm not significant" to element k in "result".
    result[i]="lm not significant"
  }
  else #if not, linear model is significant, so proceeds with calculations.
  {
    slope[i]=(summary(response.lm))$coefficients[2]
    slope.se=(summary(response.lm))$coefficients[4]
    tvalue=abs((slope[i]-1)/slope.se)
    pvalue=(1-pt(tvalue,nrow(na.omit(cbind(reference,response[,i])))-2))*2
    if (pvalue > 0.05) #verifies if "pvalue" is bigger than "alpha".
    {
      result[i]="Isometry" #if yes, the slope from response.lm is not significantly different from the expected slope, so it attributes "Isometry" to element k in "result".
    }
    else #verifies if "pvalue" is not bigger than alpha.
    { 
      result[i]="Allometry" #if it is not bigger, the slope from response.m is significantly different from the expected.slope, so it attributes "Allometry" to element k in "result".
    }
  }
}
slope #resulting slopes (NA where lm not significant)
result #result of isometry test
slope.true=slope[-c(2,10:11,13,15:18,23,27:28,34)] #removing NAs (features in which lm was not significant or that are nonetheless unusable for transformations) as well as mostly invariable or biased (f3a, f19, f33) variables.

#Exploratory checking of response variables
response2=response[,c(11,13,15:18,27:28,34)] # with variables that had no significanct lm
response2.test=response[,11] #change column here for each of the above to test
plot(log(response2.test+1)~log(reference+1))
summary(lm(log(response2.test+1)~log(reference+1)))
plot(lm(log(response2.test+1)~log(reference+1)))

#Removing mostly invariable or biased (f3a, f19, f33) and uninformative (lm not significant, from above) variables - f3a, f19, f20, f23, f25a, f26, f27, f28a, f33, f37, f38, f44
aurelia.raw2=aurelia.raw1[,-c(4,12:13,15,17:20,25,29:30,36)]

#Estimations of missing values based on the genus-level slope and "country" group (morphogroup) averages (if possible/plausible)
aurelia.raw3=aurelia.raw2[,-c(25:28)] #removing data that will not be transformed, including locality info (some missing data there from aquariums)

groups=split(aurelia.raw3[,2:24], aurelia.raw3[,25]) #defining groups based on "countries" (morphogroups)

for (i in 1:length(groups)){ #loop to enter each country group
  groups.temp=groups[[i]] #define the group entered in the loop as group.temp
  temp=which(is.na(groups.temp), arr.ind=T) #what is the location of NAs in the group
  test=which(table(temp[,2])==nrow(groups.temp)) #which of the columns have only NAs (no estimation possible)
  test2=as.numeric(names(test)) #gets the number of these columns
  if (length(temp) == 0) next #no NAs in the group to be estimated, go to the next group
  if (length(test2) > 0) { #new.temp has the columns with NAs that can be estimated
    new.temp=temp[!(temp[,2] %in% test2),]
  }
  if (length(test2) == 0) {
    new.temp=temp
  }
  if (length(new.temp) == 0) next #NAs in the group but no other values to estimate it
  if (length(new.temp) == 2) { #issue to correct when only one NA value to be estimated, new.temp becomes integer, so corrected below
    new.temp=t(as.matrix(new.temp))
  }
  #kept only features that have some value in the group for estimations
  for (j in 1:nrow(new.temp)){ 
    temp2=apply(na.omit(cbind(groups.temp[,1],groups.temp[,new.temp[j,2]])), 2, mean) #uses means of response and reference variables, excluding NA rows, as reference for estimations
    groups.temp[new.temp[j,1],new.temp[j,2]]=temp2[2]*((groups.temp[new.temp[j,1],1]/temp2[1])^slope.true[new.temp[j,2]-1]) #formula to calculate based on Lleonart et al. (2000)
  }
  groups[[i]]=groups.temp #group with estimations sent back to the groups object
}

groups
aurelia.estimated=unsplit(groups, aurelia.raw3[,25]) #grouping again based on the morphogroup split
aurelia.estimated=cbind(aurelia.estimated,aurelia.raw2[,27]) #adding countries/regions (not morphogroups)
colnames(aurelia.estimated)[24]<-"Region"
rownames(aurelia.estimated)<-aurelia.raw3$Specimen
aurelia.estimated
#Aurelias that were divided in a country/region in which estimations were not possible, were associated with other specimens in Morphogroups, based on similarity distribution and/or morphology

#Size corrections based on Lleonart et al. (2000) and full standardization
aurelia.corrected=matrix(data=NA, nrow=length(rownames(aurelia.estimated)), ncol=length(colnames(aurelia.estimated))-2)
for (i in 1:length(slope.true)){
  aurelia.corrected[,i]=aurelia.estimated[,i+1]*((mean(aurelia.estimated[,1])/aurelia.estimated[,1])^slope.true[i]) #mean f1 is 10.14988
}
aurelia.corrected=as.data.frame(aurelia.corrected)

aurelia.corrected.01=matrix(NA,173,22) #scale full standardization
for (i in 1:ncol(aurelia.corrected)) {
  for (k in 1:nrow(aurelia.corrected)) {
    aurelia.corrected.01[k,i]=(aurelia.corrected[k,i]-min(aurelia.corrected[,i]))/(max(aurelia.corrected[,i])-min(aurelia.corrected[,i]))
  }
}
colnames(aurelia.corrected.01)=colnames(aurelia.estimated[,-c(1,24)])
rownames(aurelia.corrected.01)=rownames(aurelia.estimated)

#Multivariate analysis
par(cex=0.5)
mds=cmdscale((vegdist(aurelia.corrected.01, "gower")))
mds=mds*-1 #not necessary, but used to display MDS differently
features.wa<-wascores(mds, aurelia.corrected.01, expand=T) #with weighted averages of morphological features
plot(rbind(mds,features.wa), pch = 15, type="n", xlab = "Dim1", ylab = "Dim2")
text(features.wa[,1], features.wa[,2], rownames(features.wa), cex=1.3, col="red")
text(mds[,1], mds[,2], aurelia.estimated[,24], cex=1)

#Mantel test (correlation between two distance matrices)
aurelia.gd=na.omit(aurelia.raw2[,c(1,25,26)])
colnames(aurelia.gd)=c("Specimen", "Lat","Lon")
gd.dist=dist(cbind(aurelia.gd$Lon, aurelia.gd$Lat))
md.matrix=as.matrix(vegdist(aurelia.corrected.01, "gower"))
md.matrix=md.matrix[-c(134:135,161:168),-c(134:135,161:168)] #removing aquarium specimens
md.dist=as.dist(md.matrix)
mantel.rtest(md.dist, gd.dist, nrepet=9999)
