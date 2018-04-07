##############################################################################
# Purpose:                                                                 
# Author: R Epi Team                          
# Input:                                                                     
# Output:                                                               
# Remarks:                                                                    
###############################################################################

######################################## 
#########      Imports      ############
########################################
library(tidyverse)
# https://stackoverflow.com/questions/37755037/how-to-add-code-folding-to-output-chunks-in-rmarkdown-html-documents/37839683#37839683



######################################## 
#########      Simulate      ############
########################################

n <- 1000 # sample size

df <- data.frame(dz=rep(NA, n),
                 exp=c(rep(1, n/2), rep(0,n/2)),
                 c=rbinom(n,1,0.5)

)




p1 <- 0.5 # p(D = 1 | E = 1)
p0 <- 0.1 # p(D = 1 | E = 0)
e <- 0.3 # this is the "error" introduced by confounding or the increased association given differning levels of the confounder

df$dz[df$exp == 1 & df$c == 1] <- rbinom(sum(df$exp == 1 & df$c == 1),1,p1 + e)
df$dz[df$exp == 0 & df$c == 1] <- rbinom(sum(df$exp == 0 & df$c == 1),1,p0 + e)
df$dz[df$exp == 1 & df$c == 0] <- rbinom(sum(df$exp == 1 & df$c == 0),1,p1)
df$dz[df$exp == 0 & df$c == 0] <- rbinom(sum(df$exp == 0 & df$c == 0),1,p0)

mean(df$dz)
mean(df$exp)
mean(df$c)
mean(df$dz[df$exp==1])
mean(df$dz[df$exp==0])

fisher.test(table(df$exp,df$dz))
fisher.test(table(df$c,df$dz))
fisher.test(table(df$c,df$exp))

# crude 
prop.table(xtabs(~df$exp + df$dz), margin = 1)
mean(df$dz[df$exp==1])
mean(df$dz[df$exp==0])
prop.table(xtabs(~df$c + df$dz), margin = 1)
mean(df$dz[df$c==1])
mean(df$dz[df$c==0])
cruderd = prop.table(xtabs(~df$exp + df$dz), margin = 1)[2,2] - prop.table(xtabs(~df$exp + df$dz), margin = 1)[1,2]
# stratified 
prop.table(xtabs(~df$exp + df$dz + df$c))



######## this totally gets screwed up when we condition



#### try and get exposure 

df <- data.frame(dz=rep(NA, n),
                 exp=c(rep(1, n/2), rep(0,n/2)),
                 c=rep(NA, n)
                 
)

p1 <- 0.5 # p(D = 1 | E = 1)
p0 <- 0.1 # p(D = 1 | E = 0)
e <- 0.3 # this is the "error" introduced by confounding or the increased association given differning levels of the confounder


df$c[df$exp==1] <- rbinom(sum(df$exp==1),1,0.8)
df$c[df$exp==0] <- rbinom(sum(df$exp==1),1,0.3)


df$dz[df$exp == 1 & df$c == 1] <- rbinom(sum(df$exp == 1 & df$c == 1),1,p1 + e)
df$dz[df$exp == 0 & df$c == 1] <- rbinom(sum(df$exp == 0 & df$c == 1),1,p0 + e)
df$dz[df$exp == 1 & df$c == 0] <- rbinom(sum(df$exp == 1 & df$c == 0),1,p1)
df$dz[df$exp == 0 & df$c == 0] <- rbinom(sum(df$exp == 0 & df$c == 0),1,p0)

mean(df$dz)
mean(df$exp)
mean(df$c)
mean(df$dz[df$exp==1])
mean(df$dz[df$exp==0])

fisher.test(table(df$exp,df$dz))
fisher.test(table(df$c,df$dz))
fisher.test(table(df$c,df$exp))


prop.table(xtabs(~df$exp + df$dz), margin = 1)
cruderd = prop.table(xtabs(~df$exp + df$dz), margin = 1)[2,2] - prop.table(xtabs(~df$exp + df$dz), margin = 1)[1,2]
# stratified 
prop.table(xtabs(~df$exp + df$dz + df$c), margin = 1)
stratrdz0 = prop.table(xtabs(~df$exp + df$dz + df$c)[,,1], margin = 1)[2,2] - prop.table(xtabs(~df$exp + df$dz + df$c)[,,1], margin = 1)[1,2]
stratrdz1 = prop.table(xtabs(~df$exp + df$dz + df$c)[,,2], margin = 1)[2,2] - prop.table(xtabs(~df$exp + df$dz + df$c)[,,2], margin = 1)[1,2]
realrd <- p1-p0








