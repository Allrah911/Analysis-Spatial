install.packages("ape") #Analyses of Phylogenetics and Evolution
library(ape)

Weight=matrix(0,nrow=5,ncol=5)
Weight[1,c(2,3)]=1
Weight[2,c(1,3,4)]=1
Weight[3,c(1,2,4)]=1
Weight[4,c(2,3,5)]=1
Weight[5,c(4)]=1

unemploy=c(8,6,6,3,2)
I_employ=Moran.I(unemploy, Weight) #ada dalam package ape
z=c(unemploy-mean(unemploy))
zt=t(z)

MI=5/12*(t(z)%*%Weight%*%z)/(t(z)%*%z)

#menggunakan bobot terstandarisasi
Weight_sd=Weight/rowSums(Weight)

MI_sd=(t(z)%*%Weight_sd%*%z)/(t(z)%*%z)
I_employ_sd=Moran.I(unemploy, Weight_sd)

library(readxl)
artikel=read_excel("E:/Dokumen Penting/Materi Kuliah/Analisis Spatial/TUGAS/Artikel1.xlsx", 
                   col_types = c("text", "numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric")) 

summary(artikel)
View(artikel)
attach(artikel)
regOLS=lm(Kemiskinan~+APK+APM+APS)
summary(regOLS)

library(lmtest)
dwtest(regOLS)
bptest(regOLS)

library(rgdal)
library(spdep)
library(tmap)
library(rgeos)

ptjatim=readOGR(dsn=".",layer="jatim") #dari package rgdal
plot(ptjatim)

class(ptjatim)
View(ptjatim)
coordinates(ptjatim)


colfunc=colorRampPalette(c("green", "yellow","red"))

color=colfunc(16)

ptjatim$Kemiskinan=artikel$Kemiskinan
ptjatim$APM=artikel$APM
ptjatim$APS=artikel$APS
ptjatim$APK=artikel$APK

names(ptjatim)
summary(ptjatim)

spplot(ptjatim, "Kemiskinan", col.regions=color,
       main="Peta Sebaran Kemiskinan Jatim 2020")

spplot(ptjatim, "APM", col.regions=color, colorkey=TRUE,
       main="Peta Sebaran APM Jatim 2020")

spplot(ptjatim, "APK", col.regions=color,
       main="Peta Sebaran APK Jatim 2020")

spplot(ptjatim, "APS", col.regions=color,
       main="Peta Sebaran APS Jatim 2020")

#REGRESI KLASIK
Y=matrix(c(Kemiskinan))
X=matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,APS),nrow=38,ncol=2)
Xt_X=t(X)%*%X
inv_XtX=solve(Xt_X)
Xt_Y=t(X)%*%Y
beta_topi=inv_XtX%*%Xt_Y
Ytopi=X%*%beta_topi
res=Y-Ytopi

var_res=(t(res)%*%res)/(38-1-1)
rerataY=mean(Y)
SST=sum((Y-rerataY)^2)
SSRes=sum((Y-Ytopi)^2)
SSreg=sum((Ytopi-rerataY)^2)

Fhit=(SSreg/1)/(SSreg/(38-1-1))

thit=(beta_topi[2,1])/(sqrt(var_res)*sqrt(inv_XtX[2,2]))

bobot=matrix(c(),nrow=38,ncol=38)

#bobot dari package spdep
#queen contiguity
wq=poly2nb(ptjatim)
#wq=poly2nb(ptjatim,queen = TRUE)

Wq=nb2listw(wq)
plot(ptjatim)
plot(Wq,coordinates(ptjatim),add=T)

#Dependensi
#Morans
pembilang=t(res)%*%Wq%*%res
penyebut=t(res)%*%res
Morans=pembilang/penyebut
#Uji signifikansi
Hat_matrix=X%*%inv_XtX%*%t(X)
M=diag(5)-Hat_matrix
#E(I)=tr(MW)/(n-k-1)
library(matrixcalc)
EXP_I=matrix.trace(M%*%Wq)/(5-1-1)
MWM_tW=M%*%Wq%*%M%*%t(Wq)
MWMW=M%*%Wq%*%M%*%Wq
MW=M%*%Wq
a1=matrix.trace(MWM_tW)
a2=matrix.trace(MWMW)
a3=(matrix.trace(MW))^2
Exp_II=(a1+a2+a3)/((5-1-1)*(5-1-1+2))
Var_I=Exp_II-(EXP_I)^2
Z_I=(Morans-EXP_I)/(sqrt(Var_I))

#Gagal Tolak H0 karena |Zhit|<1,96

Wq0=nb2listw(wq, style="B")

#rook contiguity
wr=poly2nb(ptjatim,queen = FALSE)
Wr=nb2listw(wr)

#


####Global Autokorelasi
moran(Kemiskinan, Wq, n=length(Wq$neighbours), S0=Szero(Wq))
moran.test(Kemiskinan, Wq, alternative="two.sided")

moran.test(Kemiskinan, Wq0, alternative="two.sided")
moran.plot(Kemiskinan,Wq,labels=`Lokasi`,main="Morans Plot")
?moran.plot

#####local autokorelasi
local_M=localmoran(Kemiskinan,Wq)
moran.map=cbind(ptjatim, local_M)
ptjatim$Ii=moran.map$Ii

library(raster)

quadrant=vector(mode="numeric",length=nrow(local_M))
m.qualification=Kemiskinan- mean(Kemiskinan) 
m.local_M=local_M[,1] - mean(local_M[,1]) 
signif=0.1
quadrant[m.qualification >0 & m.local_M>0] <- 4  
quadrant[m.qualification <0 & m.local_M<0] <- 1      
quadrant[m.qualification <0 & m.local_M>0] <- 2
quadrant[m.qualification >0 & m.local_M<0] <- 3
quadrant[local_M[,5]>signif] <- 0   
brks=c(0,1,2,3,4)
colors=c("white","blue",rgb(0,0,1,alpha=0.4),rgb(1,0,0,alpha=0.4),"red")
plot(ptjatim,border="lightgray",col=colors[findInterval(quadrant,brks,all.inside=FALSE)],main="Local Indicator LISA")
box()
legend("bottomleft", legend = c("insignificant","low-low","low-high","high-low","high-high"),
       fill=colors,bty="n")
####################################

regOLS=lm(Pengangguran~AMH+RataLamaSekolah+APK)
summary(regOLS)

lm.morantest(regOLS, Wq, alternative = "two.sided")
lm.LMtests(regOLS, Wq,test = c("all"))
lm.LMtests(regOLS, Wq,test = c("LMerr", "LMlag","SARMA"))

library(spatialreg)

#rho, theta, lambda

###SLX lag X, y=xb+WX(teta)+e
SLX=lmSLX(regOLS,data=ptjatim,Wq)
summary(SLX)
impacts(SLX, Wq)
summary(impacts(SLX, Wq,R=500),zstats=TRUE)

#SLX manual
X=model.matrix(regOLS)
head(X)
lagX1=create_WX(X,Wq,prefix = "lagx")
head(lagX1)
SLXdata=cbind(ptjatim,lagX1)
head(SLXdata)
names(SLXdata)
regSLX=lm(IPM~PDRB+AKB+lagx.PDRB+lagx.AKB,data=SLXdata)
summary(regSLX)

###SAR, y=(rho)Wy+xB+e
lag=lagsarlm(regOLS, listw=Wq, data=IPMjatim)
summary(lag)
impacts(lag,listw=Wq)
summary(impacts(lag, listw=Wq,R=500),zstats=TRUE)


###SEM, y=xB+u, u= (lambda)Wu+e
err=errorsarlm(regOLS, listw=Wq, data=IPMjatim) 
summary(err)
impacts(err,listw = Wq)

### SDEM, y=xB+WX(teta)+u, u=LWu+e
SDEM=errorsarlm(regOLS, data=IPMjatim, listw=Wq, etype="emixed")
summary(SDEM)
impacts(SDEM, listw=Wq)
summary(impacts(SDEM, listw=Wq,R=500),zstats=TRUE)

###SDM, y=(rho)Wy+xB+WX(teta)+e
SDM=lagsarlm(regOLS, data=IPMjatim, listw=Wq, type="mixed")
summary(SDM)
impacts(SDM, listw=Wq)
summary(impacts(SDM, listw=Wq,R=500),zstats=TRUE)

###Manski, all model, y=(rho)Wy+XB+WX(teta)+u,   u=(lambda)Wu+e
manski=sacsarlm(regOLS, data=IPMjatim, listw=Wq, type="sacmixed") 
summary(manski)
impacts(manski, listw=Wq)
summary(impacts(manski, listw=Wq,R=500),zstats=TRUE)

###SARMA, all model, y=(rho)Wy+XB+u,   u=(lambda)Wu+e
sarma=sacsarlm(regOLS, listw=Wq, data=IPMjatim)
summary(sarma)
impacts(sarma,listw=Wq)


##############


