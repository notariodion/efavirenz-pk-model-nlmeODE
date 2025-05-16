setwd("D:\\GitHub\\Farmakokinetik Efavirenz")

# call packages
library(nlme)
library(nlmeODE)
library(readxl)

# data preparation
efavirenz <- read_excel("efavirenz.xlsx")
data01 <- groupedData(conc ~ Time | Subject,
                      data = as.data.frame(efavirenz))
data01$Dose[data01$Time!=0] <- 0
data01$Cmt <- rep(1,dim(data01)[1])
data01$conc

# plotting data
plot(data01, ylab='plasma concentration (mcg/mL)', xlab = "time after dosing (h)")

# two compartment model
twoComp <- list(DiffEq=list(               
      dTAdt = ~ -ka*TA,     
      dD1dt = ~ ka*TA-k10*D1-k12*D1+k21*D2,
      dD2dt = ~ k12*D1-k21*D2),
      ObsEq=list(                
          TA = ~ 0,
          C1 = ~ D1/V1,
          C2 = ~ 0),
     Parms=c("ka","k10", "k12", "k21", "V1"),   
     States=c("TA","D1","D2"),       
     Init=list(0,0,0))

# fitting model to data
efa.model <- nlmeODE(twoComp, data01, LogParms = T)
efa.nlme01 <- nlme(conc ~ efa.model(ka, k10, k12, k21, V1, Time, Subject),
                 data = data01, fixed=ka+k10+k12+k21+V1~1, random = pdDiag(ka+k10+k12+k21+V1~1), 
                 start=c(ka = -0.8, k10 = -3.0, k12 = -0.5, k21 = -2.5, V1=3.5),
                 control=list(returnObject=TRUE,msVerbose=TRUE),
                 verbose=TRUE)

# reduce random effect
efa.nlme02 <- update(efa.nlme01, random = pdDiag(ka+k12~1))
efa.nlme03 <- update(efa.nlme01, random = pdDiag(ka~1))

# add covariate
efa.nlme04 <- update(efa.nlme01, fixed = list(ka+k12+k21+k10~1, V1~Wt), 
                  start = c(ka = -0.5, k10 =0.2, k12 = 1.5, 
                            k21 = -2.5, v1=5, Wt = -0.1))
efa.nlme05 <- update(efa.nlme01, fixed = list(ka+k12+k21+V1~1, k10~Age), 
                     start = c(ka = 1.0, k10 =-3.0, k12 = -1.5, 
                               k21 = -2.5, v1=10, Age = 0.1))

# check the best model
anova(efa.nlme01, efa.nlme02, efa.nlme03, efa.nlme04)
anova(efa.nlme01, efa.nlme05, efa.nlme03, efa.nlme04)

# diagnostic plot
plot(augPred(efa.nlme01,level=0:1),
      xlab="time after dosing (h)",
      ylab="plasma concentration (mcg/mL)")
plot(efa.nlme01, conc~fitted(.,0), abline =
         c(0,1), xlab='predicted concentration',
       ylab='measured concentration',  main='Population')
plot(efa.nlme01, conc~fitted(.,1), abline =
         c(0,1), xlab='predicted concentration',
       ylab='measured concentration', main = "Individual")
plot(efa.nlme01, resid(.,0,type ='n')~fitted(.,1),id = 0.05, abline = 0,
       xlab="predicted concentration",
       ylab="standardized residual", main='Population')
plot(efa.nlme01, resid(.,1,type ='n')~fitted(.,1),id = 0.05, abline = 0,
       xlab="predicted concentration",
       ylab="standardized residual", main ="Individual")

# summary result
summary(efa.nlme01)
write.csv(exp(coef(efa.nlme01)) , file= "coef_efa.csv")
