# Section 4.7.2 horseshoe crab example

crab <- read.table('C:/Users/vyabo/Desktop/AMS 573 Project/Crabs.dat',header=T)

glm1 <- glm(sat~width,family = poisson(link='log'),data=crab)
summary(glm1)$coefficients
# The proportion of residual deviance to df should be 1 if there is no overdispersion
summary(glm1)$deviance/summary(glm1)$df.residual

# Alternatively, if the model is of a poisson family, then we can use
 # an R package to check for overdispersion automatically
library(performance)
check_overdispersion(glm1)

sa <- tapply(crab$sat,crab$width,sum)
mu <- tapply(predict(glm1,type='response'),crab$width,sum)
chi.square <- sum((sa-mu)^2/mu)
phi <- chi.square/(66-2)
SE.quasi <- sqrt(phi)*summary(glm1)$coefficients[2,2]
chi.square
phi
SE.quasi

glm1.quasi <- glm(sat~width,family=quasipoisson(link = 'log'),data=crab)

confint(glm1)
confint(glm1.quasi)
CI.Length <- confint(glm1)[2,2]-confint(glm1)[2,1]
CI.Length
CI.Length.quasi <- confint(glm1.quasi)[2,2]-confint(glm1.quasi)[2,1]
CI.Length.quasi

################################################################################

# Section 4.7.4 teratology example

library(VGAM)
ter <- lirat
ter$grp <- as.factor(ter$grp)

glm2 <- glm(R/N~grp-1,weights=N,data=ter,family=binomial)
summary(glm2)$coefficients
# residual deviance to df should be 1 if no overdispersion
summary(glm2)$deviance/summary(glm2)$df.residual

pred <- unique(predict(glm2,type='response'))
SE <- sqrt(pred*(1-pred)/tapply(ter$N,ter$grp,sum))
X2 <- sum(resid(glm2,type='pearson')^2)
phi <- X2/(58-4)
pred
SE
X2
phi

SE*sqrt(phi)

glm2.quasi <- glm(R/N~grp-1,weights=N,data=ter,family=quasibinomial(link='identity'),start=pred)
summary(glm2.quasi)$coefficients

p1 <- sum(ter[ter$grp==1,]$R)/sum(ter[ter$grp==1,]$N)
p2 <- sum(ter[ter$grp==2,]$R)/sum(ter[ter$grp==2,]$N)
n1plus <- sum(ter[ter$grp==1,]$N)
n2plus <- sum(ter[ter$grp==2,]$N)
Wald_lower <- (p1-p2)-1.96*sqrt((p1*(1-p1))/n1plus+(p2*(1-p2))/n2plus)
Wald_upper <- (p1-p2)+1.96*sqrt((p1*(1-p1))/n1plus+(p2*(1-p2))/n2plus)

p1_quasi <- summary(glm2.quasi)$coefficients[1,1]
p2_quasi <- summary(glm2.quasi)$coefficients[2,1]
quasi_lower <- (p1_quasi-p2_quasi)-1.96*sqrt(phi)*sqrt((p1_quasi*(1-p1_quasi))/n1plus+(p2_quasi*(1-p2_quasi))/n2plus)
quasi_upper <- (p1_quasi-p2_quasi)+1.96*sqrt(phi)*sqrt((p1_quasi*(1-p1_quasi))/n1plus+(p2_quasi*(1-p2_quasi))/n2plus)

c(Wald_lower=Wald_lower,Wald_upper=Wald_upper)
c(Quasi_lower=quasi_lower,Quasi_upper=quasi_upper)

CI.Length <- Wald_upper-Wald_lower
CI.Length.quasi <- quasi_upper-quasi_lower
CI.Length
CI.Length.quasi

