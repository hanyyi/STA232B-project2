library(lme4)
lamb =read.csv("lamb_lab1.csv")
#preprocessing data
lamb$sire = factor(lamb$sire)
lamb$line = factor(lamb$line)
lamb$damage = factor(lamb$damage)
head(lamb)

model1 = lmer(weight ~ line + damage -1 + (1|sire),lamb,REML = TRUE)
summary(model1)

X = as.matrix(getME(model1,"X"))
Z = as.matrix(getME(model1,"Z"))
y = as.matrix(lamb$weight)

R_sigma_e = sigma(model1)^2
G_sigma_s = unlist(VarCorr(model1))
beta_hat = fixef(model1)
V = Z %*% (G_sigma_s*diag(23)) %*% t(Z) + R_sigma_e * diag(nrow(y))

#calculate eblup
eblup_rand = (G_sigma_s*diag(23)) %*% t(Z) %*% solve(V) %*% (y - X %*% beta_hat)
r1 = order(eblup_rand,decreasing = T)
sire1 = data.frame(Sire = r1,EBLUP = eblup_rand[r1])
ranef(model1)
ggplot(sire1,aes(x=Sire,y=EBLUP))+geom_point()+geom_point()+
  geom_point(data=sire1[1,],color="red",size=5)+
  geom_point(data=sire1[23,],color="green",size=5)+
  labs(x="Sire",y="EBLUP of each Sire")

#calculate eblup fo li
L = rep(NA,5)
M = 0
for (i in 1:5) {
  M = as.numeric(lamb$line == i)
  M= matrix(M/sum(M),nrow = 1)
  s = as.matrix(ranef(model1)$sire$`(Intercept)`)
  L[i] = M %*% X %*% beta_hat + M %*% Z %*% s
}

print(L)

o = order(L)
mixed_li = data.frame(order=o,Line=c(1:5),"EBLUP of mixed effect"=L,"EBLUP of Li"= as.numeric(fixef(model1))[1:5])
mixed_li[order(mixed_li$order),]

#

m = nlevels(lamb$sire)
b_phi = (R_sigma_e * G_sigma_s)/(R_sigma_e + table(lamb$sire)*G_sigma_s)

nij = as.numeric(table(lamb$sire))
fixef(model1)

s_i = matrix(NA,ncol = m, nrow = m)
s_ii = matrix(NA,ncol = m, nrow = m)
b_phi_i = matrix(NA,ncol = m, nrow = m)

for (i in 1:m) {
  lamb_i = lamb[-which(lamb$sire==i),]
  lamb_reml_i = lmer(weight ~ line + damage -1 + (1|sire), lamb_i, REML = T)
  R_sigma_ei = sigma(lamb_reml_i)^2
  G_sigma_si = unlist(VarCorr(lamb_reml_i))
  V_hat_i = G_sigma_si * Z %*% t(Z) + R_sigma_ei * diag(nrow(lamb))
  beta_i = fixef(lamb_reml_i)
  s_i[,i] = as.matrix(G_sigma_si*t(Z) %*% solve(V_hat_i) %*% (y - X %*% beta_i))
  b_phi_i[,i] = (R_sigma_ei * G_sigma_si)/(R_sigma_ei + nij * G_sigma_si)
}


MSPE = rep(NA,m)
for (j in 1:m) {
  MSPE[j] = b_phi[j] - (m-1)/m*sum(b_phi_i[j,] - b_phi[j]) + (m-1)/m*sum((s_i[j,] - s[j])^2)
}
MSPE

##
margin_of_error = sqrt(MSPE)
sire_rank = data.frame(order= c(1:23),r1,sire = eblup_rand[r1],L=eblup_rand[r1]-margin_of_error,U=eblup_rand[r1]+margin_of_error)
library(ggplot2)
ggplot(sire_rank,aes(x=order,y=sire,label=r1))+geom_pointrange(aes(ymin=L,ymax=U))+geom_text(hjust=0,vjust=-0.5)








