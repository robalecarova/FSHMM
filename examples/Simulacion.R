
# Simulacion de datos
rm(list = ls())

load("simulacionBase.RData")

library(MASS)

sigma = c(1, 15)
mean = c(0,0)
ordenEstados = c("N", "C")

estados = c("C","N","N","N")

n = length(estados)
replica = rep(0, n)
for (i in 1:length(estados)){
  estadoNumero = which.max(estados[i] == ordenEstados)
  replica[i] = rnorm(n = 1, mean[estadoNumero], sigma[estadoNumero]) / 10
}
  
replDelta = replica

offset = rnorm(1,10,1)
replica <- cumsum(c(offset, replica))

reps = 3
n = n+1
replicas = matrix(0,nrow = reps , ncol = n   )

for(i in 1:n){
  replicas[ ,i]  = rnorm(reps, replica[i], 0.5)
}

print(replDelta)
print(replica)
print(replicas)


