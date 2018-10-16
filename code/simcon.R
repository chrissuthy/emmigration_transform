source("mp.sim.R")

tmp <- mp.sim(nyears = 100)
plot(apply(tmp,2,sum)/nrow(tmp),ylim=c(0,1),type="l", col="grey90",
     ylab="Occupancy", xlab="Years", las=1)

for(i in 1:10){
  tmp <- mp.sim(nyears = 100, gamma=0)
  lines(1:ncol(tmp),apply(tmp,2,sum)/nrow(tmp),ylim=c(0,1),type="l", col=3)
}


for(i in 1:10){
  tmp <- mp.sim(nyears = 100, gamma = 1)
  lines(1:ncol(tmp),apply(tmp,2,sum)/nrow(tmp),ylim=c(0,1),type="l", col=2)
}

for(i in 1:10){
  tmp <- mp.sim(nyears = 100, gamma = 5)
  lines(1:ncol(tmp),apply(tmp,2,sum)/nrow(tmp),ylim=c(0,1),type="l", col=4)
}

for(i in 1:10){
  tmp <- mp.sim(nyears = 100, gamma = 0.00001)
  lines(1:ncol(tmp),apply(tmp,2,sum)/nrow(tmp),ylim=c(0,1),type="l", col="blue")
}
