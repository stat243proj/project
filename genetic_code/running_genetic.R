#Read the dataset
baseball = read.table(file.choose(),header=TRUE)
dataset <- baseball

Nruns <- 10
Ngen <- 100
Master <- as.list(rep(NA, 10))
best.fit <- rep(NA, Nruns)
for (i in seq_len(Nruns)) {
  Master[[i]] <- Select(dataset=dataset, 
                             response.name="salary", 
                             user.family="gaussian",
                             flag.log.scale=TRUE,
                             Niter = Ngen,
                             frac.replace = 0.2,
                             mutate.rate = 0.005,
                             plot.flag=FALSE)
  best.fit[i] <- min(Master[[i]]$fitness[,Ngen])
  cat("Completed generation: ", i, "\n")
}