x <- list.files("/backup/Course Work/Spatial Dynamics/Project/data/",full.names = TRUE)
datasets <- list()

for(i in 1:length(x))
{
  datasets[[i]] <- read.csv(x[i])
}

plot(datasets[[1]][,2],datasets[[1]][,3],main = "Flocks of 200 paticles each,Density = 2.00",xlab = "Noise",ylab = "OrderParameter",xlim = c(0,7),ylim = c(0,1))
points(datasets[[8]][,2],datasets[[8]][,3],pch = 7)
points(datasets[[11]][,2],datasets[[11]][,3],pch = 18)

legend("topright",legend = c(expression(paste(delta,"r=0.0",sep = "")),expression(paste(delta,"r=0.7",sep = "")),expression(paste(delta,"r=1.0",sep = ""))),pch = c(1,7,18))


Coefficient of assortment (ask jaideep ) - between group variation in comparison to overall variation.
Is there a binomial distribution.
