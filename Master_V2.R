#Working code made by Jim as of 12/4
#Robert modifications 12/4 


#Read the dataset
baseball.dat = read.table(file.choose(),header=TRUE)

# -------------------------------------------------------------------
# PRE-DEFINE FUNCTIONS BEFORE MAIN CODE
# -------------------------------------------------------------------
# function to get separate response variable from predictors
ExtractResponseVariable <- function(dataset,name) {
  
  #Takes a dataframe, dataset and a name of a response variable
  #Extracts the response variable and dataframe of predictors, outputs these as members 
  #of a list
  
  if ( name %in% colnames(dataset)) {
    
    name <- as.character(name)
    
    #Get matrix of predictors
    predictors <- dataset
    predictors[name] <- NULL
    
    #Get response variable
    response <- dataset[name]
    
    return(list(response,predictors))
  }
  
  else {
    
    print(paste("Name ",name," not found in dataset",sep=''))
    return(list(0L,0L))
  }
}

# -------------------------------------------------------------------
#Evaluate the fitness of some model, output from lm or glm
#The userfunc should take a fitted model and output a scalar
#fitness value

FitnessFunction <- function(model,userfunc=FALSE){
  
  if (userfunc == FALSE) {
    fitness.value <- extractAIC(model)[2]
  }
  
  else {
    
    fitness.value <- userfunc(model)
  }
  
  return(fitness.value)
}



# -------------------------------------------------------------------
# function that determines 'fitness' of an invidivudal based on the quality 
# of the LS fit. The default for determining fitness is the Aikake Criteria Index
# but the user can supply their own custom-made fitness function
# **may be worth it to treat 'predictors' as global variable or object

AssessFitness <- function(individual, response, predictors, userfunc=FALSE){
  
  #Evaluate the fitness of some model, output from lm or glm
  #The userfunc should take a fitted model and output a scalar
  #fitness value
  
  #RMS simplified the following line
  predictors.individual <- predictors[,individual==1]
  
  model.out <- lm(response[,1]~., predictors.individual)
  fitness.value <- FitnessFunction(model.out,userfunc=userfunc)

  return(fitness.value)
}


# -------------------------------------------------------------------
# Example of user-supplied fitness function only for internal testing
# A test - this does exactly the same as the AIC function,
# but its user-define so can be used to test the fitness_function   #useage

TestUserFunc <- function (fit, scale = 0, k = 2) {
  
  n <- length(fit$residuals)
  edf <- n - fit$df.residual
  RSS <- deviance.lm(fit)
  dev <- if (scale > 0) 
    RSS/scale - n
  else n * log(RSS/n)
  return(dev + k * edf)
}


# -------------------------------------------------------------------
# Function that breeds P new children based on parents' genome and fitness

Breed <- function(generation, fitness.vec, predictors, prob.mute) {
  
  # generation is a list with each element containing the genome of an individual
  # fitness.vec is a vector
  prob.reproduction <- 2*rank(-fitness.vec)/(P*(P+1))
  parent.index.list <- lapply(1:P, function(x) sample(P,2,prob = prob.reproduction,replace=FALSE))
                              
  children <- lapply(parent.index.list, function(x) CrossOverMutate(generation, x, prob.mute))
  
  # return P children to be considered for selection
  # also return fitness evaluation
  return(children)
}


# -------------------------------------------------------------------
# Function that produces a single child from two chosen parents
# and allows for the possibility of mutation

CrossOverMutate <- function(generation, parent.index, prob.mutate){
  
  #Create child individual with half of its genetic material from parent1 and the other half from parent2
  #The generic material is chosen at random using sample
  parent1 <- generation[[parent.index[1]]]
  parent2 <- generation[[parent.index[2]]]
  
  child <- parent1 
  #generate locations of genetic information to swap
  pos <- sample(1:length(parent2),as.integer(length(parent2)/2),replace=FALSE)
  child[pos] <- parent2[pos]
  
  #generate mutation vector
  mutate = rbinom(length(child),1,prob.mutate)
  #do the mutation - this will ensure that if a 2 is produced, 
  #set to zero. If not, keeps as 1.
  child = (child+mutate)%%2
  
  return(child)
}

# -------------------------------------------------------------------
# MAIN PROGRAM
# -------------------------------------------------------------------

## Put all this in a function that can be called by user on the dataset
                     
# Define response and predictor variables
subsets <- ExtractResponseVariable(baseball.dat,"salary")

# Choose to scale or reject bad data based on boolean flag
flag.log.scale <- 1 
if (flag.log.scale) {
  response <- log(subsets[[1]])
} else {
  response <- subsets[[1]]
}
predictors <- subsets[[2]]

# Define/create key variables a priori
C <- length(predictors) #Get the number of predictors (GLOBAL)
P <- as.integer(C*1.5) #number of individuals in a given generation (GLOBAL)
Niter <- 60 #number of generation iterations to carry out (GLOBAL)
prob.mutate <- 1.0/(P*sqrt(C)) #mutation rate (should be about 1%) Formula suggested by G&H

                     
fitness <- matrix(0,P,Niter) #evolution of the fitness values over model run
frac.replace <- 0.2 # % of individuals in child/adult population selected/replaced

# Define first generation (without FOR loops, lists are preferred)
generation.old <- lapply(1:P, function(x) {rbinom(C,1,0.5)}) # list of individual genomes
                     
#assess fitness of the first generation                     
fitness[,1] <- sapply(generation.old, AssessFitness, response = response, predictors = predictors, userfunc = FALSE)


# -------------------------------------------------------------------
# MAIN LOOP for genetic algorithm
# put this in a loop function 
# Loop through generations and apply selective forces to create iterative generations
start <- Sys.time()
                     
for (n in 1:(Niter-1)) { #loop through fixed number of iterations
  
  # breed selection of P children and assess their fitness
  children <- Breed(generation.old, fitness[,n], predictors, mutation.rate)
  #generation.new <- children
  
  ## simplify so that we replace parents with children without combining the generations (for now)
  
  children.fitness <- sapply(children, AssessFitness, response = response, predictors = predictors, userfunc = FALSE)
  
  #children.best.index <- which(rank(-children.fitness)>round((1-frac.replace)*P)) # select best children to keep
  #children.best <- children[children.best.index] # vector length = # of adults to be replaced
  #children.fitness.best <- children.fitness[children.best.index] # vector length = # of adults to be replaced

  # now create new generation
  #generation.old.worst.index <- which(rank(-fitness[,n])<=round(frac.replace*P)) # select worst parent by rank
  
  generation.old <- children # keep most of prior generation
  fitness[,n+1] <- children.fitness # keep most of prior generation fitness data
  print(min(children.fitness))

}
stop <- Sys.time()
print(stop-start)

# -------------------------------------------------------------------
# Plot up envelope of fitness values
          
plot(-fitness,xlim=c(0,Niter),ylim=c(50,425),type="n",ylab="Negative AIC",
     xlab="Generation",main="AIC Values For Genetic Algorithm")
for(i in 1:Niter){points(rep(i,P),-fitness[,i],pch=20)}
           
# plot the fitness matrix to see how the entire population evolves over time             
                 
                    
# get fit for 'best' individual at end
generation.new <- generation.old

#faster way of getting the best index than using which 
best.index <- order(fitness[,Niter])[1]

best.individual <- generation.new[[best.index]]
print(best.individual)
predictors.individual <- predictors[,best.individual==1]
model.out <- lm(response[,1]~., predictors.individual)

summary(model.out)
