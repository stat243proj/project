#Working code made by Jim as of 12/4

# read in baseball data to test
library(dplyr)
baseball.dat = read.table("baseball.dat", header=TRUE)

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
# function that determines 'fitness' of an invidivudal based on the quality 
# of the LS fit. The default for determining fitness is the Aikake Criteria Index
# but the user can supply their own custom-made fitness function
# **may be worth it to treat 'predictors' as global variable or object

AssessFitness <- function(individual, response, predictors, userfunc=FALSE){
  
  #Evaluate the fitness of some model, output from lm or glm
  #The userfunc should take a fitted model and output a scalar
  #fitness value
  predictors.individual <- subset(predictors, select = which(individual==1))
  model.out <- lm(response[,1]~., predictors.individual)
  
  if (userfunc == FALSE) {
    fitness.value <- extractAIC(model.out)[2]
  } else {
    fitness.value <- userfunc(model.out)
  }
  
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
  parent.index.list <- lapply(1:P, function(x) sample(P, 2, prob = prob.reproduction))
  children <- lapply(parent.index.list, function(x) Mate(generation, x, prob.mute))
  
  # return P children to be considered for selection
  # also return fitness evaluation
  return(children)
}


# -------------------------------------------------------------------
# Function that produces a single child from two chosen parents
# and allows for the possibility of mutation

Mate <- function(generation, parent.index, prob.mute) {
  
  # first mate two parents with equal likelihood of genetic contribution
  parent1 <- generation[[parent.index[1]]]
  parent2 <- generation[[parent.index[2]]]
  parent1.index <- sample(C, C/2, replace = FALSE)
  parent2.index <- c(1:C)[-parent1.index]
  child <- rep(NA, C)
  child[parent1.index] <- parent1[parent1.index]
  child[parent2.index] <- parent2[parent2.index]
  
  # now allow for mutation (switch gene sign)
  mutation.index <- which(rbinom(C, 1, prob.mute)==1)
  child[mutation.index] <- !child[mutation.index]
  
  # return child genome
  return(child)
}

# -------------------------------------------------------------------
# MAIN PROGRAM
# -------------------------------------------------------------------

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
Ngen <- 100 #number of generation iterations to carry out (GLOBAL)
mutation.rate <- 1.0/(P*sqrt(C)) #mutation rate (should be about 1%) Formula suggested by G&H
generation.old <- matrix(0,P,C) # matrix of the factor choices (P rows, m cols)
generation.new <- matrix(0,P,C) #matrix of the improved factor choices
fitness <- matrix(0,P,Ngen) #evolution of the fitness values over model run
frac.replace <- 0.2 # % of individuals in child/adult population selected/replaced

# Define first generation (without FOR loops, lists are preferred)
generation.old <- lapply(1:P, function(x) {rbinom(C,1,0.5)}) # list of individual genomes
fitness[,1] <- sapply(generation.old, AssessFitness, response = response, predictors = predictors, userfunc = FALSE)


# -------------------------------------------------------------------
# MAIN LOOP for genetic algorithm
# Loop through generations and apply selective forces to create iterative generations
start <- Sys.time()
for (n in seq_len(Ngen-1)) { #loop through fixed number of iterations
  
  # breed selection of P children and assess their fitness
  children <- Breed(generation.old, fitness[,n], predictors, mutation.rate)
  children.fitness <- sapply(children, AssessFitness, response = response, predictors = predictors, userfunc = FALSE)
  children.best.index <- which(rank(-children.fitness)>round((1-frac.replace)*P)) # select best children to keep
  children.best <- children[children.best.index] # vector length = # of adults to be replaced
  children.fitness.best <- children.fitness[children.best.index] # vector length = # of adults to be replaced

  # now create new generation
  generation.old.worst.index <- which(rank(-fitness[,n])<=round(frac.replace*P)) # select worst parent by rank
  generation.new <- generation.old # keep most of prior generation
  fitness[,n+1] <- fitness[,n] # keep most of prior generation fitness data
  generation.new[generation.old.worst.index] <- children.best
  fitness[generation.old.worst.index,n+1] <- children.fitness.best
}
stop <- Sys.time()
print(stop-start)

# -------------------------------------------------------------------
# Plot up envelope of fitness values
fit.min <- apply(fitness, min, MARGIN = 2)
fit.max <- apply(fitness, max, MARGIN = 2)
x.gen <- seq_len(Ngen)
plot(x.gen, fit.min, type = "l")
lines(x.gen, fit.max, lty = 2, type = "l")

# get fit for 'best' individual at end
best.index <- which(rank(fitness[,Ngen])==1)
best.individual <- generation.new[[best.index]]
predictors.individual <- subset(predictors, select = which(best.individual==1))
model.out <- lm(response[,1]~., predictors.individual)