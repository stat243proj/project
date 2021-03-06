#CHANGES: Jim Falter by 5pm on F Dec 8
# 1) changed default input to "AIC" for FitnessFunction()
# 2) minor editing to clean up comments
# 3) removed code at end 'testing' algorithm.  Should now be done outside of code.
# 4) Changed Breed() so it always generates P children so main code doesn't have to adjust
# 5) Changed FitnessFunction() so it either uses default AIC or user submitted
#    function which it firsts tests is a function object and returns error otherwise
# 6) Changed AssessFitness() and GeneticAlgorithmFit() so defulat userfunc="AIC" to
#    be submitted to downstream functions.
# 7) Replaced 'plot' variable with 'plot.flag' since plot() is a common R function
#    This is used to decide whether we want to show intermediate output
# 8) Changed name of main function to Select()


# -------------------------------------------------------------------
# PRE-DEFINE FUNCTIONS BEFORE MAIN CODE

# -------------------------------------------------------------------
# ExtractResponseVariable()
# function to get separate response variable from predictors
ExtractResponseVariable <- function(dataset,name) {

  #Takes a dataframe, dataset and a name of a response variable
  #Extracts the response variable and dataframe of predictors, outputs these as members
  #of a list
  if (any(is.na(dataset[,1]))==TRUE){
    print("The response variable has missing values")
    dataset <- dataset[!is.na(dataset[,1])]
  }
  else{
    if (name %in% colnames(dataset)) {

      name <- as.character(name)

      #Get matrix of predictors
      predictors <- dataset
      predictors[name] <- NULL

      #Get response variable
      response <- dataset[name]
      return(list(response,predictors))
    } else {
      print(paste("Name ",name," not found in dataset",sep=''))
      return(list(0L,0L))
    }
  }
}


# -------------------------------------------------------------------
# FitnessFunction()
# Evaluate the fitness of some model, output from lm or glm
# based on three choices of fitness function
# 1) Aikake Information Crieter  or "AIC" (default)
# 2) Bayesian Information Criteria or "BIC"
# 3) Residual sum of squares or "Residual"

FitnessFunction <- function(model, userfunc){
  
  if (userfunc == "AIC") { # default case of Aikake Info Criteria
    
    fitness.value <- extractAIC(model)[2]
    return(fitness.value)
    
  } else {
    if (class(userfunc) == "function") { # test if submitted function is a function object
      
      userfunc(model)
      return(fitness.value) 
      
    } else {
      print(paste0("WARNING: ", deparse(substitute(unserfunc)), "() is not a recognized function"))
    }
  }
}


# -------------------------------------------------------------------
# AssessFitness()
# function that determines 'fitness' of an invidivudal based on the quality
# of the LS fit. The default for determining fitness is the Aikake Criteria Index
# but the user can supply their own custom-made fitness function
# **may be worth it to treat 'predictors' as global variable or object
AssessFitness <- function(individual, response, user.family="gaussian", predictors, userfunc="AIC"){

  #Evaluate the fitness of some model, output from glm
  #The userfunc should take a fitted model and output a scalar fitness value

  #ensure that there is more than one '1' in the vector
  while (sum(individual) <= 1){

    insertonepos = sample(1:length(individual),1)
    individual[insertonepos] = 1
  }
  #print(individual)

  predictors.individual <- predictors[,individual==1]

  #Check distribution family of glm()
  if(!user.family %in% c("binomial", "gaussian", "Gamma", "poisson", "inverse.gaussian")){
    print(paste("WARNING: User defined distribution family ", user.family, " does not exist", sep=''))
    stop()
    geterrmessage()
  }
  else{
    model.out <- glm(response[,1]~., family=user.family, data=predictors.individual)
    #model.out <- lm(response[,1]~., predictors.individual)
    fitness.value <- FitnessFunction(model.out, userfunc)
  }
  
  return(fitness.value)
}


# -------------------------------------------------------------------
# Breed()
# Function that breeds P new children based on parents' genome and fitness

Breed <- function(generation, fitness.vec, predictors, prob.mutate) {

  # generation is a list with each element containing the genome of an individual
  # fitness.vec is a vector
  prob.reproduction <- 2*rank(-fitness.vec)/(P*(P+1))
  parent.index.list <- lapply(1:P, function(x) sample(P,2,prob = prob.reproduction,replace=FALSE))

  child.pairs <- lapply(parent.index.list, function(x) CrossOverMutate(generation, x, prob.mutate))
  # collate list of P child pairs into a single generation (list) of P children
  children <- unlist(child.pairs, recursive = FALSE)[1:P/2]
  # children <- lapply(child.pairs, function(x) x[[1]])

  return(children)
}


# -------------------------------------------------------------------
# CrossOverMutate()
# Function that produces two children with mirroring genomes from two chosen parents
# and simultaneously allows for the possibility for random mutation
CrossOverMutate <- function(generation, parent.index, prob.mutate=0.005){

  # Create child individual with half of its genetic material from parent1 and the other half from parent2
  # The generic material is chosen at random using sample
  parent1 <- generation[[parent.index[1]]]
  parent2 <- generation[[parent.index[2]]]
  child1 <- parent1
  child2 <- parent2
  
  # generate locations of genetic information to swap
  pos <- sample(1:length(parent2),as.integer(length(parent2)/2),replace=FALSE)
  child1[pos] <- parent2[pos]
  child2[pos] <- parent1[pos]

  # generate mutation vector
  mutate1 = rbinom(length(child1),1,prob.mutate)
  mutate2 = rbinom(length(child2),1,prob.mutate)
  
  # do the mutation - this will ensure that if a 2 is produced,
  # set to zero. If not, keeps as 1.
  child1 = (child1+mutate1)%%2
  child2 = (child2+mutate2)%%2
  child.pair = list(child1, child2)

  return(child.pair)
}


# -------------------------------------------------------------------
# RemoveClones()
# function that removes any clones from a given generation and replaces them
# with individuals randomly created from the entire genome clones are predicted
# to exist when fitness of two individuals are exactly the same
# this is highly unlikely unless they share the same exact genome

ReplaceClones <- function(generation, fitness.vec, C) {

  # function that removes any clones from a given generation and replaces them
  # with individuals randomly created from the entire genome
  # clones are predicted to exist when fitness of two individuals are exactly the same
  # this is highly unlikely unless they share the same exact genome
  clone.index <- which(duplicated(fitness.vec))
  N.clones <- length(clone.index)
  replacements <- lapply(1:N.clones, function(x) {rbinom(C,1,0.5)}) # list of new individual genomes
  generation[clone.index] <- replacements
  # the followin is to avoid computing fitness for the majority of non-clones
  fitness.replacements <- sapply(replacements, AssessFitness, response = response, predictors = predictors)
  dim(fitness.replacements)
  fitness.vec[clone.index] <- fitness.replacements
  output <- list("generation" = generation, "fitness" = fitness.vec)
  return(output)
}


# -------------------------------------------------------------------
# ExtractBestIndividual()
# Function that extracts the best individual and its corresponding fitness, and prints a set of
# summary statistics
ExtractBestIndividual <- function(generation, fitnessmatrix, plot.flag=FALSE){

  #Extract the best individual and its corresponding fitness, and print a set of
  #summary statistics

  best.index <- order(fitnessmatrix[,Niter])[1]

  best.individual <- generation[[best.index]]
  predictors.individual <- predictors[,best.individual==1]
  best.model <- lm(response[,1]~., predictors.individual)
  if (plot.flag == TRUE) {
    print(summary(best.model))
    print(best.individual)
  }
  return(best.model)
}


# -------------------------------------------------------------------
#MAIN FUNCTION() TO APPLY GENETIC ALGORITHM
Select <- function(dataset, response.name, userfunc="AIC", user.family="gaussian", flag.log.scale=TRUE, 
                                frac.replace=0.2, Niter=100, mutate.rate=FALSE, plot.flag=TRUE){

  #User can define a fitness function, log scale flag, fraction of children to replace with parents, number of iterations and a mutation rate. If these are not provided they are set to dafault

  # Define response and predictor variables
  subsets <- ExtractResponseVariable(dataset, response.name)

  # Choose to scale or reject bad data based on boolean flag

  if (flag.log.scale == TRUE) {
    response <<- log(subsets[[1]])
  } else {
    response <<- subsets[[1]]
  }
  predictors <<- subsets[[2]]

  # Define/create key variables a priori
  # These variables are accessible to all genalg functions

  C <- length(predictors) #Get the number of predictors (GLOBAL)
  Niter <<- Niter #number of iterations
  P <- as.integer(C*1.5) #number of individuals in a given generation (GLOBAL)
  P <<- 2*ceiling(P/2) # Force P to be even due to a future need to 'split' a generation in two

  #Set the mutation rate
  if (mutate.rate == FALSE) {
    prob.mutate <<- 1.0/(P*sqrt(C)) #mutation rate (should be about 1%) Formula suggested by G&H
  }
  else {
    prob.mutate <<- mutate.rate
  }

  fitness <<- matrix(0,P,Niter) #evolution of the fitness values over model run

  # Define first generation
  generation.old <<- lapply(1:P, function(x) {rbinom(C, 1, 0.5)}) # list of individual genomes

  #assess fitness of the first generation
  fitness[,1] <<- sapply(generation.old, AssessFitness, response = response, user.family, predictors = predictors, userfunc)

  # -------------------------------------------------------------------
  # MAIN LOOP for genetic algorithm
  # Loop through generations and apply selective forces to create iterative generations

  start <- Sys.time()

  for (n in 1:(Niter-1)) { #Niter -1 because we've already made the first generation

    # breed selection of P children and assess their fitness
    children <- Breed(generation.old, fitness[,n], predictors, prob.mutate)

    children.fitness <- sapply(children, AssessFitness, response = response, user.family, predictors = predictors, userfunc)
    number.children.keep <- round((1-frac.replace)*P)
    number.parents.keep <- P - number.children.keep

    # If we do want to keep parents in the new generation, then figure out the parents that we want to
    # keep. Otherwise, just replace all the parents with the children.
    if (number.parents.keep > 1){

      parents.fitness <- sapply(generation.old, AssessFitness, response = response, user.family, predictors = predictors, userfunc)

      children.best.index <- order(children.fitness)[1:number.children.keep] # select best children to keep
      children.best <- children[children.best.index] # select the children to keep
      children.fitness.best <- children.fitness[children.best.index] # select fitness of best children

      parents.best.index <- order(parents.fitness)[1:number.parents.keep] # get indices of best parents
      parents.best <- generation.old[parents.best.index] # select the parents to keep
      parents.fitness.best <- parents.fitness[parents.best.index] # select the fitness of the best parents

      #Create nre generation and new generation fitness by concatinating the vectors we made
      generation.new.fitness <- c(children.fitness.best,parents.fitness.best)
      generation.new <- c(children.best,parents.best)
    }

    else{
      generation.new <- children
      generation.new.fitness <- children.fitness
    }

    # check next generation for clones and replace them with random individuals if necessary
    clones.removed <- ReplaceClones(generation.new, generation.new.fitness, C)
    generation.new <- clones.removed$generation
    generation.new.fitness <- clones.removed$fitness

    #generation.old.worst.index <- which(rank(-fitness[,n])<=round(frac.replace*P)) # select worst parent by rank
    generation.old <- generation.new # keep most of prior generation
    fitness[,n+1] <- generation.new.fitness # keep most of prior generation fitness data

  }
  stop <- Sys.time()

  best.model <- ExtractBestIndividual(generation.old,fitness, plot.flag)


  #If user wants to plot the evolution of the population over time, do so

  if (plot.flag == TRUE) {
    
    print(min(generation.new.fitness))
    plot(-fitness,xlim=c(0,Niter),ylim=c(min(fitness), max(fitness)),type="n",ylab="Negative fitness value",
         xlab="Generation",main="Fitness values For Genetic Algorithm")
    for(i in 1:Niter){points(rep(i,P),-fitness[,i],pch=20)}
  }


  # show run time
  cat("Algorithm runtime: ", round(as.numeric(stop-start),2), " seconds\n")
  # print(stop-start)

  # RETURN OUTPUT VARIABLES
  output <- list("LastGen" = generation.new, "fitness" = fitness, "bestModel" = best.model)
  return(output)

}













