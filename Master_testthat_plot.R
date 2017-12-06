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
  
  #Replace mate with cross_over_mutate from RMS code 
  children <- lapply(parent.index.list, function(x) CrossOverMutate(generation, x, prob.mute))
  
  # return P children to be considered for selection
  # also return fitness evaluation
  return(children)
}


# -------------------------------------------------------------------
# Function that produces a single child from two chosen parents
# and allows for the possibility of mutation

CrossOverMutate <- function(generation, parent.index, prob.mutate, flag.seed=FALSE){
  if(flag.seed){
    set.seed(1)
  }
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

####test that
#test the assess fitness function
library(testthat)
test_that("Check type of return value of AssessFitness function",{
  expect_equal(class(AssessFitness(individual=rbinom(1:C, 1, 0.5), response=response, predictors, FALSE)), "numeric" )
  
  
})


#test Breed function
test_that("Check the dimension of return list from Breed function is P embeded lists, 
          with length of each element(child genome) is C",{
            expect_equal(length(Breed(generation.new, fitness[,1], predictors, prob.mute)), P)
            expect_equal(unique(unlist(
              lapply(Breed(generation.new, fitness[,1], predictors, prob.mute),function(x) length(x))
            )), C)
            
            
          })


#test CrossOverMutate function
test_that("Check the dimension of return value(a child genome) from Cross_over_mutate function is C",{
  expect_equal(length(CrossOverMutate(generation.new, c(1,2), prob.mutate)), C)
})


###main function
variable.selection <- function(data, response.pos, flag.log.scale = 1, Ngen = 100){
  #let user define which variable to be the response variable
  resoponse.name <- names(data)[response.pos]
  response.name <- as.character(response.name) 
  subsets <- ExtractResponseVariable(data,"response.name")
  # Choose to scale or reject bad data based on boolean flag
  flag.log.scale <- 1 
  if (flag.log.scale) {
    response <- log(subsets[[1]])
  } else {
    response <- subsets[[1]]
  }
  predictors <- subsets[[2]]
  
  # Define/create key variables a priori
  C <<- length(predictors) #Get the number of predictors (GLOBAL)
  P <<- as.integer(C*1.5) #number of individuals in a given generation (GLOBAL)
  Ngen <<- 100 #number of generation iterations to carry out (GLOBAL)
  mutation.rate <<- 1.0/(P*sqrt(C)) #mutation rate (should be about 1%) Formula suggested by G&H
  
  fitness <<- matrix(0,P,Ngen) #evolution of the fitness values over model run
  frac.replace <<- 0.2 # % of individuals in child/adult population selected/replaced
  
  # Define first generation (without FOR loops, lists are preferred)
  generation.old <- lapply(1:P, function(x) {rbinom(C,1,0.5)}) # list of individual genomes
  
  #assess fitness of the first generation                     
  fitness[,1] <- sapply(generation.old, AssessFitness, response = response, predictors = predictors, userfunc = FALSE)
  
  for (n in 1:(Ngen-1)) { #loop through fixed number of iterations
    
    # breed selection of P children and assess their fitness
    children <- Breed(generation.old, fitness[,n], predictors, mutation.rate)
    generation.new <- children
    
    ## simplify so that we replace parents with children without combining the generations (for now)
    
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

}

#plot
library(ggplot2)

max.generation <- apply(m, 2, max)
min.generation <- apply(m, 2, min)
dat <- data_frame(Max = max.generation, Min = min.generation, Generation =seq(1, 100, by =1))

ggplot(dat, aes(Generation)) + 
  geom_point(aes(y = Max, colour = "Max")) + 
  geom_point(aes(y = Min, colour = "Min"))









