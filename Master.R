#Working code 

extract_response_variable <- function(dataset,name) {
  
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

fitness_function <- function(model,userfunc=FALSE){
  
  #Evaluate the fitness of some model, output from lm or glm
  #The userfunc should take a fitted model and output a scalar
  #fitness value
  
  if (userfunc == FALSE) {
    fitness_value <- extractAIC(g)[2]
  }

  else {
    
    fitness_value <- userfunc(model)
  }

return(fitness_value)
}

testuserfunc <- function (fit, scale = 0, k = 2) 
{
  
  # A test - this does exactly the same as the AIC function,
  #but its user-define so can be used to test the fitness_function   #useage
  
  n <- length(fit$residuals)
  edf <- n - fit$df.residual
  RSS <- deviance.lm(fit)
  dev <- if (scale > 0) 
    RSS/scale - n
  else n * log(RSS/n)
  return(dev + k * edf)
  
}


#Make matices to store info about the the evolving models and their fitness values
subsets <- extract_response_variable(baseball_dat,"salary")

#Choose to scale or reject bad data.

#Assess fitness 
#Other fitness test
#Tornament matching to prevent children being 'worse' than the parents


response <- log(subsets[[1]])
predictors <- subsets[[2]]

# define key variables a priori
C <- length(predictors) #Get the number of predictors
P <- as.integer(m*1.5) #number of individuals in a given generation
Ngen <- 100 #number of generation iterations to carry out
mutate_rate <- 1.0/(P*sqrt(C)) #mutation rate (should be about 1%)
parents <- matrix(0,P,C) # matrix of the factor choices (P rows, m cols)
children <- matrix(0,P,C) #matrix of the improved factor choices
fitness <- matrix(0,P,Ngen) #evolution of the fitness values over model run
# best_fitness <- 0 #best fitness value 
best_fitness_generation = rep(0,Ngen) #evolution of best fitness values


#Starting generation : generate P random combinations of the factors, do the regression and record their fitness values
for (i in 1:P){
  
  #Geneate a vector of 1s and 0s corresponding to which predictors we want to use
  parents[i,] <- rbinom(m,1,0.5)
  
  #Select the factors we want to use in the regression
  predictors_chosen <- predictors[, parents[i,]==1]
  
  #linear regression of the response variable against the predictor   #data
  
  # Note what this means
  # lm( Y ~ A + B + C + D etc, data = run_vars)
  # where A,B,C,D etc are the names of columns within the dataframe
  # run_vars. The syntax here is just a more concise way or writing
  # it
  
  model.out = lm(response[,1]~., predictors_chosen)
  #g = glm(response[,1]~.,family=gaussian,data=run_vars)
  fitness[i,1] <- fitness_function(model.out, userfunc = testuserfunc)
} 

# examine initial fitness (debugging)
print(fitness[,1])

for (n in 1:Ngen) { #loop through fixed number of iterations
   # 1) produce P new children
      # a) rank parents and generate probability of reproduction
      # b) select 2 parents based on multinomial probability from above 
      # c) randomly select genetic material from first parent using binomial with p = 0.5
      # d) get alleles for all non-chosen genes from second parent
      # e) create child
      # f) repeat P times
  # 2) Create new generation
      # a) rank children and take top X percent
      # b) replace worst X percent of parents with best X percent of children
  # 3) Assess fitness of new generation
}

#Some other notes:

#Cross over could be replaced with random sampling without repacement from one individual and then combine with 
#half the genes from the other parent. 
#This creates P children, use tornament matching so that we keep some of the parents (helps to reduce the likelihood of children being worse than the parents)
#Allow mutation - user defined or 1/m
#Criteria for convergence (min, max and varience of the AIC). Keep track of the variance of generation i and i-1 and then stop when this difference is small enough

