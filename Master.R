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


m <- length(predictors) #Get the number of predictors
P <- as.integer(m*1.5) #number of samples of different predictors (generation size)
iter <- 100 #number of generation iterations to carry out
mutate_rate <- 1.0/(P*sqrt(m)) #mutation rate (should be about 1%)
generation_1 <- matrix(0,P,m) # matrix of the factor choices (P rows, m cols)
generation_2 <- matrix(0,P,m) #matrix of the improved factor choices
fitness_values <- matrix(0,P,1) #matrix of the fitness values for one generation
fitness_evolution <- matrix(0,P,iter) #evolution of the fitness values over model run
best_fitness <- 0 #best fitness value 
best_fitness_generation = rep(0,iter) #evolution of best fitness values


#Starting generation : generate P random combinations of the factors, do the regression and record their fitness values

for (i in 1:P){
  
  #Geneate a vector of 1s and 0s corresponding to which predictors we want to use
  generation_1[i,] <- rbinom(m,1,0.5)
  
  #Select the factors we want to use in the regression
  run_vars <- predictors[,generation_1[i,]==1]
  
  #linear regression of the response variable against the predictor   #data
  
  # Note what this means
  # lm( Y ~ A + B + C + D etc, data = run_vars)
  # where A,B,C,D etc are the names of columns within the dataframe
  # run_vars. The syntax here is just a more concise way or writing
  # it
  
  g = lm(response[,1]~.,run_vars)
  #g = glm(response[,1]~.,family=gaussian,data=run_vars)
  fitness_value <- fitness_function(g,userfunc = testuserfunc)
  fitness_values[i,] <- fitness_value
  fitness_evolution[i,1] <- fitness_value
  
  if (fitness_value < best_fitness){
    best_fitness <- fitness_value
    best_set <- generation_1[i,]
  }
  
} 

print(fitness_value)

#Some other notes:

#Cross over could be replaced with random sampling without repacement from one individual and then combine with 
#half the genes from the other parent. 
#This creates P children, use tornament matching so that we keep some of the parents (helps to reduce the likelihood of children being worse than the parents)
#Allow mutation - user defined or 1/m
#Criteria for convergence (min, max and varience of the AIC). Keep track of the variance of generation i and i-1 and then stop when this difference is small enough

