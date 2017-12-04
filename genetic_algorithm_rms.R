

#------------------------------------------------------------
#Separate the provided dataframe
#------------------------------------------------------------

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

#------------------------------------------------------------
#Generic fitness function 
#------------------------------------------------------------

fitness_function <- function(model,userfunc=FALSE){
  
  #Evaluate the fitness of some model, output from lm or glm
  #The userfunc should take a fitted model and output a scalar
  #fitness value
  
  if (userfunc == FALSE) {
    fitness_value <- extractAIC(model)[2]
  }
  
  else {
    
    fitness_value <- userfunc(model)
  }
  
  return(fitness_value)
}

#------------------------------------------------------------
#Test of generic fitness function
#------------------------------------------------------------

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

#------------------------------------------------------------
#Make data structures as global variables (best to do this within a class)
#------------------------------------------------------------

Set_up_data_structures <- function(C,iterations=100,mutate_rate=FALSE) {
  
  #Takes in the number of predictors, preallocates arrays to be filled
  #These must be set as global encironment variables
  #iter can be specified by the user, default is 100
  
  P <<- as.integer(C*1.5) #number of samples of different predictors (generation size)
  
  if (mutate_rate == FALSE){
    mutate_rate <<- 1.0/(P*sqrt(C)) #mutation rate (should be about 1%)
  }
  else {
    mutate_rate <<- mutate_rate
  }
  parents <<- matrix(0,P,C) # matrix of the factor choices (P rows, m cols)
  children <<- matrix(0,P,C) #matrix of the improved factor choices
  parents_fitness <<- rep(0,P) #vector of parent fitness 
  children_fitness <<- rep(0,P) #vector of children fitness
  predictor_vals <<- rep(0,C) #vector od values of the predictors used in one iteration 
  best_individial <<- rep(0,P) #best individual in one generation 
  iter <<- iterations #number of iterations (set by the user)
  Best_fitness_generation <<- rep(0,iter) #evolution of best fitness values
  Best_fitness <<- 0 #value of best fitness, updated every iteration 

}

#------------------------------------------------------------
#Make first generation
#------------------------------------------------------------

Generate_first_generation <- function() {
  
  #Fills the parents matrix with random individuals and evaluates their fitness
  
  for (i in 1:P){
    
    #Geneate a vector of 1s and 0s corresponding to which predictors we want to use
    parents[i,] <<- rbinom(C,1,0.5)
    
    #Select the factors we want to use in the regression
    predictor_vals <- predictors_df[,parents[i,]==1]
    
    #linear regression of the response variable against the predictor   #data
    
    # Note what this means
    # lm( Y ~ A + B + C + D etc, data = run_vars)
    # where A,B,C,D etc are the names of columns within the dataframe
    # run_vars. The syntax here is just a more concise way or writing
    # it
    
    g <- lm(response_df[,1]~.,predictor_vals)
    
    fitness_value <- fitness_function(g)
    parents_fitness[i] <- fitness_value

    if (fitness_value < Best_fitness){
      Best_fitness <<- fitness_value
      best_individual <<- parents[i,]
    }
    
  }
  
  #rank the fitness values in order of goodness and record the indices
  #of that order
  ordered_fitness <- rank(-parents_fitness)
  
  #Probability of selecting from a generation - proportional to the
  #fitness of the individals
  phi <<- 2*ordered_fitness/(P*(P+1))
  Best_fitness_generation[1] <<- Best_fitness
}

#------------------------------------------------------------
#Basic function to do cross-over in the way we talked about
#(random selection)
#------------------------------------------------------------

Cross_over_mutate <- function(parent1,parent2){
  
  #Create chld individual with half of its genetic material from parent1 and the other half from parent2
  #The generic material is chosen at random using sample
  child <- parent1 
  pos <- sample(1:length(parent2),as.integer(length(parent2)/2),replace=FALSE)
  child[pos] <- parent2[pos]
  
  #generate mutation vector
  mutate = rbinom(length(child),1,mutate_rate)
  #do the mutation - this will ensure that if a 2 is produced, 
  #set to zero. If not, keeps as 1.
  child = (child+mutate)%%2
  
  return(child)
}


############################################################
#Not working 

Assess_generation <- function(generation_matrix) {
  
  #Determine the fitness vector for a set of individials and return a vector if their ordered indices
  
  for (i in 1:P){
    
    #Select the factors we want to use in the regression
    predictor_vals <- predictors_df[,generation_matrix[i,]==1]
    
    g <- lm(response_df[,1]~.,predictor_vals)
    
    fitness_value <- fitness_function(g)
    fitness_vector[i] <- fitness_value
  
  }
  
  ordered_fitness <- order(-fitness_vector)
  return(list(fitness_vector,ordered_fitness))
    
}


Concatinate_generation <- function(parents,parents_order,children,children_order) {
  
  preserved_prop <- 0.1
  number_preserved_parents <- as.integer(preserved_prop*P)
  number_preserved_children <- P - number_preserved_parents
  
  print(number_preserved_children)
  print(number_preserved_parents)
  
  children_new <- children
  children_new_fitness <- c(parents_order[[1]][1:number_preserved_parents],children_order[[1]][1:number_preserved_children])
  children_new[1:number_preserved_parents,] <- parents[parents_order[[2]][1:number_preserved_parents],]
  children_new[(number_preserved_parents+1):P,] <- children[children_order[[2]][1:number_preserved_children],]
  
  return(list(children_new,children_new_fitness))
  
}
############################################################

#------------------------------------------------------------
#Iterations loop
#------------------------------------------------------------

mainloop <- function(){
  
  for (i in 1:(iter-1)){
    
    #Select two parents using probability phi, to create a child generation of the same size (P)
    
    for (j in 1:P){
      
      #Choose parents using fitness-derived probability
      parent_1 <- parents[sample(1:P,1,prob=phi),]
      parent_2 <- parents[sample(1:P,1,prob=phi),]
      
      #Generate child and append
      child <- Cross_over_mutate(parent_1,parent_2)
      children[j,] <<- child
  
    }
    
  #This is where we would assess the fitness of the parents and children, then combine them in an
  #80:20 ratio as spoken about. At the moment the generations are not combined: the children generation becomes the next set of parents
    
  #parent_assessment <- Assess_generation(parents)
  #children_assessment <- Assess_generation(children)
  
  #children_concat <- Concatinate_generation(parents,parent_assessment,children,children_assessment)
  #children <<- children_concat[[1]]
  #children_fitness <<- children_concat[[2]] 
    
    
  # A test : here a small number of the children are randomly replaced by 
  # a small number of parents
    
  #vals <- sample(1:P,as.integer(P/8),replace=FALSE)
  #children[vals,] <<- parents[vals,]
    
  for (k in 1:P){
    
    #Select the factors we want to use in the regression
    predictor_vals <- predictors_df[,children[k,]==1]
    
    g <- lm(response_df[,1]~.,predictor_vals)
    
    fitness_value <- fitness_function(g)
    children_fitness[k] <<- fitness_value
    
    if (fitness_value < Best_fitness) {
      Best_fitness <<- fitness_value
      best_individual <<- children[k,]
    }
  }
  
  Best_fitness_generation[i+1] <<- Best_fitness
  ordered_fitness <- rank(-children_fitness)
  #Probability of selecting from a generation - proportional to the
  #fitness of the individals
  phi <<- 2*ordered_fitness/(P*(P+1))
  parents <<- children
  print(best_individual)
  print(Best_fitness)

  }
}

#------------------------------------------------------------
#Apply algorithm
#------------------------------------------------------------

Genetic_Alg_Fit <- function(dataset,iterations=60,userfunc=FALSE,mutate_rate=FALSE) {
  
  #Main function call. Should have more user-related inputs
  
  #Extract the predictors and resonse as dataframes
  subsets <- extract_response_variable(baseball_dat,"salary")
  response_df <<- log(subsets[[1]])
  predictors_df <<- subsets[[2]]
  
  C <- length(predictors)
  
  #Create data structures in the global environment
  Set_up_data_structures(C,iterations=iterations, mutate_rate=mutate_rate)
  
  #Fill the parents matrix with the first set of individuals
  Generate_first_generation()
  #Do the algorithm loop
  mainloop()
  
}

#------------------------------------------------------------
#Apply algorithm
#------------------------------------------------------------

#Read the dataset
baseball_dat = read.table(file.choose(),header=TRUE)

#Run the algorithm
Genetic_Alg_Fit(baseball_dat,mutate_rate=FALSE)


#------------------------------------------------------------
#Plot and test
#------------------------------------------------------------

plot(-Best_fitness_generation,xlim=c(0,iter),ylim=c(400,425),ylab="Negative AIC",xlab="Generation",main="Max AIC Values For Genetic Algorithm")


library('hydroGOF')
#Actually select the features of interest
run_vars <- predictors_df[,best_individual==1]
g <- lm(response_df[,1]~.,run_vars)
predicted_response <- predict(g)

#Get the rms error in terms of the actual salary
root_mean_square_error <- rmse(exp(response),exp(predicted_response))
#Not really sure how 'good' this is
