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

#Some other notes:

#Cross over could be replaced with random sampling without repacement from one individual and then combine with 
#half the genes from the other parent. 
#This creates P children, use tornament matching so that we keep some of the parents (helps to reduce the likelihood of children being worse than the parents)
#Allow mutation - user defined or 1/m
#Criteria for convergence (min, max and varience of the AIC). Keep track of the variance of generation i and i-1 and then stop when this difference is small enough

