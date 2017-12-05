---
  title: "Generic Algorithm"
author: "Kexin Fei"
date: "12/3/2017"
output: pdf_document
---
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
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

fitness_function <- function(model, userfunc=FALSE){
  
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
P <- as.integer(C*1.5) #number of individuals in a given generation
Ngen <- 100 #number of generation iterations to carry out
mutate_rate <- 1.0/(P*sqrt(C)) #mutation rate (should be about 1%)
generation_old <- matrix(0,P,C) # matrix of the factor choices (P rows, m cols)
generation_new <- matrix(0,P,C) #matrix of the improved factor choices
fitness <- matrix(0,P,Ngen) #evolution of the fitness values over model run
# best_fitness <- 0 #best fitness value 
best_fitness_generation = rep(0,Ngen) #evolution of best fitness values


#Starting generation : generate P random combinations of the factors, do the regression and record their fitness values
for (i in 1:P){
  
  #Geneate a vector of 1s and 0s corresponding to which predictors we want to use
  generation_old[i,] <- rbinom(C,1,0.5)
  
  #Select the factors we want to use in the regression
  predictors_chosen <- predictors[, generation_old[i,]==1]
  
  #linear regression of the response variable against the predictor   #data
  
  # Note what this means
  # lm( Y ~ A + B + C + D etc, data = run_vars)
  # where A,B,C,D etc are the names of columns within the dataframe
  # run_vars. The syntax here is just a more concise way or writing
  # it
  
  #model.out = lm(response[,1]~., predictors_chosen)
  g = glm(response[,1]~.,family=gaussian,data=predictors_chosen)
  fitness[i,1] <- extractAIC(g)[2]
} 

# examine initial fitness (debugging)
print(fitness[,1])
```



```{r mutate_crossover}

mutate_crossover = function(parent_1, parent_2){
  ###Input: Two parents 
  ###Output: One child after mutation and cross over
  
  #1. cross over(do not split at specific location, but randomly select 1/2 variable from a parent, and another half is from another parent)
  parent1_index = sample(c(1:C), 0.5*C, replace = FALSE)
  parent1_index_boolean = ifelse(c(1:C) %in% parent1_index, 1, 0)
  
  #2. make about 1% mutation on each children
  mutate = rbinom(C, 1, mutate_rate)
  
  child = rep(NA, C)
  for (i in c(1:C)){
    #cross over of 2 parents to generate a child
    child[i] = as.logical(ifelse(parent1_index_boolean[i]==1, parent_1[i], parent_2[i]))
    #mutation on the child
    child[i] = as.numeric(ifelse(mutate[i]==1, !child[i], child[i]))
  }
  return(child)
}
```


```{r }
#Assess fitness
Assess <- function(generation, x){
  predictors_chosen <- predictors[, generation[x,]==1]
  g <- glm(response[,1]~., family=gaussian, data=predictors_chosen)
  fitness_value <- extractAIC(g)[2]
  return(fitness_value)
}

assess_fitness <- function(generation){
  x = dim(generation)[1]
  output = sapply(c(1:P), FUN=function(x){Assess(generation, x)})
  return(output)
}

```


```{r concatenation}
concatenation <- function(ordered_fitness_old, fitness_values_new, replace_rate = 0.1){
  #rate is the replace rate when generating new parent generation.
  #The default value is 10%, but user is able to change it.
  ordered_fitness_new <- rank(-fitness_values_new)
  
  #number of individuals that will be replaced
  num_sub <- round(length(fitness_values_old)*replace_rate)
  #obtain rows of worst 10% percent of parents and best 10% percent of children
  old_index = (ordered_fitness_old <= (length(ordered_fitness_old) - num_sub))
  remain_old <- generation_old[old_index,]
  remain_new <- generation_new[(ordered_fitness_new <= num_sub), ]
  #create new parent generation
  generation_old <- rbind(remain_old, remain_new)  
  
  #return new parent generation
  return(generation_old)
}
```


#Main Function
```{r main loop}
#?????iter-1
for (i in 1:Ngen){
  
  #Do the crossover and mutation stage, obtain a new generation
  fitness_values_old = assess_fitness(generation_old)
  ordered_fitness_old = rank(-fitness_values_old)
  
  phi = 2*ordered_fitness_old/(P*(P+1))
  
  for (j in 1:P){
    parent_1 = generation_old[sample(x=1:P, size=1, prob=phi), ]
    parent_2 = generation_old[sample(x=1:P, size=1, prob=phi), ]
    generation_new[j,] = mutate_crossover(parent_1,parent_2)
  }
  
  fitness_values_new = assess_fitness(generation_new)
  
  generation_old = concatenation(ordered_fitness_old, fitness_values_new, replace_rate = 0.1)
}  
```
#best fitness
for (j in 1:P){
  predictors_chosen <- predictors[,generation_old[j,]==1]
  
  g = glm(response[,1]~.,family=gaussian,data=predictors_chosen)
  
  fitness_value <- extractAIC(g)[2]
  
  fitness[j,i] <- fitness_value
  
  best_fitness <- min(fitness)
}

#plot
plot(-fitness,xlim=c(0,Ngen),ylim=c(50,425),type="n",ylab="Negative AIC",
     xlab="Generation",main="AIC Values For Genetic Algorithm")
for(i in 1:Ngen){points(rep(i,P),-fitness[,i],pch=20)}
