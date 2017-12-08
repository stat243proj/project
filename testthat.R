####test that

baseball.dat = read.table(file.choose(), header=TRUE)

test_all <- function(dataset, response.name){
  require(testthat)

  predictors <- dataset
  C <- length(predictors)
  Niter <<- Niter
  P <<- as.integer(C*1.5)
  prob.mutate <- 1.0/(P*sqrt(C))
  generation.old <- lapply(1:P, function(x) {rbinom(C,1,0.5)})
  fitness <- matrix(0,P,Niter)
  model.out <- lm(response[,1]~., predictors.individual)
  fitness.value <- extractAIC(model.out)[2]
  
  #test the Breed function
  test_that("Check the dimension of return list from Breed function is P embeded lists, 
            with length of each element(child genome) is C",{
              expect_equal(length(Breed(generation.old, fitness[,1], predictors, prob.mutate)), P/2)
              expect_equal(unique(unlist(
                lapply(Breed(generation.old, fitness[,1], predictors, prob.mutate),function(x) length(unlist(x))
              ))), 2*C)
              })
  #test CrossOverMutate function
  test_that("Check the dimension of return value(a child genome) from Cross_over_mutate function is C",{
    expect_equal(length(unlist(CrossOverMutate(generation.old, c(1,2), prob.mutate))), 2*C)
  })
}


test_all(baseball.dat, response.name = "salary")

#test failed *undefined columns selected
#test assessfitness function
test_that("Check type of return value of AssessFitness function",{
    expect_equal(class(
      AssessFitness(individual=rep(1, dim(baseball.dat)[2]), response=baseball.dat["salary"], user.family="gaussian", 
                    predictors=baseball.dat[, !names(baseball.dat) %in% c("salary")], userfunc=FALSE)
    ), "numeric")
  })

#test ReplaceClones function
  test_that("Check the dimension of return generation from ReplaceClones function is P",{
    expect_equal(length(ReplaceClones(generation.old, fitness.value, C)$generation), P)
  })
  test_that("Check the return fitness value of generation.old approximately equals to -210",{
    expect_equal(round(ReplaceClones(generation.old, fitness.value, C)$fitness), -210)
  })