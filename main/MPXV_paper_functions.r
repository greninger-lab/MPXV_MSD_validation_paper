# Geometric mean functions


gmean <- function(x, remove_na = FALSE) { 
  result <- 10^mean(log10(x), na.rm = remove_na) 
}


gsd <- function(x, remove_na = FALSE) { 
  result <- 10^sd(log10(x), na.rm = remove_na) 
  
}

gcv <- function(x, remove_na = FALSE) { 
  result <- sqrt(exp(sd(log(x), na.rm = remove_na)^2) - 1)
  
}

# Function for loading packages

check_and_load <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
