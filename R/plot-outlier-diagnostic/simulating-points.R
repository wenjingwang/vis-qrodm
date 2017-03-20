## Simulating points in space

simvar <- function(x, n = 10, method = "grid") UseMethod("simvar")
simvar.factor <- function(x, n = 10, method = "grid"){
  switch(method,
         random = x[sample(length(x), n, replace = TRUE)],
         factor(levels(x), levels = levels(x))
  )
}

simvar.numeric <- function(x, n = 10, method = "grid"){
  rng <- range(x)
  switch(method,
         random = runif(n, rng[1], rng[2]),
         seq(rng[1], rng[2], length = n))
}

generate_data <- function(data, n = 1000, method = "grid"){
  if(method != "random"){
    n <- floor(n ^ (1/ncol(data)))
    df <- data.frame(expand.grid(lapply(data, simvar, n = n,
                                        method = "grid")))
    if(method == "nonaligned"){
      cont <- !sapply(df, is.factor)
      ranges <- lapply(df[, cont], function(x) diff(range(x)))
      df[,cont] <- df[,cont] +
        do.call(cbind, lapply(ranges, function(rng)
          runif(-rng/(2*n), rng/(2*n), n=nrow(df))))
    }
    df
  }else{
    data.frame(sapply(data, simvar, n=n, method=method))
  }
}