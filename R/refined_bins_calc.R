
# this is a code to aggregate the bins
# should generate a new binspec and generate lower and upper bounds
# this new binspec contains only from and to 
# lower and upper bounds are the same as those generated from bounds_from_param

automatic_bins <- function(pre_calc_hist, threshold = 100){
  a <- aggregate(pre_calc_hist$counts, list(upper = pre_calc_hist$upper),sum)
  sumbefore <- sum(a$x)
  # test if threshold 2 works
  a_prime = a
  counter = 0
  while (any(a_prime$x < threshold)) {
    idx <- vector() # declare an empty vector for reducing row purposes
    counter = counter + 1
    for (row in (1:nrow(a_prime))){
      if (a_prime$x[row] < threshold) {
        # combine the value to the value in the next row
        a_prime$x[row+1] <- a_prime$x[row] + a_prime$x[row+1]
        # delete that row
        idx = c(idx, row)
      }
    }
    # after looping over, you have a vector of indices, use these indices to reduce the rows
    a_prime <- a_prime[-idx, ]
  }
  sumafter <- sum(a_prime$x)
  if (sumbefore != sumafter) {
    stop("before and after doesn't equal")
  }
  # compute the lower bounds
  lower <- vector() # declare an empty vector
  for (row in (1:nrow(a_prime))){
    if (row == 1) {
      lower <- c(lower, pre_calc_hist$lower[1])
    } else{
    lower <- c(lower, a_prime$upper[row-1])
    }
  }
  # bind it to a_prime dataframe
  a_prime$lower <- lower
  # reorder to put lower bound in the front
  a_prime <- a_prime[, c(3,1,2)]
  # returns upper and lower bounds
  list(lower=a_prime$lower, upper=a_prime$upper)
}
# I think it works
# this code should return two arrays of upper and lower bounds