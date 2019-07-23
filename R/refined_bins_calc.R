
# this is a code to aggregate the bins
# should generate a new binspec and generate lower and upper bounds
# this new binspec contains only from and to 
# lower and upper bounds are the same as those generated from bounds_from_param

automatic_bins <- function(pre_calc_hist, threshold_ratio = 0.005){
  a <- aggregate(pre_calc_hist$counts, list(upper = pre_calc_hist$upper),sum)
  sumbefore <- sum(a$x)
  # to fulfill Randy's comments on fully automating this 
  # I chose to change the threshold input from a value to a ratio
  # ratio = percentage of the total counts in all histograms
  # say we have 2000+ grids, total counts ~= 221M
  # We then let the ratio = 0.01, this means that the threshold_value = 2.2M
  # this strategy is a little better because it doesn't depend on the number of input grids 
  # add a error message to warn the user
  if (threshold_ratio >= 1) {
    stop("Error: threshold_ratio cannot be greater than or equal to 1! Do you want just one bin?")
  }
  threshold_value = threshold_ratio * sumbefore
  # test if threshold 2 works
  a_prime = a
  counter = 0
  while (any(a_prime$x < threshold_value)) {
    idx <- vector() # declare an empty vector for reducing row purposes
    counter = counter + 1
    for (row in (1:nrow(a_prime))){
      if (a_prime$x[row] < threshold_value) {
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
    stop("Error: before and after doesn't equal")
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