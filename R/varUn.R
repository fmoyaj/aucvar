# Kendall Tau example -> provide U-stat?????


uStat <- function(sample1, sample2, n1, n2, k1, k2, phi){
  stat <-  1/(VeryLargeIntegers::binom(n1,k1)*VeryLargeIntegers::binom(n2,k2))
  sum <- 0

  s1_label <- utils::combn(n1,k1)
  s2_label <- utils::combn(n2,k2)

  # for loop

  return (stat)
}


exactVar <- function(sample1, sample2, n1, n2, k1, k2, phi){

  N0 <- VeryLargeIntegers::binom(n1, k1)*VeryLargeIntegers::binom(n1-k1,k1)*
    VeryLargeIntegers::binom(n2,k2)*VeryLargeIntegers::binom(n2-k2,k2)

  if (N0 > 10^5){
    stop("The exhaustive number of non-overlapping pairs of data subsets is greater than 10^5.
         It is computationally expensive to compute the exact U-statistic variance.")
  }

  Qk <- uStat(sample1, sample2, n1, n2, k1, k2, phi)^2

  Q0 <- 0
  id_n1_choose_2k1 <- utils::combn(n1,2*k1)
  id_n2_choose_2k2 <- utils::combn(n2,2*k2)
  id_2k1_choose_k1 <- utils::combn(2*k1,k1)
  id_2k2_choose_k2 <- utils::combn(2*k2,k2)

  for(i in 1:VeryLargeIntegers::binom(n1,2*k1))
  {
    for(j in 1:VeryLargeIntegers::binom(n2,2*k2))
    {
      sample_2k1 <- sample1[id_n1_choose_2k1[i,]]
      sample_2k2 <- sample2[id_n2_choose_2k2[j,]]
      for(s in 1:VeryLargeIntegers::binom(2*k1,k1))
      {
        for(t in 1:VeryLargeIntegers::binom(2*k2,k2))
        {
          Q0 <- Q0 +(1/N0)*phi(c(sample_2k1[id_2k1_choose_k1[s,]], sample_2k2[id_2k2_choose_k2[t,]]) ) * phi(c(sample_2k1[-id_2k1_choose_k1[s,]], sample_2k2[-id_2k2_choose_k2[t,]]) )
        }
      }
    }
  }

  result <- Qk - Q0
  return (result)
}



#' Compute the variance of a general 2-sample U-statistic
#'
#' @param sample1 The first sample.
#' @param sample2 The second sample.
#' @param k1 The number of observations from sample 1.
#' @param k2 The number of observations from sample 2.
#' @param phi The kernel function for the U-statistic
#' @param B The number of random partitions.
#'
#' @return The variance of the 2-sample U-statistic
#' @export
#'
#' @examples
#' 'TODO'


varUn <- function(sample1, sample2, k1, k2, phi , B = Inf){

  # Check parameters
  if (typeof(phi) != "closure"){
    stop("the argument 'phi' is not a function.")
  }

  # Length of array
  n1 <- base::nrow(sample1)
  n2 <- base::nrow(sample2)

  # Length of vector
  if (base::is.null(n1) || base::is.null(n2)){
    n1 <- base::length(sample1)
    n2 <- base::length(sample2)
  }

  if (k1 > n1/2 || k2 > n2/2){
    stop("k1, k2 cannot be greater than half of the length of the corresponding sample n1 or n2")
  }


  # Maximum number of partitions we can get
  m <- base::min(base::floor(n1/k1), base::floor(n2/k2))

  for (b in 1:B){
    # For each sampled partition
    kernel <-base::array(NA,m)

    id1.sub <- base::array(NA,m)
    id2.sub <- base::array(NA,m)

    # Shuffle each of the samples
    id1.shuffle <- sample(1:n1, n1, replace=F)
    id2.shuffle <- sample(1:n2, n2, replace=F)


    # Calculate the kernel function for each block of paired data subsets
    # Based on the shuffled samples, can consider blocks of size k1 (or k2) consecutively as the subsample(s) in the kernel function
    for (block in 1:m)
    {
      kernel[block] <- phi(sample1[id1.shuffle[((block-1)*k1):(block*k1)], ], sample2[id2.shuffle[((block-1)*k2):(block*k2)], ])
    }

    phi_bar_b[b] <- base::mean(kernel)


    # Compute the within-partition sum of squares
    WPSS[b] <- (1/(m-1)) * base::sum(kernel-phi_bar_b[b])^2
  }

  # Compute the between-partition sum of squares
  BPSS <- (1/B) * base::mean((phi_bar_b - mean(phi_bar_b))^2)

  # Return the partition-resampling variation of the AUC variance
  return (1/m)*base::mean(WPSS) - BPSS
}
