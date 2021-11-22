# This function computes the posterior distribution of abilities for examinees
# given their item response data, item parameters, and population distribution
posterior <- function(likehd, weights, idx.std = NULL) {
  
  # check if MG-calibration is implemented
  if(is.null(idx.std)) {
    
    # count the number of quad points
    ntheta <- nrow(weights)
    
  } else {
    
    # count the number of quad points
    ntheta <- nrow(weights[[1]])
    
    # divide the likelihood matrix into several groups if idx.std is not NULL
    # this is only for MG-calibration
    likehd.gr <- purrr::map(.x = idx.std, ~{likehd[.x, ]})      
    
  }
  
  # create the joint likelihood matrix of likelihoods and population distribution
  joint_like <- array(0, c(nrow(likehd), ntheta))
  if(is.null(idx.std)) {
    for(i in 1:ntheta) {
      joint_like[, i] <- likehd[, i] * weights[i, 2]
    }
  } else {
    for(g in 1:length(weights)) {
      for(i in 1:ntheta) {
        joint_like[idx.std[[g]], i] <- likehd.gr[[g]][, i] * weights[[g]][i, 2]
      }
    }
  }
  
  # denominator of the posterior distribution
  denom <- rowSums(joint_like)
  
  # compute the posterior distribution of examinees across the node values
  # a row and column indicate the examinees and nodes, respectively
  posterior <- joint_like / denom
  
  # return results
  posterior
  
}


