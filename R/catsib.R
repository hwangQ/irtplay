#' CATSIB DIF detection procedure
#'   
#' @description This function analyze DIF on items using CATSIB procedure (Nandakumar & Roussos, 2004), which is a modified 
#' version of SIBTEST (Shealy & Stout, 1993).
#'   
#' @param x A data frame containing the item metadata (e.g., item parameters, number of categories, models ...). 
#' This metadata is required to estimate latent ability parameters when \code{score = NULL} or \code{purify = TRUE}. Default is NULL. 
#' See \code{\link{est_irt}}, \code{\link{irtfit}}, \code{\link{test.info}} or \code{\link{simdat}} for more details about the item metadata.
#' @param data A matrix containing examinees' response data for the items in the argument \code{x}. A row and column indicate
#' the examinees and items, respectively. 
#' @param score A vector of examinees' ability estimates. If the abilities are not provided, \code{\link{catsib}} function computes the ability 
#' estimates before computing CATSIB statistic. See \code{\link{est_score}} for more details about scoring methods. Default is NULL. 
#' @param se A vector of the standard errors of the ability estimates. The standard errors should be ordered in accordance with the order of 
#' the ability estimates specified in the \code{score} argument. Default is NULL.
#' @param group A numeric or character vector indicating group membership of examinees. The length of vector should the same as the number of rows 
#' in the response data matrix. 
#' @param focal.name A single numeric or character indicating the level of group which corresponds to the focal group. 
#' For example, if \code{group = c(0, 1, 0, 1, 1)} and '1' indicates the focal group, then \code{focal.name = 1}. 
#' @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
#' Default is 1.
#' @param n.bin A vector of two positive integer values to set the maximum and minimum numbers of bins (or intervals) on the ability scale. 
#' The first and second values indicates the maximum and minimum numbers, respectively. the \code{\link{catsib()}} function automatically starts
#' with the maximum number of the bins and to then monitor how many examinees would be eliminated due to the throwing out of sparse cells. 
#' If more than 7.5% of either the reference or focal group examinees are eliminated, the function decreases the number of cells until the number of
#' examinees eliminated from each group becomes less than or equal to 7.5%. Even if the number of examinees eliminated from one or both of the groups
#' exceeded 7.5%, the number of bins cannot be lower than the minimum number of bins. Default is c(80, 10). 
#' @param min.binsize A positive integer value to set the minimum size of each bin. To ensure stable statistical estimation, each bin is required 
#' to have a minimum number of certain examinees from both the reference and focal groups if it was to be included in the calculation of 
#' the CATSIB statistic. All bins with fewer than this minimum number are not used. Default is 3.
#' @param max.del A numerical value to set the maximum permissible proportion of examinees to be deleted from either the reference or focal groups
#' when automatically determining the number of bins on the ability scale. Default is 0.075.
#' @param weight.group A single character string indicating a group over which the expectation of CATSIB statistic and the corresponding standard error are computed.  
#' @param alpha A numeric value to specify significance \eqn{\alpha}-level of the hypothesis test using the RDIF fit statistics.
#' Default is .05.
#' @param missing A value indicating missing values in the response data set. Default is NA.
#' @param purify A logical value indicating whether a purification process will be implemented or not. Default is FALSE.
#' @param max.iter An positive integer value specifying the maximum number of iterations for the purification process. Default is 10. 
#' @param min.resp An positive integer value specifying the minimum number of item responses for an examinee when scores are computed. 
#' Default is NULL. See details below for more information.  
#' @param method A character string indicating a scoring method. Available methods are "MLE" for the maximum likelihood estimation, 
#' "MAP" for the maximum a posteriori estimation, and "EAP" for the expected a posteriori estimation. Default method is "MLE".
#' @param range A numeric vector of two components to restrict the range of ability scale for the MLE. Default is c(-4, 4).
#' @param norm.prior A numeric vector of two components specifying a mean and standard deviation of the normal prior distribution.
#' These two parameters are used to obtain the gaussian quadrature points and the corresponding weights from the normal distribution. Default is
#' c(0,1). Ignored if \code{method} is "MLE".
#' @param nquad An integer value specifying the number of gaussian quadrature points from the normal prior distribution. Default is 41.
#' Ignored if \code{method} is "MLE" or "MAP".
#' @param weights A two-column matrix or data frame containing the quadrature points (in the first column) and the corresponding weights
#' (in the second column) of the latent variable prior distribution. The weights and quadrature points can be easily obtained
#' using the function \code{\link{gen.weight}}. If NULL and \code{method} is "EAP", default values are used (see the arguments
#' of \code{norm.prior} and \code{nquad}). Ignored if \code{method} is "MLE" or "MAP".
#' @param ncore The number of logical CPU cores to use. Default is 1. See \code{\link{est_score}} for details.
#' @param verbose A logical value. If TRUE, the progress messages of purification procedure are suppressed. Default is TRUE.
#' 
#' @details See Nandakumar and Roussos (2004) for more detail about the CATSIB method. 
#'
#' @return This function returns a list of four internal objects. The four objects are: 
#' \item{no_purify}{A list of several sub-objects containing the results of DIF analysis without a purification procedure. The sub-objects are: 
#'     \describe{
#'       \item{dif_stat}{A data frame containing the results of CATSIB statistics across all evaluated items. From the first column, each column 
#'        indicates item's ID, CATSIB (\emph{beta}) statistic, standard error of the \emph{beta}, standardized \emph{beta}, p-value of the \emph{beta}, 
#'        sample size of the reference group, sample size of the focal group, and total sample size, respectively.} 
#'       \item{dif_item}{A numeric vector showing potential DIF items flagged by CATSIB statistic.}
#'       \item{contingency}{A contingency table of each item used to compute CATSIB statistic.}
#'    }
#' }
#' \item{purify}{A logical value indicating whether the purification process was used.} 
#' \item{with_purify}{A list of several sub-objects containing the results of DIF analysis with a purification procedure. The sub-objects are:
#'     \describe{
#'       \item{dif_stat}{A data frame containing the results of CATSIB statistics across all evaluated items. From the first column, each column 
#'        indicates item's ID, CATSIB (\emph{beta}) statistic, standard error of the \emph{beta}, standardized \emph{beta}, p-value of the \emph{beta}, 
#'        sample size of the reference group, sample size of the focal group, and total sample size, and \emph{n}th iteration where the CATSIB statistic 
#'        was computed, respectively.} 
#'       \item{dif_item}{A numeric vector showing potential DIF items flagged by CATSIB statistic.}
#'       \item{n.iter}{A total number of iterations implemented for the purification.}
#'       \item{complete}{A logical value indicating whether the purification process was completed. If FALSE, it means that the purification process 
#'        reached the maximum iteration number but it was not complete.}
#'       \item{contingency}{A contingency table of each item used to compute CATSIB statistic.}
#'     }
#' }
#' \item{alpha}{A significance \eqn{\alpha}-level used to compute the p-values of RDIF statistics.} 
#' 
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#' 
#' @seealso \code{\link{rdif}}, \code{\link{est_item}}, \code{\link{test.info}}, \code{\link{simdat}}, 
#' \code{\link{shape_df}}, \code{\link{gen.weight}}, \code{\link{est_score}}
#' 
#' @references
#' Nandakumar, R., & Roussos, L. (2004). Evaluation of the CATSIB DIF procedure in a pretest setting. 
#' \emph{Journal of Educational and Behavioral Statistics, 29}(2), 177-199.
#' Shealy, R., & Stout, W. (1993). A model-based standardization approach that separates true bias/DIF from group 
#' ability differences and detects test bias/DTF as well as item bias/DIF. \emph{Psychometrika}, 58, 159-194.
#' 
#' @examples
#' \donttest{
#' # call library
#' library("dplyr")
#' 
#' ## Uniform DIF detection 
#' ###############################################
#' # (1) manipulate true uniform DIF data
#' ###############################################
#' # import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")
#' 
#' # select 36 of 3PLM items which are non-DIF items 
#' par_nstd <- 
#'   bring.flexmirt(file=flex_sam, "par")$Group1$full_df %>% 
#'   dplyr::filter(.data$model == "3PLM") %>% 
#'   dplyr::filter(dplyr::row_number() %in% 1:36) %>% 
#'   dplyr::select(1:6)
#' par_nstd$id <- paste0("nondif", 1:36)
#' 
#' # generate four new items to inject uniform DIF
#' difpar_ref <- 
#'   shape_df(par.dc=list(a=c(0.8, 1.5, 0.8, 1.5), b=c(0.0, 0.0, -0.5, -0.5), g=0.15), 
#'            item.id=paste0("dif", 1:4), cats=2, model="3PLM")
#' 
#' # manipulate uniform DIF on the four new items by adding constants to b-parameters 
#' # for the focal group
#' difpar_foc <- 
#'   difpar_ref %>% 
#'   dplyr::mutate_at(.vars="par.2", .funs=function(x) x + rep(0.7, 4))
#' 
#' # combine the 4 DIF and 36 non-DIF items for both reference and focal groups
#' # thus, the first four items have uniform DIF 
#' par_ref <- rbind(difpar_ref, par_nstd)
#' par_foc <- rbind(difpar_foc, par_nstd)
#' 
#' # generate the true thetas
#' set.seed(123)
#' theta_ref <- rnorm(500, 0.0, 1.0)
#' theta_foc <- rnorm(500, 0.0, 1.0)
#' 
#' # generate the response data
#' resp_ref <- simdat(par_ref, theta=theta_ref, D=1)
#' resp_foc <- simdat(par_foc, theta=theta_foc, D=1)
#' data <- rbind(resp_ref, resp_foc)
#' 
#' ###############################################
#' # (2) estimate the item and ability parameters 
#' #     using the aggregate data
#' ###############################################
#' # estimate the item parameters 
#' est_mod <- est_irt(data=data, D=1, model="3PLM")
#' est_par <- est_mod$par.est
#' 
#' # estimate the ability parameters using MLE
#' theta_est <- est_score(x=est_par, data=data, method="MLE")
#' score <- theta_est$est.theta
#' se <- theta_est$se.theta
#' 
#' ###############################################
#' # (3) conduct DIF analysis
#' ###############################################
#' # create a vector of group membership indicators 
#' # where '1' indicates the focal group 
#' group <- c(rep(0, 500), rep(1, 500))
#' 
#' # (a)-1 compute SIBTEST statistic by providing scores, 
#' #       and without a purification 
#' dif_1 <- catsib(x=NULL, data=data, D=1, score=score, se=se, group=group, focal.name=1,
#'  weight.group="comb", alpha=0.05, missing=NA, purify=FALSE)
#' print(dif_1)
#' 
#' # (a)-2 compute SIBTEST statistic by providing scores, 
#' #       and with a purification 
#' dif_2 <- catsib(x=est_par, data=data, D=1, score=score, se=se, group=group, focal.name=1,
#'  weight.group="comb", alpha=0.05, missing=NA, purify=TRUE)
#' dif_2$no_purify$dif_stat
#' 
#' }
#' 
#' @export
#' 
catsib <- function(x=NULL, data, score=NULL, se=NULL, group, focal.name, D=1, n.bin=c(80, 10), min.binsize=3, max.del=0.075, 
                   weight.group = c("foc", "ref", "comb"), alpha=0.05, missing=NA, purify=FALSE, max.iter=10, min.resp=NULL, method="MLE", 
                   range=c(-4, 4), norm.prior=c(0, 1), nquad=41, weights=NULL, ncore=1, verbose=TRUE) {
  
  # match.call
  cl <- match.call()
  
  ##----------------------------------
  ## (1) prepare DIF analysis
  ##----------------------------------
  # when purify = TRUE, the item metadata should be provided. 
  # if not, stop. 
  if(purify & is.null(x)) {
    stop("To implement a purification process, the item metadata must be provided in the argument 'x'.", call.=FALSE)
  }
  
  # clean the data frame of the item metadata 
  if(!is.null(x)) {
    
    # give column names
    x <- data.frame(x)
    colnames(x) <- c("id", "cats", "model", paste0("par.", 1:(ncol(x) - 3)))
    
    # add par.3 column when there is no par.3 column (just in case that all items are 2PLMs)
    if(ncol(x[, -c(1, 2, 3)]) == 2) {
      x <- data.frame(x, par.3=NA)
    }
    
    # clear the item metadata set
    x <- back2df(metalist2(x))
    
    # change the model to the character string 
    x$model <- as.character(x$model)
    
    # stop when the model includes any polytomous model
    if(any(x$model %in% c("GRM", "GPCM")) | any(x$cats > 2)) {
      stop("The current version only supports dichotomous response data.", call.=FALSE)
    }
    
  }
  
  # transform the response data to a matrix form
  data <- data.matrix(data)
  
  # re-code missing values
  if(!is.na(missing)) {
    data[data == missing] <- NA
  }
  
  # stop when the model includes any polytomous response data
  if(any(data > 1, na.rm=TRUE)) {
    stop("The current version only supports dichotomous response data.", call.=FALSE)
  }
  
  # create an item id 
  if(is.null(x)) {
    item.id <- paste0("item.", 1:ncol(data))
  } else {
    item.id <- x$id
  }
  
  # compute the score if score = NULL
  if(!is.null(score)) {
    # transform scores to a vector form
    if(is.matrix(score) | is.data.frame(score)) {
      score <- as.numeric(data.matrix(score))
    }
    # transform scores to a vector form
    if(is.matrix(se) | is.data.frame(se)) {
      score <- as.numeric(data.matrix(se))
    }
  } else {
    
    # if min.resp is not NULL, find the examinees who have the number of responses
    # less than specified value (e.g., 5). Then, replace their all responses with NA
    if(!is.null(min.resp)) {
      n_resp <- rowSums(!is.na(data))
      loc_less <- which(n_resp < min.resp & n_resp > 0)
      data[loc_less, ] <- NA
    }
    score_rst <- est_score(x=x, data=data, D=D, method=method, range=range, norm.prior=norm.prior, 
                           nquad=nquad, weights=weights, ncore=ncore)
    score <- score_rst$est.theta
    se <- score_rst$se.theta
  }
  
  # confirm the target ability density function
  weight.group <- match.arg(weight.group)
  weight.group <- switch(weight.group, 
                         foc = "foc", 
                         ref = "ref", 
                         comb = "comb")
  
  # set the maximum & minimum number of bins (intervals)
  max.bin <- n.bin[1]
  min.bin <- n.bin[2]
  
  # a) when no purification is set
  # do only one iteration of DIF analysis:
  # compute the beta statistic and its SE for all items using
  # corrected theta scores
  dif_rst <- 
    catsib_one(data=data, group=group, focal.name=focal.name, score=score, se=se, range=range,
               item.id=item.id, max.bin=max.bin, min.bin=min.bin, min.binsize=min.binsize, 
               max.del=max.del, weight.group=weight.group, alpha=alpha)
  
  # create two empty lists to contain the results
  no_purify <- list(dif_stat=NULL, dif_item=NULL, contingency=NULL)
  with_purify <- list(dif_stat=NULL, dif_item=NULL, n.iter=NULL, 
                      complete=NULL, contingency=NULL)
  
  # record the first DIF detection results into the no purification list  
  no_purify$dif_stat <- dif_rst$dif_stat
  no_purify$dif_item <- dif_rst$dif_item
  no_purify$contingency <- dif_rst$contingency
  
  # when purification is used
  if(purify) {
    
    # create an empty vector and empty data frames 
    # to contain the detected DIF items, statistics, and contingency tables
    dif_item <- NULL
    dif_stat <- 
      data.frame(id=rep(NA_character_, nrow(x)), beta=NA, se=NA,
                 z.beta=NA, p.val=NA, n.ref=NA, n.foc=NA, n.total=NA, n.iter=NA, 
                 stringsAsFactors=FALSE)
    contingency <- vector('list', nrow(x))
    names(contingency) <- item.id 
    
    # extract the first DIF analysis results
    # and check if at least one DIF item is detected
    dif_item_tmp <- dif_rst$dif_item
    dif_stat_tmp <- dif_rst$dif_stat
    contingency_tmp <- dif_rst$contingency
    
    # copy the response data and item meta data
    x_puri <- x
    data_puri <- data
    
    # start the iteration if any item is detected as an DIF item
    if(!is.null(dif_item_tmp)) {
      
      # record unique item numbers
      item_num <- 1:nrow(x)
      
      # in case when at least one DIF item is detected from the no purification DIF analysis
      # in this case, the maximum number of iteration must be greater than 0.   
      # if not, stop and return an error message
      if(max.iter < 1) stop("The maximum iteration (i.e., max.iter) must be greater than 0 when purify = TRUE.", call.=FALSE)
      
      # print a message
      if(verbose) {
        cat("Purification started...", '\n')
      }
      
      for(i in 1:max.iter) {
        
        # print a message
        if(verbose) {
          cat("\r", paste0("Iteration: ", i))
        }
        
        # a flagged item which has the largest significant DIF statistic
        flag_max <- which.max(abs(dif_stat_tmp$z.beta)) 
        
        # check an item that is deleted
        del_item <- item_num[flag_max]
        
        # add the deleted item as the DIF item
        dif_item <- c(dif_item, del_item)
        
        # add the DIF statistics and moments for the detected DIF item
        dif_stat[del_item, 1:8] <- dif_stat_tmp[flag_max, ]
        dif_stat[del_item, 9] <- i - 1
        contingency[del_item] <- contingency_tmp[flag_max]
        
        # refine the leftover items
        item_num <- item_num[-flag_max]
        
        # remove the detected DIF item data which has the largest statistic from the item metadata
        x_puri <- x_puri[-flag_max, ]
        
        # remove the detected DIF item data which has the largest statistic from the response data
        data_puri <- data_puri[, -flag_max]
        
        # if min.resp is not NULL, find the examinees who have the number of responses
        # less than specified value (e.g., 5). Then, replace their all responses with NA
        if(!is.null(min.resp)) {
          n_resp <- rowSums(!is.na(data_puri))
          loc_less <- which(n_resp < min.resp & n_resp > 0)
          data_puri[loc_less, ] <- NA
        }
        
        # compute the updated ability estimates after deleting the detected DIF item data
        score_rst_puri <- est_score(x=x_puri, data=data_puri, D=D, method=method, range=range, norm.prior=norm.prior, 
                                    nquad=nquad, weights=weights, ncore=ncore)
        score_puri <- score_rst_puri$est.theta
        se_puri <- score_rst_puri$se.theta
        
        # do DIF analysis using the updated ability estimates
        item.id <- x_puri$id
        dif_rst_tmp <- 
          catsib_one(data=data_puri, group=group, focal.name=focal.name, score=score_puri, se=se_puri, range=range,
                     item.id=item.id, max.bin=max.bin, min.bin=min.bin, min.binsize=min.binsize, 
                     max.del=max.del, weight.group=weight.group, alpha=alpha)
        
        # extract the first DIF analysis results
        # and check if at least one DIF item is detected
        dif_item_tmp <- dif_rst_tmp$dif_item
        dif_stat_tmp <- dif_rst_tmp$dif_stat
        contingency_tmp <- dif_rst_tmp$contingency
        
        # check if a further DIF item is flagged
        if(is.null(dif_item_tmp)) {
          
          # add no additional DIF item
          dif_item <- dif_item
          
          # add the DIF statistics for rest of items
          dif_stat[item_num, 1:8] <- dif_stat_tmp
          dif_stat[item_num, 9] <- i
          contingency[item_num] <- contingency_tmp
          
          break
          
        } 
      }
      
      # print a message
      if(verbose) {
        cat("", "\n")
      }
      
      # record the actual number of iteration 
      n_iter <- i
      
      # if the iteration reached out the maximum number of iteration but the purification was incomplete, 
      # then, return a warning message
      if(max.iter == n_iter & !is.null(dif_item_tmp)) {
        warning("The iteration reached out the maximum number of iteration before purification is complete.", call.=FALSE)
        complete <- FALSE
        
        # add flagged DIF item at the last iteration
        dif_item <- c(dif_item, item_num[dif_item_tmp])
        
        # add the DIF statistics for rest of items
        dif_stat[item_num, 1:8] <- dif_stat_tmp
        dif_stat[item_num, 9] <- i
        contingency[item_num] <- contingency_tmp
        
      } else {
        complete <- TRUE
        
        # print a message
        if(verbose) {
          cat("Purification is finished.", '\n')
        }
      }
      
      # record the final DIF detection results with the purification procedure
      with_purify$dif_stat <- dif_stat
      with_purify$dif_item <- sort(dif_item)
      with_purify$n.iter <- n_iter
      with_purify$complete <- complete
      with_purify$contingency <- contingency
      
      
    } else {
      
      # in case when no DIF item is detected from the first DIF analysis results
      with_purify$dif_stat <- cbind(no_purify$dif_stat, n.iter=0)
      with_purify$n.iter <- 0
      with_purify$complete <- TRUE
      with_purify$contingency <- no_purify$contingency
      
    }
    
  }
  
  
  # summarize the results
  rst <- list(no_purify=no_purify, purify=purify, with_purify=with_purify, alpha=alpha) 
  
  # return the DIF detection results
  class(rst) <- "catsib"
  rst$call <- cl
  rst
  
}


# This function performs a regression correction for ability estimates and, then
# computes the beta statistic and its SE for all items 
catsib_one <- function(data, group, focal.name, score, se, range,
                       item.id, max.bin, min.bin, min.binsize=3, max.del=0.075, 
                       weight.group, alpha=0.05) {
  
  ##------------------------------------------
  ## prepare data sets
  ##------------------------------------------
  # count the number of items
  nitem <- length(item.id)
  
  # find the location of examinees for the reference and the focal groups 
  loc_ref <- which(group != focal.name)
  loc_foc <- which(group == focal.name)
  
  # divide the response data into the two group data 
  resp_ref <- data[loc_ref, , drop=FALSE]
  resp_foc <- data[loc_foc, , drop=FALSE]
  
  # divide the thetas into the two group data 
  score_ref <- score[loc_ref]
  score_foc <- score[loc_foc]
  
  # divide the SEs into the two group data 
  se_ref <- se[loc_ref]
  se_foc <- se[loc_foc]
  
  ##------------------------------------------
  ## a regression correction for theta scores
  ##------------------------------------------
  # compute the mean and variance of scores for each group
  # note that, the bound scores are excluded.
  lb_score <- range[1]
  up_score <- range[2]
  # mu_ref <- mean(score_ref, na.rm=TRUE)
  # mu_foc <- mean(score_foc, na.rm=TRUE)
  # sigma2_ref <- stats::var(score_ref, na.rm=TRUE)
  # sigma2_foc <- stats::var(score_foc, na.rm=TRUE)
  mu_ref <- mean(score_ref[score_ref > lb_score & score_ref < up_score], na.rm=TRUE)
  mu_foc <- mean(score_foc[score_foc > lb_score & score_foc < up_score], na.rm=TRUE)
  sigma2_ref <- stats::var(score_ref[score_ref > lb_score & score_ref < up_score], na.rm=TRUE)
  sigma2_foc <- stats::var(score_foc[score_foc > lb_score & score_foc < up_score], na.rm=TRUE)
  
  # compute the error variance of scores for each group
  errvar_ref <- stats::var(se_ref[score_ref > lb_score & score_ref < up_score], na.rm=TRUE)
  errvar_foc <- stats::var(se_foc[score_foc > lb_score & score_foc < up_score], na.rm=TRUE)
  
  # compute the squared correlation (a.k.a. reliability) between theta estimate and true theta 
  # rho_ref <- suppressWarnings(sqrt(1 - se_ref^2 / sigma2_ref))
  # rho_foc <- suppressWarnings(sqrt(1 - se_foc^2 / sigma2_foc))
  rho_ref2 <- suppressWarnings(1 - errvar_ref / sigma2_ref)
  rho_foc2 <- suppressWarnings(1 - errvar_foc / sigma2_foc)
  
  # apply a regression correction to the ability estimates
  crscore_ref <- mu_ref + rho_ref2 * (score_ref - mu_ref)
  crscore_foc <- mu_foc + rho_foc2 * (score_foc - mu_foc)
  
  ##------------------------------------------
  ## compute CATSIB statistic
  ##------------------------------------------
  # conduct a DIF analysis
  catsib_dif <- purrr::map(.x=1:nitem,
                           .f=function(i) {
                             catsib_item(crscore_ref=crscore_ref, crscore_foc=crscore_foc, 
                                         resp.ref=resp_ref[, i], resp.foc=resp_foc[, i], 
                                         max.bin=max.bin, min.bin=min.bin, min.binsize=min.binsize, 
                                         max.del=max.del, weight.group=weight.group)
                           })
  
  # extract the DIF analysis results
  dif_stat <- 
    purrr::map(.x=catsib_dif, "dif_stat") %>% 
    dplyr::bind_rows() %>% 
    dplyr::mutate_if("is.numeric", "round", digits=4) %>% 
    dplyr::mutate(id=item.id) %>% 
    dplyr::relocate(.data$id, .before="beta")
  contingency <- purrr::map(.x=catsib_dif, "contingency")
  names(contingency) <- item.id
  dif_item <- as.numeric(which(dif_stat$p.val <= alpha))
  if(length(dif_item) == 0) dif_item <- NULL
  
  # return the results
  rst <- list(dif_stat=dif_stat, dif_item=dif_item, contingency=contingency)
  rst
  
}


# This function computes the beta statistic and its SE for an item
catsib_item <- function(crscore_ref, crscore_foc, resp.ref, resp.foc, 
                        max.bin, min.bin, min.binsize=3, max.del=0.075,
                        weight.group) {
  
  # combine all corrected theta scores 
  crscore <- c(crscore_ref, crscore_foc)
  
  # set the range of the ability scale
  min.crscore <- min(crscore, na.rm=TRUE)
  max.crscore <- max(crscore, na.rm=TRUE)
  
  # decide the number of bins and create an initial frequency table
  for(num.bin in max.bin:min.bin) {
    
    # compute the cut-scores to divide the theta scale into the bins
    cutscore <- seq(from=min.crscore, to=max.crscore, length.out=num.bin + 1)
    
    # assign a group variable to each score
    bin_ref <- cut(crscore_ref, breaks=cutscore, include.lowest=TRUE, dig.lab=7)
    bin_foc <- cut(crscore_foc, breaks=cutscore, include.lowest=TRUE, dig.lab=7)
    
    # create a temporary data frame of bin frequency for both groups
    tmp_df <- 
      data.frame(resp=resp.ref, bin=bin_ref) %>%
      tidyr::drop_na() %>% 
      dplyr::group_by(.data$bin, .drop=FALSE) %>% 
      dplyr::summarise(n.ref=dplyr::n()) %>% 
      dplyr::left_join(y=data.frame(resp=resp.foc, bin=bin_foc) %>%
                         tidyr::drop_na() %>% 
                         dplyr::group_by(.data$bin, .drop=FALSE) %>% 
                         dplyr::summarise(n.foc=dplyr::n()), by="bin") %>% 
      dplyr::ungroup()
    
    # check if the counts of remaining sample is greater than equal to minimum a criterion 
    isok_ref <- (sum(tmp_df$n.ref[tmp_df$n.ref >= min.binsize & tmp_df$n.foc >= min.binsize]) / 
                   sum(tmp_df$n.ref)) >= 1 - max.del
    isok_foc <- (sum(tmp_df$n.foc[tmp_df$n.ref >= min.binsize & tmp_df$n.foc >= min.binsize]) / 
                   sum(tmp_df$n.foc)) >= 1 - max.del
    
    # if the criterion is met, then break out the loop
    if(all(isok_ref, isok_foc)) {
      break
    }
    
  }
  
  # final data frame containing all components to compute the beta statistic
  item_df <- 
    data.frame(resp=resp.ref, bin=bin_ref) %>%
    tidyr::drop_na() %>% 
    dplyr::group_by(.data$bin, .drop=FALSE) %>% 
    dplyr::summarise(n.ref=dplyr::n(), prop.ref=sum(.data$resp == 1, na.rm=TRUE)/.data$n.ref, 
                     var.ref=stats::var(.data$resp, na.rm=TRUE)) %>% 
    dplyr::left_join(y=data.frame(resp=resp.foc, bin=bin_foc) %>%
                       tidyr::drop_na() %>% 
                       dplyr::group_by(.data$bin, .drop=FALSE) %>% 
                       dplyr::summarise(n.foc=dplyr::n(), 
                                        prop.foc=sum(.data$resp == 1, na.rm=TRUE)/.data$n.foc, 
                                        var.foc=stats::var(.data$resp, na.rm=TRUE)), by="bin") %>% 
    dplyr::ungroup() %>%
    dplyr::filter(.data$n.ref >= 3 & .data$n.foc >= 3) %>% 
    dplyr::mutate(n.total=.data$n.ref + .data$n.foc, 
                  weight=dplyr::case_when(
                    weight.group == "foc" ~ .data$n.foc / sum(.data$n.foc),
                    weight.group == "ref" ~ .data$n.ref / sum(.data$n.ref),
                    weight.group == "comb" ~ .data$n.total / sum(.data$n.total))) %>%
    dplyr::mutate(
      beta=(.data$prop.ref - .data$prop.foc) * .data$weight, 
      var.beta=(.data$var.ref / .data$n.ref + .data$var.foc / .data$n.foc) * 
        .data$weight^2) 
  
  # compute the beta statistic and its SE 
  beta <- sum(item_df$beta)
  se_beta <- sqrt(sum(item_df$var.beta))
  z_beta <- beta / se_beta
  pval_beta <- 2 * stats::pnorm(q=abs(z_beta), mean=0, sd=1, lower.tail=FALSE)  
  n.ref=sum(item_df$n.ref)
  n.foc=sum(item_df$n.foc)
  stat_df <- data.frame(beta=beta, se=se_beta, z.beta=z_beta, p.val=pval_beta, 
                        n.ref=n.ref, n.foc=n.foc, n.total=n.ref + n.foc)
  
  # round the numbers of the data frame
  item_df2 <- 
    janitor::adorn_totals(dat=item_df, where="row", fill=NA,,, .data$n.ref, .data$n.foc, .data$n.total, 
                          .data$beta, .data$var.beta) %>% 
    dplyr::mutate_at(.vars=c(3, 4, 6, 7), "round", digits=4)
  
  # return the results
  list(dif_stat=stat_df, contingency=item_df2)
  
}

