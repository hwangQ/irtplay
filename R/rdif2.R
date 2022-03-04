#' IRT residual-based DIF (RDIF) detection framework
#' 
#' @description This function computes three RDIF statistics (Lim, Choe, & Han, 2022; Lim, Choe, Han, Lee, & Hong, 2021), 
#' which are \eqn{RDIF_{R}}, \eqn{RDIF_{S}}, and \eqn{RDIF_{RS}}, for each item. \eqn{RDIF_{R}} primarily 
#' captures the typical contrast in raw residual pattern between two groups caused by uniform DIF whereas 
#' \eqn{RDIF_{S}} primarily captures the typical contrast in squared residual pattern between two groups caused
#' by nonuniform DIF. \eqn{RDIF_{RS}} can reasonably capture both types of DIF.   
#' 
#' @param x A data frame containing the item metadata (e.g., item parameters, number of categories, models ...), an object of class \code{\link{est_item}}
#' obtained from the function \code{\link{est_item}}, or an object of class \code{\link{est_irt}} obtained from the function \code{\link{est_irt}}.
#' The data frame of item metadata can be easily obtained using the function \code{\link{shape_df}}. See \code{\link{est_irt}}, \code{\link{irtfit}}, 
#' \code{\link{test.info}} or \code{\link{simdat}} for more details about the item metadata.
#' @param data A matrix containing examinees' response data for the items in the argument \code{x}. A row and column indicate
#' the examinees and items, respectively.
#' @param score A vector of examinees' ability estimates. If the abilities are not provided, \code{\link{rdif}} function estimates the abilities before 
#' computing RDIF statistics. See \code{\link{est_score}} for more details about scoring methods. Default is NULL. 
#' @param se A vector
#' @param reg.cr A logical value
#' @param group A numeric or character vector indicating group membership of examinees. The length of vector should the same with the number of rows 
#' in the response data matrix. 
#' @param focal.name A single numeric or character indicating the level of group which corresponds to the focal group. 
#' For example, if \code{group = c(0, 1, 0, 1, 1)} and '1' indicates the focal group, then \code{focal.name = 1}. 
#' @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
#' Default is 1.
#' @param alpha A numeric value to specify significance \eqn{\alpha}-level of the hypothesis test using the RDIF fit statistics.
#' Default is .05.
#' @param missing A value indicating missing values in the response data set. Default is NA.
#' @param purify A logical value indicating whether a purification process will be implemented or not. Default is FALSE.   
#' @param purify.by A character string specifying a RDIF statistic with which the purification is implemented. Available statistics 
#' are "rdif_rs" for \eqn{RDIF_{RS}}, "rdif_r" for \eqn{RDIF_{R}}, and "rdif_s" for \eqn{RDIF_{S}}. 
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
#' @details The RDIF framework (Lim et al., 2022; Lim et al., 2021) consists of three IRT residual-based statistics: \eqn{RDIF_{R}}, \eqn{RDIF_{S}}, 
#' and \eqn{RDIF_{RS}}. Under the null hypothesis that a test contains no DIF items, \eqn{RDIF_{R}} and \eqn{RDIF_{S}} follow 
#' normal distributions asymptotically. \eqn{RDIF_{RS}} is a based on a bivariate normal distribution of \eqn{RDIF_{R}} and 
#' \eqn{RDIF_{S}} statistics. Under the null hypothesis of no DIF items, it follows a \eqn{\chi^{2}} distribution asymptotically 
#' with 2 degrees of freedom. See Lim et al. (2022) for more details about RDIF framework. 
#' 
#' The \code{\link{rdif}} function computes all three RDIF statistics of \eqn{RDIF_{R}}, \eqn{RDIF_{S}}, and \eqn{RDIF_{RS}}. The current 
#' version of \code{\link{rdif}} function only supports dichotomous item response data. To compute the three statistics, the \code{\link{rdif}} function 
#' requires (1) item parameter estimates obtained from aggregate data regardless of group membership, (2) examinees' ability estimates 
#' (e.g., MLE), and (3) examinees' item response data. Note that the ability estimates need to be computed using the aggregate data-based 
#' item parameter estimates. The item parameter estimates should be provided in the \code{x} argument, the ability estimates should 
#' be provided in the \code{score} argument, and the response data should be provided in the \code{data} argument. When the abilities 
#' are not given in the \code{score} argument (i.e., \code{score = NULL}), the \code{\link{rdif}} function estimates examinees' abilities 
#' automatically using the scoring method specified in the \code{method} argument (e.g., \code{method = "MLE"}). 
#' 
#' The \code{group} argument accepts a vector of either two distinct numeric or character variables. Between two distinct variable, one is to 
#' represent the reference group and another one is to represent the focal group. The length of the vector should be the same with the number
#' of rows in the response data and each value in the vector should indicate each examinee of the response data. Once the \code{gruop} is 
#' specified, a single numeric or character value needs to be provided in the \code{focal.name} argument to define which group variable in 
#' the \code{group} argument represents the focal group.   
#' 
#' As other DIF detection approaches, an iterative purification process can be implemented for the RDIF framework. 
#' When \code{purify = TRUE}, the purification process is implemented based on one of RDIF statistics specified in the \code{purify.by} 
#' argument (e.g, \code{purify.by="rdif_rs"}). At each iterative purification, examinees' latent abilities are computed using purified items and 
#' scoring method specified in the \code{method} argument. The iterative purification process stops when no further DIF items are found or 
#' the process reaches a predetermined limit of iteration, which can be specified in the \code{max.iter} argument. See Lim et al. (2022) 
#' for more details about the purification procedure.
#' 
#' Scoring with a few items entails large standard errors which in turn could compromise DIF detection with RDIF framework. 
#' The \code{min.resp} argument can be used to avoid using scores with large standard errors when computing the RDIF statistics, espeically 
#' during the purification process. For example, if \code{min.resp} is not NULL (e.g., \code{min.resp=5}), item responses of examinees 
#' whose tally of item responses are less than the specified minimum number are treated as missing values (i.e., NA). Accordingly, 
#' their ability estimates become missing values and are not used for computing the RDIF statistics. If \code{min.resp=NULL}, 
#' an examinee's score will be computed as long as there exists, at least, 1 item response for the examinee. 
#' 
#' 
#' @return This function returns a list of four internal objects. The four objects are: 
#' \item{no_purify}{A list of several sub-objects containing the results of DIF analysis without a purification procedure. The sub-objects are: 
#'     \describe{
#'       \item{dif_stat}{A data frame containing the results of three RDIF statistics across all evaluated items. From the first column, each column 
#'        indicates item's ID, \eqn{RDIF_{R}} statistic, standardized \eqn{RDIF_{R}}, \eqn{RDIF_{S}} statistic, standardized, \eqn{RDIF_{S}}, 
#'        \eqn{RDIF_{RS}} statistic, p-value of the \eqn{RDIF_{R}}, p-value of the \eqn{RDIF_{S}}, p-value of the \eqn{RDIF_{RS}}, sample size of 
#'        the reference group, sample size of the focal group, and total sample size, respectively. Note that \eqn{RDIF_{RS}} does not have its standardized 
#'        value because it is a \eqn{\chi^{2}} statistic.} 
#'       \item{moments}{A data frame containing the moments of three RDIF statistics. From the first column, each column indicates item's ID, 
#'        mean of \eqn{RDIF_{R}}, standard deviation of \eqn{RDIF_{R}}, mean of \eqn{RDIF_{S}}, standard deviation of \eqn{RDIF_{S}}, and
#'        covariance of \eqn{RDIF_{R}} and \eqn{RDIF_{S}}, respectively.}
#'       \item{dif_item}{A list of three numeric vectors showing potential DIF items flagged by each of the RDIF statistics. Each of the numeric vector 
#'        means the items flagged by \eqn{RDIF_{R}}, \eqn{RDIF_{S}}, and \eqn{RDIF_{RS}}, respectively.}
#'       \item{score}{A vector of ability estimates used to compute the RDIF statistics.}
#'    }
#' }
#' \item{purify}{A logical value indicating whether the purification process was used.} 
#' \item{with_purify}{A list of several sub-objects containing the results of DIF analysis with a purification procedure. The sub-objects are:
#'     \describe{
#'       \item{purify.by}{A character string indicating which RDIF statistic is used for the purification. "rdif_r", "rdif_s", and "rdif_rs" refers to
#'        \eqn{RDIF_{R}}, \eqn{RDIF_{S}}, and \eqn{RDIF_{RS}}, respectively.}
#'       \item{dif_stat}{A data frame containing the results of three RDIF statistics across all evaluated items. From the first column, each column 
#'        indicates item's ID, \eqn{RDIF_{R}} statistic, standardized \eqn{RDIF_{R}}, \eqn{RDIF_{S}} statistic, standardized, \eqn{RDIF_{S}}, 
#'        \eqn{RDIF_{RS}} statistic, p-value of the \eqn{RDIF_{R}}, p-value of the \eqn{RDIF_{S}}, p-value of the \eqn{RDIF_{RS}}, sample size of 
#'        the reference group, sample size of the focal group, total sample size, and \emph{n}th iteration where the RDIF statistics were computed, 
#'        respectively.}
#'       \item{moments}{A data frame containing the moments of three RDIF statistics. From the first column, each column indicates item's ID, 
#'        mean of \eqn{RDIF_{R}}, standard deviation of \eqn{RDIF_{R}}, mean of \eqn{RDIF_{S}}, standard deviation of \eqn{RDIF_{S}}, covariance 
#'        of \eqn{RDIF_{R}} and \eqn{RDIF_{S}}, and \emph{n}th iteration where the RDIF statistics were computed, respectively.}
#'       \item{dif_item}{A list of three numeric vectors showing potential DIF items flagged by each of the RDIF statistics. Each of the numeric vector 
#'        means the items flagged by \eqn{RDIF_{R}}, \eqn{RDIF_{S}}, and \eqn{RDIF_{RS}}, respectively.}
#'       \item{n.iter}{A total number of iterations implemented for the purification.}
#'       \item{score}{A vector of final purified ability estimates used to compute the RDIF statistics.}
#'       \item{complete}{A logical value indicating whether the purification process was completed. If FALSE, it means that the purification process 
#'        reached the maximum iteration number but it was not complete.}
#'     }
#' }
#' \item{alpha}{A significance \eqn{\alpha}-level used to compute the p-values of RDIF statistics.}  
#' 
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso \code{\link{est_item}}, \code{\link{test.info}}, \code{\link{simdat}}, \code{\link{shape_df}}, 
#' \code{\link{gen.weight}}, \code{\link{est_score}}
#'
#' @references
#' Lim, H., Choe, E. M., & Han, K. T. (2022). A residual-based differential item functioning detection framework in 
#' item response theory. \emph{Journal of Educational Measurement}. 
#' 
#' Lim, H., Choe, E. M., Han, K. T., Lee, S., & Hong, M. (2021, June). \emph{IRT residual approach 
#' to detecting DIF.} Paper presented at the Annual Meeting of the National Council on Measurement 
#' in Education. Online. 
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
#' theta_ref <- rnorm(1000, 1.0, 1.0)
#' theta_foc <- rnorm(1000, 0.0, 1.0)
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
#' score_rst <- est_score(x=est_par, data=data, method="MLE")
#' score <- score_rst$est.theta
#' se <- score_rst$se.theta
#' 
#' ###############################################
#' # (3) conduct DIF analysis
#' ###############################################
#' # create a vector of group membership indicators 
#' # where '1' indicates the focal group 
#' group <- c(rep(0, 1000), rep(1, 1000))
#' 
#' # (a)-1 compute RDIF statistics by providing scores, 
#' #       and without a purification 
#' dif_nopuri_1 <- rdif(x=est_par, data=data, score=score, 
#'                      group=group, focal.name=1, D=1, alpha=0.05, purify=FALSE, purify.by="rdif_r")
#' print(dif_nopuri_1)
#'
#' dif_nopuri_1_cr <- rdif2(x=est_par, data=data, score=score, se=se, reg.cr=TRUE,
#'                          group=group, focal.name=1, D=1, alpha=0.05, purify=FALSE, purify.by="rdif_r")
#' print(dif_nopuri_1_cr)
#'   
#' # (a)-2 compute RDIF statistics by not providing scores 
#' #       and without a purification 
#' dif_nopuri_2 <- rdif(x=est_par, data=data, score=NULL, 
#'                      group=group, focal.name=1, D=1, alpha=0.05, 
#'                      method="MLE")
#' print(dif_nopuri_2)
#' 
#' # (b)-1 compute RDIF statistics with a purification 
#' #       based on \eqn{RDIF_{R}}
#' dif_puri_r <- rdif(x=est_par, data=data, score=score, 
#'                    group=group, focal.name=1, D=1, alpha=0.05, 
#'                    purify=TRUE, purify.by="rdif_r")
#' print(dif_puri_r)
#' 
#' # (b)-2 compute RDIF statistics with a purification 
#' #       based on \eqn{RDIF_{S}}
#' dif_puri_s <- rdif(x=est_par, data=data, score=score, 
#'                    group=group, focal.name=1, D=1, alpha=0.05, 
#'                    purify=TRUE, purify.by="rdif_s")
#' print(dif_puri_s)
#' 
#' # (b)-3 compute RDIF statistics with a purification 
#' #       based on \eqn{RDIF_{RS}}
#' dif_puri_rs <- rdif(x=est_par, data=data, score=score, 
#'                     group=group, focal.name=1, D=1, alpha=0.05, 
#'                     purify=TRUE, purify.by="rdif_rs")
#' print(dif_puri_rs)
#' }
#' 
#' @export
rdif2 <- function(x, data, score=NULL, se=NULL, reg.cr=FALSE, level.rel="group", group, 
                  focal.name, D=1, alpha=0.05, missing=NA, purify=FALSE, 
                  purify.by=c("rdif_rs", "rdif_r", "rdif_s"), max.iter=10, min.resp=NULL, method="MLE", 
                  range=c(-4, 4), norm.prior=c(0, 1), nquad=41, weights=NULL, ncore=1, verbose=TRUE) {
  
  
  # match.call
  cl <- match.call()
  
  ##----------------------------------
  ## (1) prepare DIF analysis
  ##----------------------------------
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
  
  # compute the score if score = NULL
  if(!is.null(score)) {
    # transform scores to a vector form
    if(is.matrix(score) | is.data.frame(score)) {
      score <- as.numeric(data.matrix(score))
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
                           nquad=nquad, weights=weights, ncore=ncore)$est.theta
    score <- score_rst$est.theta
    se <- score_rst$se.theta
  }
  
  # a) when no purification is set
  # do only one iteration of DIF analysis
  dif_rst <- rdif_one2(x=x, data=data, score=score, se=se, reg.cr=reg.cr, level.rel=level.rel, 
                       range=range, group=group, focal.name=focal.name, D=D, alpha=alpha)
  
  # create two empty lists to contain the results
  # no_purify <- list(dif_stat=NULL, effect_size=NULL, moments=NULL, dif_item=NULL)
  # with_purify <- list(purify.by=NULL, dif_stat=NULL, effect_size=NULL, moments=NULL, 
  #                     dif_item=NULL, n.iter=NULL, score=NULL, complete=NULL)
  no_purify <- list(dif_stat=NULL, moments=NULL, dif_item=NULL, score=NULL)
  with_purify <- list(purify.by=NULL, dif_stat=NULL, moments=NULL, 
                      dif_item=NULL, n.iter=NULL, score=NULL, complete=NULL)
  
  # record the first DIF detection results into the no purification list  
  no_purify$dif_stat <- dif_rst$dif_stat
  # no_purify$effect_size <- dif_rst$effect_size
  no_purify$dif_item <- dif_rst$dif_item
  no_purify$moments <- data.frame(id=x$id, 
                                  dif_rst$moments$rdif_r[, c(1, 3)],
                                  dif_rst$moments$rdif_s[, c(1, 3)], 
                                  dif_rst$covariance, stringsAsFactors=FALSE)
  names(no_purify$moments) <- c("id", "mu.rdif_r", "sigma.rdif_r", "mu.rdif_s", "sigma.rdif_s", "covariance")
  no_purify$score <- score
  
  # when purification is used
  if(purify) {
    
    # verify the criterion for purification
    purify.by <- match.arg(purify.by)
    
    # create an empty vector and empty data frames 
    # to contain the detected DIF items, statistics, and moments
    dif_item <- NULL
    dif_stat <- 
      data.frame(id=rep(NA_character_, nrow(x)), rdif_r=NA, z.rdif_r=NA,
                 rdif_s=NA, z.rdif_s=NA, rdif_rs=NA, p.val.rdif_r=NA, p.val.rdif_s=NA, p.val.rdif_rs=NA, 
                 n.ref=NA, n.foc=NA, n.total=NA, n.iter=NA, stringsAsFactors=FALSE)
    # efs_df <- 
    #   data.frame(id=rep(NA_character_, nrow(x)), he_rdif_r=NA, he_rdif_s=NA, 
    #              gl_rdif_r=NA, gl_rdif_s=NA, n.iter=NA, stringsAsFactors=FALSE)
    mmt_df <- 
      data.frame(id=rep(NA_character_, nrow(x)), mu.rdif_r=NA, sigma.rdif_r=NA, 
                 mu.rdif_s=NA, sigma.rdif_s=NA, covariance=NA, n.iter=NA, stringsAsFactors=FALSE)
    
    # extract the first DIF analysis results
    # and check if at least one DIF item is detected
    dif_item_tmp <- dif_rst$dif_item[[purify.by]]
    dif_stat_tmp <- dif_rst$dif_stat
    # efs_df_tmp <- dif_rst$effect_size
    mmt_df_tmp <- no_purify$moments
    
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
        flag_max <- 
          switch(purify.by,
                 rdif_r=which.max(abs(dif_stat_tmp$z.rdif_r)), 
                 rdif_s=which.max(abs(dif_stat_tmp$z.rdif_s)), 
                 rdif_rs=which.max(dif_stat_tmp$rdif_rs))
        
        # check an item that is deleted
        del_item <- item_num[flag_max]
        
        # add the deleted item as the DIF item
        dif_item <- c(dif_item, del_item)
        
        # add the DIF statistics and moments for the detected DIF item
        dif_stat[del_item, 1:12] <- dif_stat_tmp[flag_max, ]
        dif_stat[del_item, 13] <- i - 1
        # efs_df[del_item, 1:5] <- efs_df_tmp[flag_max, ]
        # efs_df[del_item, 6] <- i - 1
        mmt_df[del_item, 1:6] <- mmt_df_tmp[flag_max, ]
        mmt_df[del_item, 7] <- i - 1
        
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
        score_puri_rst <- est_score(x=x_puri, data=data_puri, D=D, method=method, range=range, norm.prior=norm.prior, 
                                    nquad=nquad, weights=weights, ncore=ncore)
        score_puri <- score_puri_rst$est.theta
        se_puri <- score_puri_rst$se.theta
        
        # do DIF analysis using the updated ability estimates
        dif_rst_tmp <- rdif_one2(x=x_puri, data=data_puri, score=score_puri, se=se_puri, reg.cr=reg.cr, 
                                 level.rel=level.rel, range=range, group=group, focal.name=focal.name, D=D, alpha=alpha)
        
        # extract the first DIF analysis results
        # and check if at least one DIF item is detected
        dif_item_tmp <- dif_rst_tmp$dif_item[[purify.by]]
        dif_stat_tmp <- dif_rst_tmp$dif_stat
        # efs_df_tmp <- dif_rst_tmp$effect_size
        mmt_df_tmp <- data.frame(id=dif_rst_tmp$dif_stat$id, 
                                 dif_rst_tmp$moments$rdif_r[, c(1, 3)],
                                 dif_rst_tmp$moments$rdif_s[, c(1, 3)], 
                                 dif_rst_tmp$covariance, stringsAsFactors=FALSE)
        names(mmt_df_tmp) <- c("id", "mu.rdif_r", "sigma.rdif_r", "mu.rdif_s", "sigma.rdif_s", "covariance")
        
        # check if a further DIF item is flagged
        if(is.null(dif_item_tmp)) {
          
          # add no additional DIF item
          dif_item <- dif_item
          
          # add the DIF statistics for rest of items
          dif_stat[item_num, 1:12] <- dif_stat_tmp
          dif_stat[item_num, 13] <- i
          # efs_df[item_num, 1:5] <- efs_df_tmp
          # efs_df[item_num, 6] <- i
          mmt_df[item_num, 1:6] <- mmt_df_tmp
          mmt_df[item_num, 7] <- i
          
          break
          
        } 
      }
      
      # print a message
      if(verbose) {
        cat("", "\n")
      }
      
      # record the actual number of iteration 
      n_iter <- i
      
      # if the iteration reached out the maximum number of iteration but the purification is incomplete, 
      # then, return a warning message
      if(max.iter == n_iter & !is.null(dif_item_tmp)) {
        warning("The iteration reached out the maximum number of iteration before purification is complete.", call.=FALSE)
        complete <- FALSE
        
        # add flagged DIF item at the last iteration
        dif_item <- c(dif_item, item_num[dif_item_tmp])
        
        # add the DIF statistics for rest of items
        dif_stat[item_num, 1:12] <- dif_stat_tmp
        dif_stat[item_num, 13] <- i
        # efs_df[item_num, 1:5] <- efs_df_tmp
        # efs_df[item_num, 6] <- i
        mmt_df[item_num, 1:6] <- mmt_df_tmp
        mmt_df[item_num, 7] <- i
        
        
      } else {
        complete <- TRUE
        
        # print a message
        if(verbose) {
          cat("Purification is finished.", '\n')
        }
        
      }
      
      # record the final DIF detection results with the purification procedure
      with_purify$purify.by <- purify.by
      with_purify$dif_stat <- dif_stat
      # with_purify$effect_size <- efs_df
      with_purify$moments <- mmt_df
      with_purify$dif_item <- sort(dif_item)
      with_purify$n.iter <- n_iter
      with_purify$score <- score_puri
      with_purify$se <- se_puri
      with_purify$complete <- complete
      
    } else {
      
      # in case when no DIF item is detected from the first DIF analysis results
      with_purify$purify.by <- purify.by
      with_purify$dif_stat <- cbind(no_purify$dif_stat, n.iter=0)
      # with_purify$effect_size <- cbind(no_purify$effect_size, n.iter=0)
      with_purify$moments <- cbind(no_purify$moments, n.iter=0)
      with_purify$n.iter <- 0
      with_purify$complete <- TRUE
      
    }
    
  }
  
  # summarize the results
  rst <- list(no_purify=no_purify, purify=purify, with_purify=with_purify, alpha=alpha) 
  
  # return the DIF detection results
  class(rst) <- "rdif"
  rst$call <- cl
  rst
  
  
}

# This function conducts one iteration of DIF analysis using the IRT residual based statistics
rdif_one2 <- function(x, data, score, se, reg.cr=FALSE, level.rel="group", 
                      range, group, focal.name, D=1, alpha=0.05) {
  
  # listrize the item metadata
  meta <- metalist2(x)
  
  # check the unique score categories
  cats <- x$cats
  
  # check the number of items
  nitem <- nrow(x)
  
  ##---------------------------------
  # compute the two statistics 
  ##---------------------------------
  # find the location of examinees for the reference and the focal groups 
  loc_ref <- which(group != focal.name)
  loc_foc <- which(group == focal.name)
  
  # divide the response data into the two group data 
  resp_ref <- data[loc_ref, , drop=FALSE]
  resp_foc <- data[loc_foc, , drop=FALSE]
  
  # check sample size
  n_ref <- colSums(!is.na(resp_ref))
  n_foc <- colSums(!is.na(resp_foc))
  
  # check if an item has all missing data for either of two groups
  all_miss <- sort(unique(c(which(n_ref == 0), which(n_foc == 0))))
  
  # divide the thetas into the two group data 
  score_ref <- score[loc_ref]
  score_foc <- score[loc_foc]
  
  # apply the regression correction for theta scores
  if(reg.cr) {
    
    # divide the SEs into the two group data 
    se_ref <- se[loc_ref]
    se_foc <- se[loc_foc]
    
    # compute the mean and variance of scores for each group
    # note that, the bound scores are excluded.
    lb_score <- range[1]
    up_score <- range[2]
    mu_ref <- mean(score_ref[score_ref > lb_score & score_ref < up_score], na.rm=TRUE)
    mu_foc <- mean(score_foc[score_foc > lb_score & score_foc < up_score], na.rm=TRUE)
    sigma2_ref <- stats::var(score_ref[score_ref > lb_score & score_ref < up_score], na.rm=TRUE)
    sigma2_foc <- stats::var(score_foc[score_foc > lb_score & score_foc < up_score], na.rm=TRUE)
    
    # compute the error variance of scores for each group
    level.rel <- tolower(level.rel)
    if(level.rel == "indiv") {
      
      # use individual level error variance      
      errvar_ref <- se_ref^2
      errvar_foc <- se_foc^2
      
      # compute the squared correlation (a.k.a. reliability) between theta estimate and true theta 
      rho_ref <- suppressWarnings(sqrt(1 - errvar_ref / sigma2_ref))
      rho_foc <- suppressWarnings(sqrt(1 - errvar_foc / sigma2_foc))
      rho_ref2 <- rho_ref^2
      rho_foc2 <- rho_foc^2
      rho_ref2[is.nan(rho_ref2)] <- NA_real_
      rho_foc2[is.nan(rho_foc2)] <- NA_real_
      
    } else if(level.rel == "group") {
      
      # use group level error variance      
      errvar_ref <- mean((se_ref^2)[score_ref > lb_score & score_ref < up_score], na.rm=TRUE)
      errvar_foc <- mean((se_foc^2)[score_foc > lb_score & score_foc < up_score], na.rm=TRUE)
      
      # compute the squared correlation (a.k.a. reliability) between theta estimate and true theta 
      rho_ref2 <- suppressWarnings(1 - errvar_ref / sigma2_ref)
      rho_foc2 <- suppressWarnings(1 - errvar_foc / sigma2_foc)      
    } 
    
    # apply a regression correction to the ability estimates
    score_ref <- mu_ref + rho_ref2 * (score_ref - mu_ref)
    score_foc <- mu_foc + rho_foc2 * (score_foc - mu_foc)
    
  }
  
  # compute the model-predicted probability of answering correctly (a.k.a. model-expected item score)
  extscore_ref <- trace(meta=meta, theta=score_ref, D=D)$trace_df
  extscore_foc <- trace(meta=meta, theta=score_foc, D=D)$trace_df
  
  # compute the model probability of score categories
  prob_ref <- trace4(meta=meta, theta=score_ref, D=D)$trace
  prob_foc <- trace4(meta=meta, theta=score_foc, D=D)$trace
  
  # replace NA values into the missing data location
  extscore_ref[is.na(resp_ref)] <- NA
  extscore_foc[is.na(resp_foc)] <- NA
  
  # compute the raw residuals   
  resid_ref <- resp_ref - extscore_ref
  resid_foc <- resp_foc - extscore_foc
  
  # compute the residual-based DIF statistic
  rdif_r <- colMeans(resid_foc, na.rm=TRUE) - colMeans(resid_ref, na.rm=TRUE) # the difference of mean raw residuals (DMRR)
  rdif_s <- colMeans(resid_foc^2, na.rm=TRUE) - colMeans(resid_ref^2, na.rm=TRUE) # the difference of mean squared residuals (DMSR)
  
  # compute the means and variances of the two statistics for the hypothesis testing
  moments <- resid_moments(p_ref=prob_ref, p_foc=prob_foc, n_ref=n_ref, n_foc=n_foc, 
                           resp_ref=resp_ref, resp_foc=resp_foc, cats=cats)
  moments_rdif_r <- moments$rdif_r
  moments_rdif_s <- moments$rdif_s
  covar <- moments$covariance
  poolsd <- moments$poolsd
  
  # compute the chi-square statistics
  chisq <- c()
  for(i in 1:nitem) {
    
    if(i %in% all_miss) {
      chisq[i] <- NaN
    } else{
      
      # create a var-covariance matrix between rdif_r and rdif_s
      cov_mat <- array(NA, c(2, 2))
      
      # replace NAs with the analytically computed covariance
      cov_mat[col(cov_mat) != row(cov_mat)] <- covar[i]
      
      # replace NAs with the analytically computed variances
      diag(cov_mat) <- c(moments_rdif_r[i, 2], moments_rdif_s[i, 2])
      
      # create a vector of mean rdif_r and mean rdif_s
      mu_vec <- cbind(moments_rdif_r[i, 1], moments_rdif_s[i, 1])
      
      # create a vector of rdif_r and rdif_s
      est_mu_vec <- cbind(rdif_r[i], rdif_s[i])
      
      # compute the chi-square statistic
      inv_cov <- suppressWarnings(tryCatch({solve(cov_mat, tol=1e-200)}, error = function(e) {NULL}))
      if(is.null(inv_cov)) {
        inv_cov <- suppressWarnings(tryCatch({solve(cov_mat + 1e-15, tol=1e-200)}, error = function(e) {NULL}))
        if(is.null(inv_cov)) {
          inv_cov <- suppressWarnings(tryCatch({solve(cov_mat + 1e-10, tol=1e-200)}, error = function(e) {NULL}))
        }
      }
      chisq[i] <- as.numeric((est_mu_vec - mu_vec) %*% inv_cov %*% t(est_mu_vec - mu_vec)) 
    }
    
  }
  
  # standardize the two statistics of rdif_r and rdif_s
  z_stat_rdif_r <- (rdif_r - moments_rdif_r$mu) / moments_rdif_r$sigma
  z_stat_rdif_s <- (rdif_s - moments_rdif_s$mu) / moments_rdif_s$sigma
  
  # calculate p-values for all three statistics
  pval_rdif_r <- round(2 * stats::pnorm(q=abs(z_stat_rdif_r), mean=0, sd=1, lower.tail=FALSE), 4)
  pval_rdif_s <- round(2 * stats::pnorm(q=abs(z_stat_rdif_s), mean=0, sd=1, lower.tail=FALSE), 4)
  pval_rdif_rs <- round(stats::pchisq(chisq, df=2, lower.tail=FALSE), 4)
  
  # compute three effect size for rdif_r and rdif_s
  rdif_stats <- list(x=rdif_r, y=rdif_s)
  n_total <- n_foc + n_ref
  
  # Cohen's d
  # efs_cohen <- 
  #   purrr::map2(.x=rdif_stats, .y=poolsd$cohen, .f=function(x, y) round(x / y, 4)) %>%
  #   data.frame() %>% 
  #   dplyr::rename_all(.f=function(x) c("co_rdif_r", "co_rdif_s"))
  
  # Hedge's g 
  efs_hedge <- 
    purrr::map2(.x=rdif_stats, .y=poolsd$hedge, 
                .f=function(x, y) {round((x / y) * (1 - (3 / (4 * (n_foc + n_ref) - 9))), 4)}) %>% 
    data.frame() %>% 
    dplyr::rename_all(.f=function(x) c("he_rdif_r", "he_rdif_s"))
  
  # Glass's delta
  efs_glass <- 
    purrr::map2(.x=rdif_stats, .y=poolsd$glass, .f=function(x, y) round(x / y, 4)) %>% 
    data.frame() %>% 
    dplyr::rename_all(.f=function(x) c("gl_rdif_r", "gl_rdif_s"))
  
  # combine all effect size
  effect_size <- data.frame(id=x$id, efs_hedge, efs_glass, stringsAsFactors=FALSE)
  
  # create a data frame to contain the results
  stat_df <- 
    data.frame(id=x$id, rdif_r=round(rdif_r, 4), z.rdif_r=round(z_stat_rdif_r, 4),
               rdif_s=round(rdif_s, 4), z.rdif_s=round(z_stat_rdif_s, 4), rdif_rs=round(chisq, 4), 
               p.val.rdif_r=pval_rdif_r, p.val.rdif_s=pval_rdif_s, p.val.rdif_rs=pval_rdif_rs, 
               n.ref=n_ref, n.foc=n_foc, n.total=n_total, stringsAsFactors=FALSE)
  rownames(stat_df) <- NULL
  
  # find the flagged items
  dif_item_rdif_r <- as.numeric(which(pval_rdif_r <= alpha))
  dif_item_rdif_s <- as.numeric(which(pval_rdif_s <= alpha))
  dif_item_rdif_rs <- which(pval_rdif_rs <= alpha)
  if(length(dif_item_rdif_r) == 0) dif_item_rdif_r <- NULL
  if(length(dif_item_rdif_s) == 0) dif_item_rdif_s <- NULL
  if(length(dif_item_rdif_rs) == 0) dif_item_rdif_rs <- NULL
  
  # summarize the results
  rst <- list(dif_stat=stat_df, effect_size=effect_size, 
              dif_item=list(rdif_r=dif_item_rdif_r, rdif_s=dif_item_rdif_s, rdif_rs=dif_item_rdif_rs), 
              moments=list(rdif_r=moments_rdif_r, rdif_s=moments_rdif_s), covariance=covar, alpha=alpha) 
  
  # return the results  
  rst
  
}

