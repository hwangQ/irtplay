#' @export
print.catsib <- function(x, digits = max(2L, getOption("digits") - 5L), ...) {
  
  call.expr <- deparse(x$call)
  cat("\nCall:\n", paste(call.expr, sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  
  # re-organize the dit stats data.frame
  dif_stat_nopurify <- 
    x$no_purify$dif_stat %>% 
    dplyr::select(.data$id, .data$n.ref, .data$n.foc, .data$n.total, 
                  .data$beta, .data$se, .data$z.beta, 
                  .data$p.val) %>% 
    dplyr::mutate_at(.vars = 5:8, round, digits = 3) %>% 
    dplyr::mutate(" "=stats::symnum(.data$p.val, c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                    symbols = c("***", "**", "*", ".", ""))) 
  
  # check if purification is used
  purify <- x$purify
  if(purify) {
    complete <- x$with_purify$complete
    n.iter <- x$with_purify$n.iter
    
    # re-organize the dit stats data.frame
    dif_stat_purify <- 
      x$with_purify$dif_stat %>% 
      dplyr::select(.data$id, .data$n.iter, .data$n.ref, .data$n.foc, .data$n.total, 
                    .data$beta, .data$se, .data$z.beta, .data$p.val) %>% 
      dplyr::mutate_at(.vars = 6:9, round, digits = 3) %>% 
      dplyr::mutate(" "=stats::symnum(.data$p.val, c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                      symbols = c("***", "**", "*", ".", "")))
    
  }
  
  ## print the results 
  cat("DIF analysis using CATSIB method", "\n\n")
  
  cat(" 1. Without purification \n\n")
  cat("  - Potential DIF Items: \n")
  cat("   ", paste(x$no_purify$dif_item, collapse=", "), "\n")
  cat("  - Test Statistic: \n\n")
  print(dif_stat_nopurify, digits=3, print.gap=NULL, quote=FALSE)
  cat("\n")
  cat("'***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1 ", 
      "\n")
  cat("Significance level:", x$alpha, "\n\n\n")
  
  cat(" 2. With purification \n\n")
  if(!purify) {
    cat("  - Purification was not implemented.", "\n\n")
  } else {
    cat("  - Completion of purification: ", complete,  "\n", sep="")
    cat("  - Number of iterations: ", n.iter,  "\n", sep="")
    cat("  - Potential DIF Items: \n", sep="")
    cat("   ", paste(x$with_purify$dif_item , collapse=", "), "\n")
    cat("  - Test Statistic: \n\n")
    print(dif_stat_purify, digits=3, print.gap=NULL, quote=FALSE)
    cat("\n")
    cat("'***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1 ", 
        "\n")
    cat("Significance level:", x$alpha, "\n\n")
  }
  
  invisible(x)
  
}

#' @export
print.summary.est_mg <- function(x, digits = max(2L, getOption("digits") - 5L), ...) {
  
  cat("\nCall:\n", paste(x$call.expr, sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Summary of the Data \n")
  cat(" Number of Items: \n", sep="")
  cat("  Overall: ", x$nitem$overall," unique items \n", sep="")
  cat("  By group: ", paste(paste(x$nitem$group, "(", x$group.name, ")", sep=""), collapse = ", "), "\n", sep="")
  cat(" Number of Cases: \n", sep="")
  cat("  Overall: ", x$ncase$overall,"\n", sep="")
  cat("  By group: ", paste(paste(x$ncase$group, "(", x$group.name, ")", sep=""), collapse = ", "), "\n\n", sep="")
  
  cat("Summary of Estimation Process \n")
  cat(" Maximum number of EM cycles: ", x$MaxE, "\n", sep="")
  cat(" Convergence criterion of E-step: ", x$Etol, "\n", sep="")
  cat(" Number of rectangular quadrature points: ", nrow(x$weights[[1]]), "\n", sep="")
  cat(" Minimum & Maximum quadrature points: ", x$weights[[1]][1, 1], ", ", -x$weights[[1]][1, 1], "\n", sep="")
  cat(" Number of free parameters: ", x$npar.est, "\n", sep="")
  cat(" Number of fixed items: \n", sep="")
  cat("  Overall: ", length(x$fix.loc$overall),"\n", sep="")
  cat("  By group: ", paste(paste(purrr::map_dbl(.x=x$fix.loc$group, length), "(", x$group.name, ")", sep=""), collapse = ", "), "\n", sep="")
  cat(" Number of E-step cycles completed: ", x$niter, "\n", sep="")
  cat(" Maximum parameter change: ", x$maxpar.diff, "\n\n", sep="")
  
  cat("Processing time (in seconds) \n")
  cat(" EM algorithm: ", x$EMtime, "\n", sep="")
  cat(" Standard error computation: ", x$SEtime, "\n", sep="")
  cat(" Total computation: ", x$TotalTime, "\n\n", sep="")
  
  cat("Convergence and Stability of Solution \n")
  cat(" First-order test: ", x$test.1, "\n", sep="")
  cat(" Second-order test: ", x$test.2, "\n", sep="")
  cat(" Computation of variance-covariance matrix: \n", "  ", x$var.note, "\n\n", sep="")
  
  cat("Summary of Estimation Results \n")
  cat(" -2loglikelihood: \n", sep="")
  cat("  Overall: ", round(-2 * x$loglikelihood$overall, 3),"\n", sep="")
  cat("  By group: ", paste(paste(round(-2 * unlist(x$loglikelihood$group), 3), "(", x$group.name, ")", sep=""), collapse = ", "), "\n\n", sep="")
  cat(" Akaike Information Criterion (AIC): ", round(x$aic, 3), "\n", sep="")
  cat(" Bayesian Information Criterion (BIC): ", round(x$bic, 3), "\n", sep="")
  cat(" Item Parameters (Overall): \n")
  item.par <- purrr::modify_if(.x=x$estimates$overall, .p=is.numeric, .f=round, digits=digits)
  print(item.par, print.gap=2, quote=FALSE)
  cat(" Group Parameters: \n")
  group.par <- 
    purrr::modify(.x=x$group.par, .f=round, digits=digits) %>% 
    dplyr::bind_rows()
  rownames(group.par) <- paste(rep(c("estimate", "se"), x$ngroup), "(", rep(x$group.name, each=2), ")", sep="")
  print(group.par, print.gap=2, quote=FALSE)
  cat("\n")
  invisible(x)
  
}


#' @export
print.est_mg <- function(x, digits = max(2L, getOption("digits") - 5L), ...) {
  
  call.expr <- deparse(x$call)
  cat("\nCall:\n", paste(call.expr, sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  
  cat("Multiple-Group Item parameter estimation using MMLE-EM. \n")
  cat(x$niter, " E-step cycles were completed using ", nrow(x$weights[[1]]), " quadrature points.", "\n", sep="")
  cat("First-order test: ", x$test.1, "\n", sep="")
  cat("Second-order test: ", x$test.2, "\n", sep="")
  cat("Computation of variance-covariance matrix: \n", "  ", x$var.note, "\n\n", sep="")
  cat("Log-likelihood: ", (x$loglikelihood$overall), "\n", sep="")
  
  cat("\n")
  invisible(x)
  
}

#' @export
print.rdif <- function(x, digits = max(2L, getOption("digits") - 5L), ...) {
  
  call.expr <- deparse(x$call)
  cat("\nCall:\n", paste(call.expr, sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  
  # re-organize the dit stats data.frame
  dif_stat_nopurify <- 
    x$no_purify$dif_stat %>% 
    dplyr::select(.data$id, .data$n.ref, .data$n.foc, .data$rdif_r, 
                  .data$p.val.rdif_r, .data$rdif_s, .data$p.val.rdif_s, 
                  .data$rdif_rs, .data$p.val.rdif_rs) %>% 
    dplyr::mutate_at(.vars = 4:9, round, digits = 3) %>% 
    dplyr::mutate(" "=stats::symnum(.data$p.val.rdif_r, c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                    symbols = c("***", "**", "*", ".", "")), 
                  "  "=stats::symnum(.data$p.val.rdif_s, c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                     symbols = c("***", "**", "*", ".", "")), 
                  "   "=stats::symnum(.data$p.val.rdif_rs, c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                      symbols = c("***", "**", "*", ".", ""))) %>% 
    dplyr::relocate(" ", .after=.data$p.val.rdif_r) %>% 
    dplyr::relocate("  ", .after=.data$p.val.rdif_s) %>% 
    dplyr::relocate("   ", .after=.data$p.val.rdif_rs) 
  
  # check if purification is used
  purify <- x$purify
  if(purify) {
    purify.by <- x$with_purify$purify.by
    complete <- x$with_purify$complete
    n.iter <- x$with_purify$n.iter
    if(purify.by == "rdif_r") {
      purify.stat <- "RDIF(R)"
    }
    if(purify.by == "rdif_s") {
      purify.stat <- "RDIF(S)"
    }
    if(purify.by == "rdif_rs") {
      purify.stat <- "RDIF(RS)"
    }
    
    # re-organize the dit stats data.frame
    dif_stat_purify <- 
      x$with_purify$dif_stat %>% 
      dplyr::select(.data$id, .data$n.iter, .data$n.ref, .data$n.foc, 
                    .data$rdif_r, .data$p.val.rdif_r, .data$rdif_s, .data$p.val.rdif_s, 
                    .data$rdif_rs, .data$p.val.rdif_rs) %>% 
      dplyr::mutate_at(.vars = 5:10, round, digits = 3) %>% 
      dplyr::mutate(" "=stats::symnum(.data$p.val.rdif_r, c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                      symbols = c("***", "**", "*", ".", "")), 
                    "  "=stats::symnum(.data$p.val.rdif_s, c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                       symbols = c("***", "**", "*", ".", "")), 
                    "   "=stats::symnum(.data$p.val.rdif_rs, c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                        symbols = c("***", "**", "*", ".", ""))) %>% 
      dplyr::relocate(" ", .after=.data$p.val.rdif_r) %>% 
      dplyr::relocate("  ", .after=.data$p.val.rdif_s) %>% 
      dplyr::relocate("   ", .after=.data$p.val.rdif_rs) 
    
  }
  
  ## print the results 
  cat("DIF analysis using three RDIF statistics", "\n\n")
  
  cat(" 1. Without purification \n\n")
  cat("  - DIF Items identified by RDIF(R): \n")
  cat("   ", paste(x$no_purify$dif_item$rdif_r, collapse=", "), "\n")
  cat("  - DIF Items identified by RDIF(S): \n")
  cat("   ", paste(x$no_purify$dif_item$rdif_s, collapse=", "), "\n")
  cat("  - DIF Items identified by RDIF(RS): \n")
  cat("   ", paste(x$no_purify$dif_item$rdif_rs, collapse=", "), "\n")
  cat("  - RDIF Statistics: \n\n")
  print(dif_stat_nopurify, digits=3, print.gap=NULL, quote=FALSE)
  cat("\n")
  cat("'***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1 ", 
      "\n")
  cat("Significance level:", x$alpha, "\n\n\n")
  
  cat(" 2. With purification \n\n")
  if(!purify) {
    cat("  - Purification was not implemented.", "\n\n")
  } else {
    cat("  - Completion of purification: ", complete,  "\n", sep="")
    cat("  - Number of iterations: ", n.iter,  "\n", sep="")
    cat("  - RDIF statistic used for purification: ", purify.stat, "\n", sep="")
    cat("  - DIF Items identified by ", purify.stat, ": \n", sep="")
    cat("   ", paste(x$with_purify$dif_item , collapse=", "), "\n")
    cat("  - RDIF Statistics: \n\n")
    print(dif_stat_purify, digits=3, print.gap=NULL, quote=FALSE)
    cat("\n")
    cat("'***'p < 0.001 '**'p < 0.01 '*'p < 0.05 '.'p < 0.1 ' 'p < 1 ", 
        "\n")
    cat("Significance level:", x$alpha, "\n\n")
  }
  
  invisible(x)
  
}


#' @export
print.irtfit <- function(x, ...) {

  call.expr <- deparse(x$call)
  cat("\nCall:\n", paste(call.expr, sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Significance level for chi-square fit statistic:", x$ancillary$alpha, "\n\n")
  cat("Item fit statistics: \n")
  print(x$fit_stat, print.gap=2, quote=FALSE)
  cat("\n")
  cat("Caution is needed in interpreting infit and outfit statistics for non-Rasch models. \n")
  invisible(x)

}


#' @export
print.est_irt <- function(x, digits = max(2L, getOption("digits") - 5L), ...) {

  call.expr <- deparse(x$call)
  cat("\nCall:\n", paste(call.expr, sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  cat("Item parameter estimation using MMLE-EM. \n")
  cat(x$niter, " E-step cycles were completed using ", nrow(x$weights), " quadrature points.", "\n", sep="")
  cat("First-order test: ", x$test.1, "\n", sep="")
  cat("Second-order test: ", x$test.2, "\n", sep="")
  cat("Computation of variance-covariance matrix: \n", "  ", x$var.note, "\n\n", sep="")
  cat("Log-likelihood: ", (x$loglikelihood), "\n", sep="")

  cat("\n")
  invisible(x)

}


#' @export
print.summary.est_irt <- function(x, digits = max(2L, getOption("digits") - 5L), ...) {

  cat("\nCall:\n", paste(x$call.expr, sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Summary of the Data \n")
  cat(" Number of Items: ", x$nitem, "\n", sep="")
  cat(" Number of Cases: ", x$ncase, "\n\n", sep="")

  cat("Summary of Estimation Process \n")
  cat(" Maximum number of EM cycles: ", x$MaxE, "\n", sep="")
  cat(" Convergence criterion of E-step: ", x$Etol, "\n", sep="")
  cat(" Number of rectangular quadrature points: ", nrow(x$weights), "\n", sep="")
  cat(" Minimum & Maximum quadrature points: ", x$weights[1, 1], ", ", -x$weights[1, 1], "\n", sep="")
  cat(" Number of free parameters: ", x$npar.est, "\n", sep="")
  cat(" Number of fixed items: ", length(x$fix.loc), "\n", sep="")
  cat(" Number of E-step cycles completed: ", x$niter, "\n", sep="")
  cat(" Maximum parameter change: ", x$maxpar.diff, "\n\n", sep="")

  cat("Processing time (in seconds) \n")
  cat(" EM algorithm: ", x$EMtime, "\n", sep="")
  cat(" Standard error computation: ", x$SEtime, "\n", sep="")
  cat(" Total computation: ", x$TotalTime, "\n\n", sep="")

  cat("Convergence and Stability of Solution \n")
  cat(" First-order test: ", x$test.1, "\n", sep="")
  cat(" Second-order test: ", x$test.2, "\n", sep="")
  cat(" Computation of variance-covariance matrix: \n", "  ", x$var.note, "\n\n", sep="")

  cat("Summary of Estimation Results \n")
  cat(" -2loglikelihood: ", round((-2 * x$loglikelihood), 3), "\n", sep="")
  cat(" Akaike Information Criterion (AIC): ", round(x$aic, 3), "\n", sep="")
  cat(" Bayesian Information Criterion (BIC): ", round(x$bic, 3), "\n", sep="")
  cat(" Item Parameters: \n")
  item.par <- purrr::modify_if(.x=x$estimates, .p=is.numeric, .f=round, digits=digits)
  print(item.par, print.gap=2, quote=FALSE)
  cat(" Group Parameters: \n")
  group.par <- round(x$group.par, digits=digits)
  print(group.par, print.gap=2, quote=FALSE)
  cat("\n")
  invisible(x)

}


#' @export
print.est_item <- function(x, digits = max(2L, getOption("digits") - 5L), ...) {

  call.expr <- deparse(x$call)
  cat("\nCall:\n", paste(call.expr, sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  cat("Fixed ability parameter calibration (Stocking's Method A). \n")
  cat(x$convergence, "\n\n")
  cat("Log-likelihood: ", (x$loglikelihood), "\n", sep="")

  cat("\n")
  invisible(x)

}


#' @export
print.summary.est_item <- function(x, digits = max(2L, getOption("digits") - 5L), ...) {

  cat("\nCall:\n", paste(x$call.expr, sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Summary of the Data \n")
  cat(" Number of Items in Response Data: ", x$nitem, "\n", sep="")
  cat(" Number of Excluded Items: ", x$nitem.del, "\n", sep="")
  cat(" Number of free parameters: ", x$npar.est, "\n", sep="")
  cat(" Number of Responses for Each Item: \n")
  print(data.frame(id=x$estimates$id, n=x$n.response), print.gap=2, quote=FALSE)
  cat("\n")

  cat("Processing time (in seconds) \n")
  cat(" Total computation: ", x$TotalTime, "\n\n", sep="")

  cat("Convergence of Solution \n")
  cat(" ", x$convergence, "\n\n", sep="")

  cat("Summary of Estimation Results \n")
  cat(" -2loglikelihood: ", round((-2 * x$loglikelihood), 3), "\n", sep="")
  cat(" Item Parameters: \n")
  item.par <- purrr::modify_if(.x=x$estimates, .p=is.numeric, .f=round, digits=digits)
  print(item.par, print.gap=2, quote=FALSE)
  cat("\n")
  cat(" Group Parameters: \n")
  group.par <- round(x$group.par, digits=digits)
  print(group.par, print.gap=2, quote=FALSE)
  cat("\n")
  invisible(x)

}
