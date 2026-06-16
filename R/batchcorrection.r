#' Batch Correction with SVA
#'
#' This function will do batch correction on mvalues,
#'
#' @param mval matrix of mvalues
#' @param pdata sample Sheet (dataframe)
#' @param model model to apply with sva (string) eg: "~group+gender"
#'
#' @return A matrix of batch corrected mvalues
#'
#' @importFrom sva sva
#' @importFrom stats model.matrix residuals lm
bc_sva <- function(mval, pdata, model) {

  formula1 <- stats::as.formula(tolower(model))
  design <- stats::model.matrix(formula1, data = pdata)
  if (nrow(design) > ncol(mval)) {
    stop("Missing samples in betas ! ")
  }
  if (nrow(design) < ncol(mval)) {
    stop("Missing samples in pdata ! ")
  }

  if (!grepl("\\+", model)) {
    sva_m <- sva::sva(mval, design)
  } else {
    formula0 <- stats::as.formula(gsub("~[^+]*\\+", "~1+", model))
    model0 <- stats::model.matrix(formula0, data = pdata)
    sva_m <- sva::sva(mval, design, model0)
  }

  message(paste0("number of surrogate variables", sva_m$n.sv))
  if (!sva_m$n.sv) {
    message("0 surrogate variables have been found, ",
            "batchcorrection is useless !")
  } else {
    mval <- t(stats::residuals(stats::lm(t(mval) ~ sva_m$sv)))
  }

  mval
}
