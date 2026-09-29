#' Batch Correction with SVA
#'
#' This function will do batch correction on mvalues,
#'
#' @param mval matrix of mvalues
#' @param pdata sample Sheet (dataframe)
#' @param model model to apply with sva (string) eg: "~group+gender"
#' @param model0 null model to apply with sva (string) eg: "~gender"
#'
#' @return A matrix of batch corrected mvalues
#'
#' @importFrom sva sva
#' @importFrom stats model.matrix residuals lm
bc_sva <- function(mval, pdata, model, model0 = NULL) {

  formula1 <- stats::as.formula(tolower(model))
  design <- stats::model.matrix(formula1, data = pdata)

  if (nrow(design) > ncol(mval)) {
    stop("Missing samples in betas ! ")
  }
  if (nrow(design) < ncol(mval)) {
    stop("Missing samples in pdata ! ")
  }

  # Build matrix for null model (model0) based on user input or default behavior
  if (!is.null(model0)) {
    # model explicitly provided by user
    formula0 <- stats::as.formula(tolower(model0))
    design0 <- stats::model.matrix(formula0, data = pdata)
  } else if (!grepl("\\+", model)) {
    # Model with only one variable (ex: "~ grp") -> model0 = ~ 1
    formula0 <- stats::as.formula("~ 1")
    design0 <- stats::model.matrix(formula0, data = pdata)
  } else {
    # Default fallback for simple additive formula
    formula0 <- stats::as.formula(gsub("~[^+]*\\+", "~1+", model))
    design0 <- stats::model.matrix(formula0, data = pdata)
  }

  # SVA correction with the provided design and null model
  sva_m <- sva::sva(mval, design, design0)

  message(paste0("Number of surrogate variables found: ", sva_m$n.sv))

  if (sva_m$n.sv == 0) {
    message(
      "0 surrogate variables have been found, batchcorrection is useless !"
    )
  } else {
    mval <- t(stats::residuals(stats::lm(t(mval) ~ sva_m$sv)))
  }

  mval
}
