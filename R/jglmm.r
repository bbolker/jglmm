utils::globalVariables(c("."))

#' @importFrom dplyr "%>%" mutate
#' @importFrom glue glue
#' @importFrom generics augment tidy
#' @importFrom JuliaCall julia_assign julia_command julia_eval
#' @importFrom rlang .data
#' @importFrom stats qnorm
NULL

julia_libraries <- c("MixedModels", "DataFrames", "StatsModels")

#' Set up Julia and required libraries
#'
#' @param run_glmm_model (logical): run simple GLMM so JIT compilation gets done (speed up timing of subsequent models)
#' @param quietly print progress messages?
#' @export
jglmm_setup <- function(run_glmm_model=FALSE, quietly=FALSE) {
  ## FIXME: add more messages, plus 'quiet
  JuliaCall::julia_setup(verbose=!quietly)
  for (lib in julia_libraries) {
      JuliaCall::julia_library(lib)
  }
  if (run_glmm_model) {
      if (!requireNamespace("lme4")) {
          warning("lme4 not available, skipping glmm shakedown run")
      } else {
          message("running glmm shakedown model ...")
          ## avoid NSE ...
          cbpp <- lme4::cbpp
          cbpp$prop <- cbpp$incidence / cbpp$size
          j <- jglmm(prop ~ period + (1 | herd), data = cbpp, family = "binomial",
                weights = cbpp$size)
      }
  } ## run_glmm_model
    invisible(NULL)
}


#' Fitting Generalized Linear Mixed-Effects Models in Julia
#'
#' @param formula A two-sided linear formula object describing both the
#'   fixed-effects and random-effects part of the model, with the response on
#'   the left of a ~ operator and the terms, separated by + operators, on the
#'   right. Random-effects terms are distinguished by vertical bars ("|")
#'   separating expressions for design matrices from grouping factors.
#' @param data A data frame containing the variables named in formula.
#' @param family (optional) The distribution family for the response variable
#'   (defaults to "normal").
#' @param link (optional) The model link function (defaults to "identity").
#' @param weights (optional) A vector of prior case weights.
#' @param contrasts (optional) A named list mapping column names of categorical
#'   variables in data to coding schemes (defauls to dummy coding all
#'   categorical variables).
#' @param return_val return fitted model ("jglmm") or just Julia model string ("julia_model_str")?
#'
#' @return An object of class `jglmm`.
#' @export
#'
#' @examples
#' \dontrun{
#' # linear model
#' lm1 <- jglmm(Reaction ~ Days + (Days | Subject), lme4::sleepstudy)
#'
#' # logistic model
#' cbpp <- dplyr::mutate(lme4::cbpp, prop = incidence / size)
#' gm1 <- jglmm(prop ~ period + (1 | herd), data = cbpp, family = "binomial",
#'              weights = cbpp$size)
#' gm2 <- jglmm(prop ~ period + (1 | herd), data = cbpp, family = "binomial",
#'              weights = cbpp$size, contrasts = list(period = "effects"))
#' }
jglmm <- function(formula, data, family = "normal", link = NULL, weights = NULL,
                  contrasts = NULL,
                  return_val=c("jglmm","julia_model_str")) {

  return_val <- match.arg(return_val)
  stopifnot(
    family %in% c("bernoulli", "binomial", "gamma", "normal", "poisson"),
    link %in% c("cauchit", "cloglog", "identity", "inverse", "logit", "log",
                "probit", "sqrt"),
    contrasts %in% c("dummy", "effects", "helmert")
  )

  # pass formula and data
  julia_assign("formula", formula)
  julia_assign("data", data)

  # construct model arguments
  model_args <- c("formula", "data")

  Family <- stringr::str_to_title(family)
  Link <- stringr::str_to_title(link)
  # choose between LinearMixedModel and GeneralizedLinearMixedModel
  if (family == "normal" && (is.null(link) || link == "identity")) {
    model_fun <- "MixedModels.LinearMixedModel"
  } else {
    model_fun <- "MixedModels.GeneralizedLinearMixedModel"
    model_args <- c(model_args, glue("{Family}()"))
    if (!is.null(link)) {
      model_args <- c(model_args, glue("{Link}Link()"))
    }
  }

  if (!is.null(contrasts)) {
    contrasts_args <- contrasts %>%
      purrr::map2_chr(names(.),
                      ~glue(":{.y} => {stringr::str_to_title(.x)}Coding()")) %>%
      paste(collapse = ", ")
    model_args <- c(model_args, glue("contrasts = Dict({contrasts_args})"))
  }

  if (!is.null(weights)) {
    julia_assign("weights", weights)
    model_args <- c(model_args, "wts = weights")
  }

  ## set up and fit model
  julia_model_str <- glue("fit({model_fun}, {paste(model_args, collapse = ', ')})")
  if (return_val=="julia_model_str") return(julia_model_str)
  model <- julia_eval(julia_model_str)

  results <- list(formula = formula, data = data, model = model)
  attr(results,"julia_model") <- julia_model_str
  class(results) <- "jglmm"
  return(results)

}

#' Tidying methods for jglmm models
#'
#' @param x An object of class `jglmm`, as returned by `jglmm`.
#'
#' @name jglmm_tidiers
#'
#' @examples
#' \dontrun{
#' cbpp <- dplyr::mutate(lme4::cbpp, prop = incidence / size)
#' gm <- jglmm(prop ~ period + (1 | herd), data = cbpp, family = "binomial",
#'             weights = cbpp$size)
#' tidy(gm)
#' augment(gm)
#' }
NULL

#' @rdname jglmm_tidiers
#'
#' @return `tidy` returns a tibble of fixed effect estimates
#' @param conf.int include confidence interval?
#' @param conf.level confidence level
#' @param ... extra arguments (ignored)
#' @export
tidy.jglmm <- function(x, conf.int=FALSE, conf.level=0.95, ...) {
  ## utils::globalVariables(c("estimate", "std.error")) ## induces dependence on R > 2.15.1 ..
  ## and isn't working as expected to suppress NOTEs?  
  estimate <- std.error <- NULL  
  julia_assign("model", x$model)
  julia_command("coef = coeftable(model);")
  julia_command("coef_df = DataFrame(coef.cols);")
  julia_command("rename!(coef_df, coef.colnms, makeunique = true);")
  julia_command("coef_df[!, :term] = coef.rownms;")
  r <- julia_eval("coef_df") %>%
      dplyr::as_tibble() %>%
      dplyr::select(.data$term, estimate = .data$Estimate,
                  std.error = .data$Std.Error, z.value = .data$`z value`,
                  p.value = .data$`P(>|z|)`)
  if (conf.int) {
      qn <- qnorm((1+conf.level)/2)
      r <- r %>%
          dplyr::mutate(conf.low=estimate-qn*std.error, conf.high=estimate+qn*std.error)
  }
  return(r)
}

#' @rdname jglmm_tidiers
#'
#' @return `augment` returns a tibble of the original data used to fit the model
#'   with an additional `.fitted` column containing the fitted response valuese.
#'
#' @export
augment.jglmm <- function(x) {
  julia_assign("model", x$model)
  fits <- julia_eval("fitted(model)")
  x$data$.fitted <- fits
  x$data %>% dplyr::as_tibble()
}
