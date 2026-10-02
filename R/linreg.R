#' Linear Regression Model
#'
#' Fits a linear regression model using linear algebra.
#'
#' @param formula A formula object (e.g., \code{Petal.Length ~ Species}).
#' @param data A data frame containing the variables.
#'
#' @return An object of class \code{linreg}.
#' @export
linreg <- function(formula, data){
  mf <- model.frame(formula,data)
  X <- model.matrix(formula,mf)
  Y <- data[[all.vars(formula)[1]]]

  beta<-solve(t(X)%*%X) %*% (t(X)%*%Y)
  fitted_values <- X %*% beta
  residuals <- Y - fitted_values

  n <- nrow(X)
  p <- ncol(X)
  df_residual <- n-p

  ssr <-sum(residuals^2)
  residual_variance <- ssr/df_residual

  vorc <- residual_variance * solve(t(X) %*% X)

  se_beta <- sqrt(diag(vorc))
  t_values <- as.vector(beta)/ se_beta

  p_values <- 2 * pt(-abs(t_values), df = df_residual)

  qr_decomp <- qr(X)
  Q <- qr.Q(qr_decomp)
  R <- qr.R(qr_decomp)
  R_inv <- backsolve(R, diag(ncol(X)))
  qr_beta_new <- R_inv %*% t(Q) %*% Y
  qr_vcov_matrix <- residual_variance * (R_inv %*% t(R_inv))

  out <- list(
    formula = formula,
    data_name = deparse(substitute(data)),
    coefficients = beta,
    fitted_values = fitted_values,
    residuals = residuals,
    df_residual = df_residual,
    residual_variance = residual_variance,
    sigma_sq = residual_variance,
    vcov = vorc,
    se = se_beta,
    t_values = t_values,
    p_values = p_values,
    qr_beta_new = qr_beta_new,
    qr_vcov_matrix = qr_vcov_matrix
  )

  class(out) <- "linreg"

  return(out)



}

#' @export
print.linreg <- function(x, ...) {
  cat("Call:\n")
  cat("linreg(formula = ", deparse(x$formula), ", data = ", x$data_name, ")\n\n", sep = "")

  cat("Coefficients:\n")
  coef_vector <- as.vector(x$coefficients)
  names(coef_vector) <- rownames(x$coefficients)

  print(coef_vector)
}

#' @export
#' @import ggplot2
plot.linreg <- function(x, ...) {
  fitted_vals <- as.vector(x$fitted_values)
  residuals_vals <- as.vector(x$residuals)

  std_residuals <- residuals_vals / sqrt(x$sigma_sq)
  sqrt_abs_std_res <- sqrt(abs(std_residuals))

  df_plot <- data.frame(
    Fitted = fitted_vals,
    Residuals = residuals_vals,
    Sqrt_Abs_Std_Res = sqrt_abs_std_res,
    Observation = seq_along(residuals_vals)
  )

  top_outliers <- order(abs(residuals_vals), decreasing = TRUE)[1:3]

  p1 <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Fitted, y = Residuals)) +
    ggplot2::geom_point(shape = 21, color = "#00b0f0", fill = "#00b0f0", alpha = 0.7) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_smooth(method = "loess", se = FALSE, color = "#00b5b8", formula = y ~ x) +
    ggplot2::geom_text(
      data = df_plot[top_outliers, ],
      ggplot2::aes(label = Observation),
      hjust = -0.3, vjust = 0.3, size = 3
    ) +
    ggplot2::labs(
      title = "Residuals vs Fitted",
      x = paste0("Fitted values\nlinreg(", deparse(x$formula), ")"),
      y = "Residuals"
    ) +
    ggplot2::theme_light()

  p2 <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Fitted, y = Sqrt_Abs_Std_Res)) +
    ggplot2::geom_point(shape = 21, color = "#00b0f0", fill = "#00b0f0", alpha = 0.7) +
    ggplot2::geom_smooth(method = "loess", se = FALSE, color = "#00e3a5", formula = y ~ x) +
    ggplot2::geom_text(
      data = df_plot[top_outliers, ],
      ggplot2::aes(label = Observation),
      hjust = -0.3, vjust = 0.3, size = 3
    ) +
    ggplot2::labs(
      title = "Scale-Location",
      x = paste0("Fitted values\nlinreg(", deparse(x$formula), ")"),
      y = expression(sqrt("|Standardized residuals|"))
    ) +
    ggplot2::theme_light()

  print(p1)
  print(p2)
}

#' @export
resid.linreg <- function(object, ...) {
  return(as.vector(object$residuals))
}

#' @export
pred <- function(object, ...) {
  UseMethod("pred")
}

#' @export
pred.linreg <- function(object, ...) {
  return(as.vector(object$fitted_values))
}

#' @export
coef.linreg <- function(object, ...) {
  coef_vec <- as.vector(object$coefficients)

  if (!is.null(rownames(object$coefficients))) {
    names(coef_vec) <- rownames(object$coefficients)
  }

  return(coef_vec)
}

#' @export
summary.linreg <- function(object, ...) {
  cat("Call:\n")
  cat("linreg(formula = ", deparse(object$formula), ", data = ", object$data_name, ")\n\n", sep = "")

  coef_names <- rownames(object$coefficients)

  summary_mat <- cbind(
    Estimate = as.vector(object$coefficients),
    `Std. Error` = object$se,
    `t value` = object$t_values,
    `Pr(>|t|)` = object$p_values
  )

  if (!is.null(coef_names)) {
    rownames(summary_mat) <- coef_names
  }

  cat("Coefficients:\n")
  stats::printCoefmat(summary_mat, P.values = TRUE, has.Pvalue = TRUE)

  cat("\nResidual standard error:", sqrt(object$residual_variance), "on", object$df_residual, "degrees of freedom\n")
}



