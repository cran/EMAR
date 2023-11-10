
#' Fit statistics for model assessment
#'
#' This function helps users to generate 16 fit indices for model assessment. Once a model(s) is fitted, users can apply the function [fitstat()] to get the indices.
#' @param obj fitted model of the class lm, nls, gls, gnls, lme, or nlme.
#' @returns Model.name: fitted model name, p: number of parameters in the fitted model, n: number of observation, SSR: residual sum of squares, TRE: total relative error, Bias: mean bias, MRB: mean relative bias, MAB: mean absolute bias, MAPE: mean absolute percentage error, MSE: mean squared error, RMSE: root mean square error, Percent.RMSE: percentage root mean squared error, R2: coefficient of determination, R2adj: adjusted coefficient of determination, APC: Amemiya's prediction criterion, logL: Log-likelihood, AIC: Akaike information criterion, AICc: corrected Akaike information criterion, BIC: Bayesian information criterion, HQC: Hannan-Quin information criterion.
#' @note The lower the better for the SSR, TRE, Bias, MRB, MAB, MAPE, MSE, RMSE, Percent.RMSE, APC, AIC, AICc, BIC and HQC indices. The higher the better for R2 and R2adj indices. Users can choose which indices to use to evaluate their models from the output.
#' @keywords fitstat
#' @author Ogana F.N. and Corral-Rivas S.
#' @seealso [valstat()], which gives the fit indices of the model based on the independent/validation data
#' @importFrom stats coef
#' @importFrom stats fitted
#' @importFrom stats residuals
#' @importFrom stats logLik
#' @importFrom stats AIC
#' @importFrom stats BIC
#' @export
#' @examples
#' library(EMAR)
#'
#' # sample data
#' Age <- 1:50
#' Yield <- exp(6.5 - 39.5/Age)
#' fit_data <- data.frame(Age, Yield)
#'
#' # fit your model(s)
#' Eq01 <- lm(Yield ~ Age, data=fit_data)
#' Eq02 <- nls(Yield ~ b0 * Age ^ b1, data=fit_data, start=list(b0 = 2, b1 = 1))
#'
#' # Get the fit statistics for the model(s)
#' fitstat(Eq01)
#' fitstat(Eq02)
#'
#' # with the 'rbind' function, Users can generate output for multiple models at once.
#' indices <- rbind(fitstat(Eq01), fitstat(Eq02))
#' print(indices)

fitstat <- function(obj){
  Model.name <- deparse(substitute(obj))
  n <- length(residuals(obj))
  p <- length(coef(obj))
  y <- residuals(obj) + fitted(obj)
  SSR <- sum(residuals(obj)^2)
  SST <- sum((y - mean(y))^2)
  TRE <- sum(abs(residuals(obj))) / sum(fitted(obj))
  Bias <- mean(residuals(obj))
  MRB <- mean((fitted(obj)-y)/y)
  AbsBias <- mean(abs(residuals(obj)))
  MAPE <- (mean(abs(residuals(obj))/y))*100
  MSE <- SSR / n
  RMSE <- sqrt(SSR / n)
  Percent.RMSE <- (RMSE / mean(y)) * 100
  R2 <- 1 - SSR / SST
  R2adj <- 1 - SSR / SST * (n - 1) / (n - p)
  APC <- ((n + p) / (n - p))*(1 - (R2))
  LL <- logLik(obj)
  logL <- LL[[1]]
  AIC <- AIC(obj)
  AICc <- AIC + (2*p*(p+1))/(n-p-1)
  BIC <- BIC(obj)
  HQC <- 2 * p * log(log(n)) - 2 * logL

  return(data.frame(Model.name = Model.name, p = p, n = n, SSR = SSR, TRE = TRE, Bias = Bias, MRB = MRB, MAB = AbsBias, MAPE = MAPE,
                    MSE = MSE, RMSE = RMSE, Percent.RMSE = Percent.RMSE,
                    R2 = R2, R2adj = R2adj, APC = APC, logL = logL, AIC = AIC, AICc = AICc, BIC = BIC, HQC = HQC))
}

#' Validation statistics for model assessment
#'
#' This function helps users to generate 16 fit indices for model assessment based on independent/validation data set. In addition to empirical models, the function [valstat()] can generate fit indices for AI-based models such as artificial neural network, supervise vector machine, etc.
#' @param obs.y observed values from the independent/validation data
#' @param pred.y predicted values from the model
#' @param p number of parameters in the model. This is needed to compute the 'criteria-based indices' and adjusted coefficient of determination. Users could enter any value for AI-based models with an unknown number of parameters (p) and assess their models using the indices that are invariant of p. See the section on note.
#' @returns n: number of observation in the validation data, SSR: residual sum of squares, TRE: total relative error, Bias: mean bias, MRB: mean relative bias, MAB: mean absolute bias, MAPE: mean absolute percentage error, MSE: mean squared error, RMSE: root mean squared error, Percent.RMSE: percentage root mean squared error, R2: coefficient of determination, R2adj: adjusted coefficient of determination, APC: Amemiya's prediction criterion, logL: Log-likelihood, AIC: Akaike information criterion, AICc: corrected Akaike information criterion, BIC: Bayesian information criterion, HQC: Hannan-Quin information criterion.
#' @note The lower the better for the SSR, TRE, Bias, MRB, MAB, MAPE, MSE, RMSE, Percent.RMSE, APC, logL, AIC, AICc, BIC and HQC indices. The higher the better for R2 and R2adj indices. Users can choose which indices to use to evaluate their models from the output. Given the difficulty of determining the number of parameters (p) in AI-based models, users might consider using error-based indices, and coefficients of determination (R2).
#' @keywords valstat
#' @author Ogana F.N. and Corral-Rivas S.
#' @seealso [fitstat()], which gives the fit indices of the model based on the fitting data
#' @export
#' @examples
#' library(EMAR)
#'
#' # fitting data
#' Age <- 1:50
#' Yield <- exp(6.5 - 39.5/Age)
#' dat <- data.frame(Age, Yield)
#'
#' # fit the model to the fitting data
#' Eq01 <- lm(Yield ~ Age, data=dat)
#'
#' # independent/validation data
#' test_data <- data.frame(Age=1:50, Yield=2.5*Age^1.4)
#'
#' # predict with the model i.e. Eq01, using the independent/validation data
#' test_data$pred.Yield <- predict(Eq01, test_data)
#'
#' # Evaluate the quality of the prediction using the 'valstat' function.
#' # You need the observed and predicted values. Specify the number of parameters in the model.
#'
#' valstat(test_data$Yield, test_data$pred.Yield, 2)

valstat <- function(obs.y, pred.y, p){
  Residual <- obs.y - pred.y
  n <- length(Residual)
  SSR <- sum(Residual^2)
  SST <- sum((obs.y - mean(obs.y))^2)
  TRE <- sum(abs(Residual)) / sum(pred.y)
  Bias <- mean(Residual)
  MRB <- mean((pred.y-obs.y)/obs.y)
  AbsBias <- mean(abs(Residual))
  MAPE <- (mean(abs(Residual)/obs.y))*100
  MSE <- SSR / n
  RMSE <- sqrt(SSR / n)
  Percent.RMSE <- (RMSE / mean(obs.y)) * 100
  R2 <- 1 - SSR / SST
  R2adj <- 1 - SSR / SST * (n - 1) / (n - p)
  APC <- ((n + p) / (n - p))*(1 - (R2))
  logL <- 0.5 * (-n * (log(2 * pi) + 1 - log(n) + log(SSR)))
  AIC <- 2 * (p + 1) - 2 * logL
  AICc <- AIC + (2*p*(p+1))/(n-p-1)
  BIC <- (p + 1) * log(n) - 2 * logL
  HQC <- 2 * p * log(log(n)) - 2 * logL

  return(data.frame(n = n, SSR = SSR, TRE = TRE, Bias = Bias, MRB = MRB, MAB = AbsBias, MAPE = MAPE,
                    MSE = MSE, RMSE = RMSE, Percent.RMSE = Percent.RMSE,
                    R2 = R2, R2adj = R2adj, APC = APC, logL = logL, AIC = AIC, AICc = AICc, BIC = BIC, HQC = HQC))
}
