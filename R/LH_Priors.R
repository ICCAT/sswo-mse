#' Sample Life History Parameters for the Taylor Priors
#'
#' Generates `n` samples from a truncated log-normal distribution from `Mean` and
#' `CV` values located in the CSV specified in `path`
#'
#' @param n The number of samples
#' @param truncSD The number of standard deviations to truncate
#' @param path The full file path for the CSV file with the `Mean` and `SD` values. Assumes
#' the first column are the names of the parameters
#'
#' @return A data.frame of the sampled values
#' @export
#'
Sample_Trunc_LogNormal_Taylor <- function(n=10000, truncSD=1.96,
                                          path="Resources/Taylor_Values.csv") {

  lh.pars <- read.csv(path)
  neg_ind <- which(lh.pars$Mean<0)
  lh.pars$Mean[neg_ind] <- -lh.pars$Mean[neg_ind]
  lh.pars$mu <-  log(lh.pars$Mean)-0.5* lh.pars$CV^2
  lh.pars$a <- lh.pars$mu - lh.pars$CV * truncSD
  lh.pars$b <- lh.pars$mu + lh.pars$CV * truncSD

  # reverse sign
  vals <- truncnorm::rtruncnorm(n*nrow(lh.pars), a=lh.pars$a, b=lh.pars$b, mean=lh.pars$mu, sd=lh.pars$CV)
  vals <- matrix(vals, nrow=n, ncol=nrow(lh.pars), byrow=TRUE) |>
    exp() |> as.data.frame()
  colnames(vals) <- lh.pars[,1]
  if (length(neg_ind)>0)
    vals[,neg_ind] <- -vals[,neg_ind]
  vals
}

?tmvtnorm::rtmvnorm()

#' Extract predicted values from FishLife
#'
#' @param Genus Genus
#' @param Species Species
#' @param level Prediction level to extract: 1:
#'
#' @return A list object returned from `FishLife::Plot_taxa`
#' @export
get_FishLife_Predictions <- function(Genus="Xiphias", Species="gladius", level=2) {
  swo.fishlife <-FishLife::Plot_taxa(FishLife::Search_species(
    Genus=Genus,
    Species=Species)$match_taxonomy,
    mfrow=c(3,2),
  )[[level]]
}


#' @describeIn get_FishLife_Predictions Calculate the correlation matrix
#' @param Cov_pred `Cov_pred` object returned by `get_FishLife_Predictions`
#' @export
calc_Correlation_Matrix <- function(Cov_pred) {
  cor <- stats::cov2cor(Cov_pred)   #extract the correlation matrix from the covariance matrix
  diag(cor) <- NA
  cor <- cor |>  dplyr::as_tibble() |>
    dplyr::mutate(term=colnames(cor)) |>
    dplyr::relocate(term)
  #have to typecast to use focus because swo.fishlife.cor2 is not an object derived from the corr function
  class(cor) <-c("cor_df", "tbl_df", "tbl", "data.frame")
  cor
}


#' Add correlation structure to selected values
#'
#' @param simulated.data.trunc A data.frame returned by `Sample_Trunc_LogNormal_Taylor`
#' @param cor_matrix A correlation matrix returned by `calc_Correlation_Matrix`
#' @param Var_Names A data.frame with the required variable names from
#' `calc_Correlation_Matrix` in the first column, and matching variable names from
#' `Sample_Trunc_LogNormal_Taylor` in the second column
#' @param truncSD The number of standard deviations to truncate
#'
#' @return A data.frame with values from `simulated.data.trunc` updated
#' with the correlation structure from `cor_matrix`
#' @export
Add_Correlation_Taylor <- function(simulated.data.trunc,
                                   cor_matrix=NULL,
                                   Var_Names=NULL,
                                   truncSD=1.96)  {
  red.cor.mat <- cor_matrix |>
    corrr::focus(Var_Names[,1], mirror = TRUE) |>
    dplyr::arrange(match(term, Var_Names[,1])) |>
    dplyr::select(-term) |>
    as.matrix()

  diag(red.cor.mat)<-1

  red.sim.data <- simulated.data.trunc |>
    dplyr::select(dplyr::all_of(Var_Names[,2])) |>
    log()

  red.sim.dat.means<-colMeans(red.sim.data)
  red.sim.dat.sds<-sapply(red.sim.data, sd)

  sim.data.cov <- red.cor.mat *as.matrix(red.sim.dat.sds) %*% t(as.matrix(red.sim.dat.sds))
  row.names(sim.data.cov)<-colnames(sim.data.cov)

  #draw from a multivariate distribution of M, Linf, k, L50
  lower <- red.sim.dat.means - red.sim.dat.sds*truncSD
  upper <- red.sim.dat.means + red.sim.dat.sds*truncSD
  #simulate from a truncated multivariate normal distribution in log space
  red.simdata.t.mvrnorm <- tmvtnorm::rtmvnorm(length(simulated.data.trunc[,1]),
                                              mean = red.sim.dat.means,
                                              sigma = sim.data.cov,
                                              lower=lower,
                                              upper=upper
  ) |> exp() |> as.data.frame()

  colnames(red.simdata.t.mvrnorm)<-colnames(red.sim.data)


  simulated.data.trunc_updated <- simulated.data.trunc
  for (nm in names(red.simdata.t.mvrnorm)) {
    simulated.data.trunc_updated[,nm] <-  simulated.data.trunc[,nm]
  }
  simulated.data.trunc_updated

}

#' Calculate Steepness using the Mangel (2010) method
#'
#' @param x sample number (up to `nrow(DF)`)
#' @param DF A data.frame with `nrow` samples and specifically named columns, see `Details`
#' @param plus.group Age of the plus group
#'
#' @details
#' `DF` must have the following column names:
#' - `mort`: Natural mortality rate
#' - `vLinf`:
#' - `vonK`:
#' - `vtto`:
#' - `length.weight.alpha`:
#' - `length.weight.beta`:
#' - `L50`:
#' - `L50_95`:
#' - `BFConstant`:
#' - `Bfa`:
#' - `BFb`:
#' - `ParChange`:
#' - `RealizedFecFactor`:
#' - `NS`:
#' - `Se`:
#'
#' @return the predicted steepness value
#' @export
Calculate_Steepness <- function(x, DF,
                                plus.group=16) {
  ages <- 1:plus.group
  S <- c(1,exp(-DF$mort[x]*ages[1:(plus.group-1)]))

  N <- 1* S
  N[plus.group]<-N[plus.group-1]*S[2]/(1-S[2])

  la<- DF$vLinf[x]*(1-exp(-DF$vonK[x] *(ages-DF$vtto[x])))
  wa<- DF$length.weight.alpha[x]*la^DF$length.weight.beta[x]

  a50 <- iVB(DF$vtto[x], DF$vonK[x], DF$vLinf[x], DF$L50[x])
  a95 <- iVB(DF$vtto[x], DF$vonK[x], DF$vLinf[x], DF$L50[x]+DF$L50_95[x])
  mat <- 1/(1 + exp(-log(19) * ((ages - a50)/(a95-a50))))
  eggs<-N*wa*mat

  batch.fecundity<- DF$BFConstant[x] + DF$Bfa[x] * ages^(DF$BFb[x]*DF$ParChange[x])
  realized.fec<- DF$RealizedFecFactor[x] *mat*batch.fecundity
  SPR<-sum(eggs)

  NS.E.over.W<- DF$NS[x]*sum(realized.fec)/sum(wa)
  alpha<-NS.E.over.W*DF$Se[x]
  steepness<-alpha*SPR/(4+alpha*SPR)
  steepness
}

#' @describeIn Calculate_Steepness Calculate inverse von Bert
#' @export
iVB <- function (t0, K, Linf, L) {
  max(1, ((-log(1 - L/Linf))/K + t0))
}
