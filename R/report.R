#' Determine which results are suspicious
#'
#' \lifecycle{maturing}
#' The \link{GWAS} function writes all results, both valid and
#' invalid, to a log file. This function uses heuristics to try to
#' classify rows as suspicious or unsuspicious.  The
#' \link{loadResults} function returns unsuspicious rows while
#' \link{loadSuspicious} returns suspicious rows.
#'
#' Is it not recommended to call \code{isSuspicious} directly because
#' it is already called by \link{loadResults} and \link{loadSuspicious}.
#'
#' @details OpenMx reports exceptions in the \sQuote{catch1}
#'   column. Any error message in the \sQuote{catch1} column is
#'   suspicious. Any optimizer status code besides \sQuote{OK} is
#'   suspicious. It is suspicious if the focal parameter or its
#'   standard error is \code{NA}. If \sQuote{signAdj} was requested
#'   and it is \code{NA} then suspicion is also aroused.  A mean and
#'   standard deviation are computed for each available parameter.
#' Any estimate with an absolute Z score larger than 10 is suspicious.
#'
#' @param got output passed from \link{loadResults}
#' @param pars names of the parameters available in \code{got}
#' @return
#' a vector of logicals for each row of \code{got} indicating suspicion (if TRUE)
#'
#' @family reporting
#' @export
#' @importFrom stats sd
#' @examples
#' tdir <- tempdir()
#' dir <- system.file("extdata", package = "gwsem")
#' pheno <- data.frame(anxiety=rnorm(500))
#' m1 <- buildItem(pheno, 'anxiety')
#' GWAS(m1, file.path(dir,"example.pgen"),
#'     file.path(tdir,"out.log"))
#' loadResults(file.path(tdir,"out.log"), "snp2anxiety")  # called in here
isSuspicious <- function(got, pars) {
  # minMAF? TODO
  mask <- (!is.na(got[['catch1']]) & got[['catch1']] != '') | got[['statusCode']] != 'OK'
  if (nrow(got) > 2) {
    for (p1 in pars) {
      mask <- mask | is.na(got[[p1]])
      if (!is.null(got[[paste0(p1,'SE')]])) {
        mask <- mask | is.na(got[[paste0(p1,'SE')]])
      }
    }
  }
  mask
}

#' Load GWAS results into a single data.frame
#'
#' \lifecycle{maturing}
#' The \code{signAdj} column is important and not optional for latent
#' factor models.  Loadings to factor indicators can take any sign. If
#' your focus is the regression from the SNP to the factor then this
#' regression estimate will need to be multiplied by the sign of one
#' of the factor loadings. Pick a loading associated with a strong
#' indicator of the factor.
#'
#' A1 is the reference allele and A2 is the alternate allele.
#'
#' Two columns are added, \code{Z} and \code{P}. \code{Z} is the
#' focal parameter divded by its standard error. \code{P} is the
#' unadjusted two-sided normal CDF corresponding to the absolute
#' \code{Z} score.
#'
#' Suspicious rows are excluded.
#'
#' @param path vector of paths to result files created by \link{GWAS}
#' @param focus parameter name on which to calculate a Z score and p-value
#' @param extraColumns character vector of additional columns to load
#' @param .retainSE logical. Keep a column for the SE of the focus parameter
#' @param signAdj name of column. Value of focus parameter is multiplied by the sign of the named column
#' @param moderatorLevel see details
#' @template args-dots-barrier
#' @family reporting
#' @export
#' @importFrom data.table fread
#' @importFrom stats pnorm
#' @examples
#' tdir <- tempdir()
#' dir <- system.file("extdata", package = "gwsem")
#' pheno <- data.frame(anxiety=rnorm(500))
#' m1 <- buildItem(pheno, 'anxiety')
#' GWAS(m1, file.path(dir,"example.pgen"),
#'     file.path(tdir,"out.log"))
#' loadResults(file.path(tdir,"out.log"), "snp2anxiety")
loadResults <- function(path, focus, ..., extraColumns=c(),
			.retainSE=FALSE, signAdj=NULL, moderatorLevel=NULL) {
  if (length(list(...)) > 0) stop("Rejected are any values passed in the '...' argument")
  sel <- c('MxComputeLoop1', 'CHR','BP','SNP','A1','A2','statusCode','catch1',
	   focus,paste0(focus,'SE'), extraColumns, signAdj)
  mainEffect <- c()
  if (!is.null(moderatorLevel)) {
    got <- regexpr('^snp_(\\w+)_to_(\\w+)$', focus, perl=TRUE)
    if (got == -1) {
      stop(paste("When moderatorLevel requested, focus must be of the form",
                 "snp_XXX_to_YYY not", focus))
    } else {
      cstart <- attr(got,"capture.start")
      clen <- attr(got,"capture.length")
      mod <- substr(focus, cstart[1], cstart[1]+clen[1]-1L)
      outcome <- substr(focus, cstart[2], cstart[2]+clen[2]-1L)
      mainEffect <- paste0('snp_to_', outcome)
      sel <- c(sel, mainEffect, paste0(mainEffect, 'SE'), paste0('V', focus, ':', mainEffect))
    }
  }
  got <- list()
  for (p1 in path) {
    d1 <- fread(p1, header=TRUE,
                sep="\t", check.names=FALSE, quote="", select = sel)
    d1 <- d1[!isSuspicious(d1, c(focus, signAdj, mainEffect)),]
    if (!is.null(signAdj)) {
      d1[[focus]] <- d1[[focus]] * sign(d1[[signAdj]])
      if (!(signAdj %in% extraColumns)) d1[[signAdj]] <- NULL
    }
    if (!is.null(moderatorLevel)) {
      d1[[focus]] <- d1[[mainEffect]] + moderatorLevel * d1[[focus]]
      d1[[paste0(focus,'SE')]] <-
        sqrt(d1[[paste0(mainEffect,'SE')]]^2 +
               moderatorLevel^2 * d1[[paste0(focus,'SE')]]^2 +
               2*moderatorLevel * d1[[paste0('V', focus, ':', mainEffect)]])
      for (col in c(mainEffect, paste0(mainEffect, 'SE'), paste0('V', focus, ':', mainEffect))) {
        d1[[col]] <- NULL
      }
    }
    got <- rbind(got, d1)
  }
  got$Z <- got[[focus]] / got[[paste0(focus,'SE')]]
  if (!.retainSE) got[[paste0(focus,'SE')]] <- NULL # redundent; save RAM
  got$P <- 2*pnorm(-abs(got$Z))
  attr(got, 'focus') <- focus
  class(got) <- c("gwsemResult", class(got))
  got
}

#' Load suspicious GWAS results into a single data.frame
#'
#' \lifecycle{maturing}
#' These rows are excluded by \link{loadResults}.
#'
#' @param path vector of paths to result files created by \link{GWAS}
#' @param focus parameter name on which to calculate a Z score and p-value
#' @param extraColumns character vector of additional columns to load
#' @param signAdj name of column. Value of focus parameter is multiplied by the sign of the named column
#' @param moderatorLevel see details
#' @template args-dots-barrier
#' @family reporting
#' @export
#' @importFrom data.table fread
#' @examples
#' tdir <- tempdir()
#' dir <- system.file("extdata", package = "gwsem")
#' pheno <- data.frame(anxiety=rnorm(500))
#' m1 <- buildItem(pheno, 'anxiety')
#' GWAS(m1, file.path(dir,"example.pgen"),
#'     file.path(tdir,"out.log"))
#' loadSuspicious(file.path(tdir,"out.log"), "snp2anxiety")
loadSuspicious <- function(path, focus, ..., extraColumns=c(),
			signAdj=NULL, moderatorLevel=NULL) {
  if (length(list(...)) > 0) stop("Rejected are any values passed in the '...' argument")
  sel <- c('MxComputeLoop1', 'CHR','BP','SNP','A1','A2','statusCode','catch1',
	   focus,paste0(focus,'SE'), extraColumns, signAdj)
  mainEffect <- c()
  if (!is.null(moderatorLevel)) {
    got <- regexpr('^snp_(\\w+)_to_(\\w+)$', focus, perl=TRUE)
    if (got == -1) {
      stop(paste("When moderatorLevel requested, focus must be of the form",
                 "snp_XXX_to_YYY not", focus))
    } else {
      cstart <- attr(got,"capture.start")
      clen <- attr(got,"capture.length")
      mod <- substr(focus, cstart[1], cstart[1]+clen[1]-1L)
      outcome <- substr(focus, cstart[2], cstart[2]+clen[2]-1L)
      mainEffect <- paste0('snp_to_', outcome)
      sel <- c(sel, mainEffect, paste0(mainEffect, 'SE'), paste0('V', focus, ':', mainEffect))
    }
  }
  got <- list()
  for (p1 in path) {
    d1 <- fread(p1, header=TRUE,
                sep="\t", check.names=FALSE, quote="", select = sel)
    d1 <- d1[isSuspicious(d1, c(focus, signAdj)),]
    if (!is.null(signAdj)) {
      d1[[focus]] <- d1[[focus]] * sign(d1[[signAdj]])
    }
    if (!is.null(moderatorLevel)) {
      d1[[focus]] <- d1[[mainEffect]] + moderatorLevel * d1[[focus]]
      d1[[paste0(focus,'SE')]] <-
        sqrt(d1[[paste0(mainEffect,'SE')]]^2 +
               moderatorLevel^2 * d1[[paste0(focus,'SE')]]^2 +
               2*moderatorLevel * d1[[paste0('V', focus, ':', mainEffect)]])
    }
    got <- rbind(got, d1)
  }
  got$Z <- got[[focus]] / got[[paste0(focus,'SE')]]
  got$P <- 2*pnorm(-abs(got$Z))
  attr(got, 'focus') <- focus
  got
}

#' Creates a Manhattan plot
#'
#' Uses the qqman package to create a Manhattan plot.
#'
#' @param x the result of \link{loadResults}
#' @param y an extra argument that should not be used
#' @param ... arguments forwarded to \link[qqman]{manhattan}
#' @export
#' @importFrom qqman manhattan
#' @family reporting
#' @return
#' A Manhattan plot.
#' @examples
#' tdir <- tempdir()
#' dir <- system.file("extdata", package = "gwsem")
#' pheno <- data.frame(anxiety=rnorm(500))
#' m1 <- buildItem(pheno, 'anxiety')
#' GWAS(m1, file.path(dir,"example.pgen"),
#'     file.path(tdir,"out.log"))
#' got <- loadResults(file.path(tdir,"out.log"), "snp_to_anxiety")
#' plot(got)
plot.gwsemResult <- function(x, y, ...) {
	if (!missing(y)) stop("plot does not accept a y= argument")
	x$P[is.na(x$P)] <- 1
	manhattan(x, ...)
}
