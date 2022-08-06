# x The fiiirstones arguments
# num The firrrstones arguments
# k The firssssstones arguments
noceros <- function(x, num = TRUE, k = 0) {

  nn <- length(which(x > k))

  if (num) {
    nn

  } else {
    if (nn > 0) {
      which(x > k)
    } else {
      NULL
    }
  }
}

# datos The firstones arguments
# k The secondones arguments
# noceros
sinceros <-
  function (datos, k) {
    datos = as.matrix(datos)
    datos0 <- as.matrix(datos)

    if (is.null(k)) {

      mini0 <- min(datos[noceros(datos, num = FALSE, k = 0)])

      kc <- mini0/2

      datos0[datos0 == 0] <- kc

    } else {

      datos0[datos0 == 0] <- k

    }

    datos0
  }

# object The firstones argumentss
# method The secondones argumentss
# refColumn The thirdones argumentss
# logratioTrim The fourthones argumentss
# sumTrim The fifthones argumentss
# doWeighting The sixones argumentss
# Acutoff The sevenones argumentss
# quantile The eightones argumentss
.calcNormFactors <- function(object, method=c("TMM","quantile"), refColumn=NULL,
                             logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10,
                             quantile=0.75) {

  method <- match.arg(method)

  if( is.matrix(object) ) {
    if(is.null(refColumn))
      refColumn <- 1
    data <- object
    libsize <- colSums(data)
  } else {
    stop("calcNormFactors() only operates on 'matrix' objects")
  }

  f <- switch(method,
              TMM = apply(data,2,.calcFactorWeighted,ref=data[,refColumn],
                          logratioTrim=logratioTrim, sumTrim=sumTrim, doWeighting=doWeighting,
                          Acutoff=Acutoff),
              quantile = .calcFactorQuantile(data, libsize, q=quantile))

  f <- f/exp(mean(log(f)))

  return(f)

}

# data The firstones argumentsss
# lib.size The firstones argumentssss
# q The firstones argumentssssss
.calcFactorQuantile <- function (data, lib.size, q=0.75) {
  y <- t(t(data)/lib.size)
  f <- apply(y,2,function(x) quantile(x,p=q))
  f/exp(mean(log(f)))
}

# obs The firstonesss argumentss
# ref The secondonesss argumentss
# refColumn The thirdonesss argumentsss
# logratioTrim The fourthssonesss argumentss
# sumTrim The fivethssonesss argumentsssssss
# doWeighting The sixssonesss argumentsssssss
# Acutoff The sixssonesss argumentssssssssssss
.calcFactorWeighted <- function(obs, ref, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10) {

  if( all(obs==ref) )
    return(1)

  obs <- as.numeric(obs)
  ref <- as.numeric(ref)

  nO <- sum(obs)
  nR <- sum(ref)
  logR <- log2((obs/nO)/(ref/nR))          # log ratio of expression, accounting for library size
  absE <- (log2(obs/nO) + log2(ref/nR))/2  # absolute expression
  v <- (nO-obs)/nO/obs + (nR-ref)/nR/ref   # estimated asymptotic variance

  # remove infinite values, cutoff based on A
  fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)

  logR <- logR[fin]
  absE <- absE[fin]
  v <- v[fin]

  # taken from the original mean() function
  n <- sum(fin)
  loL <- floor(n * logratioTrim) + 1
  hiL <- n + 1 - loL
  loS <- floor(n * sumTrim) + 1
  hiS <- n + 1 - loS

  #keep <- (rank(logR) %in% loL:hiL) & (rank(absE) %in% loS:hiS)
  # a fix from leonardo ivan almonacid cardenas, since rank() can return
  # non-integer values when there are a lot of ties
  keep <- (rank(logR)>=loL & rank(logR)<=hiL) & (rank(absE)>=loS & rank(absE)<=hiS)

  if (doWeighting)
    2^( sum(logR[keep]/v[keep], na.rm=TRUE) / sum(1/v[keep], na.rm=TRUE) )
  else
    2^( mean(logR[keep], na.rm=TRUE) )
}

# datos The firstonesss argumentssssssss
# long The secondonessssss argumentsssss
# lc The thirdonesssssss argumentssssss
# k The fourthssonessssss argumentsssssss
# refColumn The fivethssssonesss argumentsssssss
# logratioTrim The sixssonesss argumentsssssssfff
# sumTrim The sevessonesss argumentssssssssssss
# doWeighting The eightssonesss argumentssssssssssss
# Acutoff The ninesonesss argumentssssssssssss
# sinceros
# .calcNormFactors
tmm <- function (datos, long = 1000, lc = 0, k = 0, refColumn = 1,
                 logratioTrim = .3, sumTrim = 0.05, doWeighting = TRUE,
                 Acutoff = -1e10) {

  # lc: Length correction. Expression is divided by long^lc. lc can be any real number.

  if (!is.null(ncol(long))) {
    mynames = long[,1]
    long = long[,2]
    names(long) = mynames
  }

  L <- (long/1000)^lc
  datos = datos/L

  total <- colSums(as.matrix(datos))

  datos0 <- sinceros(datos, k)

  if (ncol(as.matrix(datos)) > 1) {

    fk <- .calcNormFactors(as.matrix(datos), refColumn = refColumn, method = "TMM",
                           logratioTrim = logratioTrim, sumTrim = sumTrim,
                           doWeighting = doWeighting, Acutoff = Acutoff)

    fk = fk * (total/mean(total))

    datos.norm <- t(t(datos0)/fk)

  } else {

    datos.norm <- datos0/L

  }

  na.omit(datos.norm)

}
