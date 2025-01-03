Bases.NonNulles <- function(F.bases) {
  F.bases.Present <- c()
  M <- dim(F.bases)[2]
  PresentBases <- c()
  for (j in 1:M) {
    fj <- F.bases[, j]
    if (sum(abs(fj)) > 10^-10) {
      PresentBases <- c(PresentBases, j)
    }
  }
  F.bases.Present <- F.bases[, PresentBases]
  return(list(Fh.P = F.bases.Present, PresentBases = PresentBases))
}

CreationBases <- function(position, keep = NULL) {
  Ff <- c()

  n <- max(position)
  Ff <- c()
  for (i in seq(1, n, 1)) {
    Ff <- cbind(Ff, sin(2 * pi * i * position / n), cos(2 * pi * i * position / n))
  }

  x <- position / n
  Fx <- c()
  rg <- seq(0.1, 2, 0.02)
  for (i in 1:length(rg)) {
    Fx <- cbind(Fx, x^rg[i])
  }

  F.tot <- cbind(Ff, Fx)

  if (!is.null(keep)) {
    F.tot <- F.tot[, keep]
  }

  F.Bases <- c()
  F.Bases.tot <- c()
  F.Bases.tot <- Bases.NonNulles(F.tot)
  Num.Bases.Pres <- F.Bases.tot$PresentBases
  F.Bases <- F.Bases.tot$Fh.P

  return(list(F.Bases = F.Bases, Num.Bases.Pres = Num.Bases.Pres))
}
