mscorev <-
function (ymat, inner = 0, trim = 3, lambda = 0.5)
{
# Similar to the mscorev function from sensitivitymv version 1.3
    stopifnot((inner>=0)&(inner<=trim))
    stopifnot((lambda>0)&(lambda<1))
    if (is.data.frame(ymat)) ymat <- as.matrix(ymat)
    stopifnot(is.matrix(ymat))
    stopifnot(2<=dim(ymat)[2])
    stopifnot(all(!is.na(as.vector(ymat[,1:2]))))
    n <- dim(ymat)[1]
    m <- dim(ymat)[2]
    out <- matrix(NA, n, m)
    one <- rep(1, m - 1)
    difs <- array(NA, c(n, m, m - 1))
    TonT <- FALSE
    qu <- lambda
    for (j in 1:m) {
        difs[, j, ] <- outer(as.vector(unlist(ymat[, j])), one,
            "*") - ymat[, -j]
    }
    ms <- as.vector(difs)
    if ((trim < Inf) | (inner > 0)) {
        hqu <- as.numeric(quantile(abs(ms), qu, na.rm = TRUE))
        if (hqu > 0) {
            ms <- ms/hqu
            if ((trim < Inf) & (inner < trim)) {
                ab <- pmin(1, pmax(0, (abs(ms) - inner))/(trim -
                  inner))
            }
            else if ((trim < Inf) & (inner == trim)) {
                ab <- 1 * (abs(ms) > inner)
            }
            else {
                ab <- pmax(0, abs(ms) - inner)
            }
            ms <- sign(ms) * ab
        }
        else {
            warning("Error: Scale factor is zero.  Increase lambda.")
        }
    }
    ms <- array(ms, c(n, m, m - 1))
    ms <- apply(ms, c(1, 2), sum, na.rm = TRUE)
    ms[is.na(ymat)] <- NA
    colnames(ms) <- colnames(ymat)
    ni <- apply(!is.na(ymat), 1, sum)
    use <- (ni >= 2) & (!is.na(ms[, 1]))
    ms <- ms[use, ]
    ni <- ni[use]
    if (TonT) {
        ms <- (ms/outer(ni - 1, rep(1, m), "*"))/(dim(ms)[1])
    }
    else {
        ms <- ms/outer(ni, rep(1, m), "*")
    }
    ms
}
