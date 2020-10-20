
#' @import Rcpp
#' @export
#' @description Main function to solve efficiently and quickly a symmetric bandlinear system. Theses systems are solved much faster than standards system, dropping from complexity O(n³) to O(0.5*nk²), where k is the number of sub diagonal.
#' @title bandsolve
#' @param A Band square matrix in rotated form. The rotated form can be obtained with the function as.rotated: it's the visual rotation by 90 degrees of the matrix, where subdiagonal are discarded.
#' @param b right hand side of the equation. Can be either a vector or a matrix. If not supplied, the function return the inverse of A.
#' @param inplace Should results overwrite pre-existing data? Default set to false.
#' @return Solution of the linear problem.
#' @examples
#'

#' A = diag(4)
#' A[2,3] = 2
#' A[3,2] = 2
#' R = mat2rot(A)
#' solve(A)
#' bandsolve(R)
#'
#' set.seed(100)
#'
#' n = 1000
#' D0 = rep(1.25, n)
#' D1 = rep(-0.5, n-1)
#' b = rnorm(n)

bandsolve <- function(A, b = NULL, inplace = FALSE) {
    if ((nrow(A) == ncol(A)) & (A[nrow(A), ncol(A)] != 0))
        stop("A should be a rotated matrix!")
    if (A[nrow(A), 2] != 0)
        stop("A should be a rotated matrix!")
    if (is.vector(b)) {
        if (length(b) != nrow(A))
            stop("Dimension problem")
        if (inplace) {
            return(bandsolve_cpp(A, as.matrix(b))$x)
        } else {
            Amem = matrix(NA, nrow(A), ncol(A))
            Amem[] = A[]
            bmem = rep(NA, length(b))
            bmem[] = b[]
            return(bandsolve_cpp(Amem, as.matrix(bmem))$x)
        }
    } else if (is.matrix(b)) {
        if (nrow(b) != nrow(A))
            stop("Dimension problem")
        if (inplace) {
            return(bandsolve_cpp(A, b)$x)
        } else {
            Amem = matrix(NA, nrow(A), ncol(A))
            Amem = A[]
            Bmem = matrix(NA, nrow(b), ncol(b))
            Bmem[] = b[]
            return(bandsolve_cpp(Amem, Bmem)$x)
        }
    } else if (is.null(b)) {
        B = diag(nrow(A))
        if (inplace) {
            return(bandsolve_cpp(A, B)$x)
        } else {
            Amem = matrix(NA, nrow(A), ncol(A))
            Amem[] = A[]
            return(x = bandsolve_cpp(Amem, B)[[2]])
        }
    } else {
        stop("b must either be a vector or a matrix")
    }
}


#' @useDynLib aspline
#'
#' @title Get back a symmetric square matrix based on his rotated row-wised version.
#' @description Get back a symmetric square matrix based on his rotated row-wised version.
#' The rotated form of the input is such each column correspond to a diagonal, where the first column is the main diagonal and next ones are the upper/lower-diagonal.
#' To match dimension, last element of these columns are discarded.
#' @param R Rotated matrix.
#' @return Band square matrix.
#' @examples
#'
#' D0 = 1:5;
#' D1 = c(0,1,0,0);
#' D2 = rep(2,3);
#'
#' A = rot2mat(cbind(D0,c(D1,0),c(D2,0,0)))
#' A
#' mat2rot(rot2mat(cbind(D0,c(D1,0),c(D2,0,0))))
#'
#' @export
rot2mat <- function(R) {
    N = nrow(R)
    l = ncol(R)
    M = matrix(0, N, N)
    if (l > 1) {
        for (i in 1:(l - 1)) {
            if (i == N - 1) {
                M[-c(1:i), -c(N:(N - i + 1))] = R[-c(N:(N - i + 1)), i + 1]
                M[-c(N:(N - i + 1)), -c(1:i)] = R[-c(N:(N - i + 1)), i + 1]
            } else {
                diag(M[-c(1:i), -c(N:(N - i + 1))]) = R[-c(N:(N - i + 1)), i + 1]
                diag(M[-c(N:(N - i + 1)), -c(1:i)]) = R[-c(N:(N - i + 1)), i + 1]
            }
        }
        diag(M) = R[, 1]
    }
    return(M)
}


#' @useDynLib aspline
#'
#' @import Rcpp
#' @export
#' @description Rotate a symmetric band matrix to get the rotated matrix associated.
#' Each column of the rotated matrix correspond to a diagonal. The first column is the main diagonal, the second one is the upper-diagonal and so on.
#' Artificial 0 are placed at the end of each column if necessary.
#' @title Rotate a band matrix to get the rotated row-wised matrix associated.
#' @param M Band square matrix or a list of diagonal.
#' @return Rotated matrix.
#' @examples
#'
#' A = diag(4)
#' A[2,3] = 2
#' A[3,2] = 2
#'
#' ## Original Matrix
#' A
#' ## Rotated version
#' R = mat2rot(A)
#' R
#'
#' rot2mat(mat2rot(A))

mat2rot <- function(M) {
  if (is.matrix(M)) {
    N = ncol(M)
    l = 0
    for (i in 1:N) {
      lprime = which(M[i, ] != 0)
      if (lprime[length(lprime)] - i > l)
        l = lprime[length(lprime)] - i
    }
    R = matrix(0, N, l + 1)
    R[, 1] = diag(M)
    if (l > 0) {
      for (j in 1:l) {
        if (j == N - 1) {
          R[, 1 + j] = c(M[-c(N:(N - j + 1)), -c(1:j)], rep(0, j))
        } else {
          R[, 1 + j] = c(diag(M[-c(N:(N - j + 1)), -c(1:j)]), rep(0, j))
        }
      }
    }
    return(R)
  } else if (is.list(M)) {
    R = matrix(0, length(M[[1]]), length(M))
    for (i in 1:length(M)) {
      R[, i] = c(M[[i]], rep(0, i - 1))
    }
    return(R)
  }
}
