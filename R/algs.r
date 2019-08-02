
# Transition matrix for the example model.
P <- matrix(c(
   0,   0,   0,   0,   1,
   1,   0,   0, 1/4,   0,
   0, 1/2, 1/2,   0,   0,
   0, 1/2, 1/4,   0,   0,
   0,   0, 1/4, 3/4,   0), nrow = 5)

n <- nrow(P)
Q <- P[-n, -n]

S <- matrix(c(
   1,   0,   0, 1/4,
   0, 1/2, 1/2,   0,
   0, 1/2, 1/4,   0,
   0,   0, 1/4, 1/2,
   0,   0,   0, 1/4), nrow = 4)


######################################################################
# Algorithms
######################################################################

# Compute the non-terminal expectation and
# variance matrices for the stochastic matrix
# P.
#
# P: a row-stochastic matrix
# return: list(N,V) where
#   N: the expected occurrence of each state
#   V: the associated variances
get_nte <- function(P) {
	n <- nrow(P)
	Q <- P[-n, -n]
	N <- solve(diag(n - 1) - Q)
	V <- N %*% (2 * diag(diag(N)) - diag(n - 1)) - N * N
	list(N, V)
}

# Compute the Perron eigenvector for the
# stochastic matrix P and return it. The
# computation is performed by computation
# of the fundamental matrix.
#
# P: a square row-stochastic matrix
# N: the fundamental matrix, if available
# return: the Perron eigenvector
get_pi <- function(P, N) {
	n <- nrow(P)
	if (nargs() < 2) {
		N <- get_nte(P)[[1]]
	}
	len <- 1 + sum(N[1, ])
	pi <- N[1, ]/len
	append(pi, c(1 - sum(pi)))
}

# Compute the Perron eigenvector for the row
# stochastic matrix P and return it.  The
# computation is performed using the power
# method.
#
# P: a row-stochastic matrix
# g: (optional) an initial guess
# return: list(y,step) where
#   y: an approximation of the Perron eigenvector
#   step: the number of iterations required
get_pi_approx <- function(P, g) {
	n <- nrow(P)
	d <- n + 1
	P <- rbind(cbind(P, 0), 0)
	P[n, d] = 1
	P[d, d] = 1/2
	P[d, 1] = 1/2
	P[n, 1] = 0
	if (missing(g)) {
		yold <- rep(1/d, d)
	} else {
		g <- c(g, g[n] * 2)
		yold <- g/sum(g)
	}
	y <- yold %*% P
	step <- 1
	limit <- 1e-06
	while (sum(abs(yold - y)) > limit) {
		step = step + 1
		yold = y
		y = y %*% P
	}
	y = y[1:5]
	y = y/sum(y)
	list(y, step)
}

# Compute the matrix of arc sensitivities.
# The first column of the returned matrix
# is the source state, the second column is
# the target state. The remaining columns
# are the changes in the occupancies.
#
# P: the state transition matrix
# pi: the (optional) long run occupancies
# return: the sensitivities matrix
get_sensitivities <- function(P, pi) {
	n <- nrow(P)
	if (missing(pi)) {
		pi <- get_pi_approx(P)
	}
	Z <- matrix(rep(c(0), (n + 2)), nrow = 1)
	m <- 1
	for (i in 1:n) {
		for (j in 1:n) {
			if (0 < P[i, j] && P[i, j] < 1) {
				t <- P[i, ]
				x <- 1 - P[i, j]
				P[i, ] = t/x * 0.05
				P[i, j] = 0.95
				ph <- get_pi_approx(P, pi)[[1]]
				P[i, ] = t/x * 0.95
				P[i, j] = 0.05
				pl <- get_pi_approx(P, pi)[[1]]
				Z = rbind(Z, 0)
				Z[m, 1:2] = c(i, j)
				Z[m, 3:(n + 2)] = (ph - pl)/0.9
				m = m + 1
				P[i, ] = t
			}
		}
	}
	Z[1:m - 1, ]
}
