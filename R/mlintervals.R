
## UTILITTY FUNCTIONS

## function to extract the coefficients of all the parabolas fitted. returns a nx3 matrix
.parabolas= function(x) {

	M= x$data
	n= length(M)
	knots= x$knots

	first= 1.5*M[1]-.5*knots[1]
	last= 1.5*M[n]-.5*knots[n-1]
	knots= c(first, knots, last)

	a= knots[-n-1]
	b= knots[-1]

	p0= a
	p1= 6*M-4*a-2*b
	p2= -6*M+3*a+3*b

	return(cbind(p0, p1, p2))
}


############################
## CLASS interv_fit

#' Interpolating the process, given an interval averaged series
#'
#' Maximum Likelihood interpolation of a Wiener process, given a series of observed integrals
#' over adiacent time intervals of equal lenght. To get point estimates of the interpolating
#' function, use method \code{\link[=predict.interv_fit]{predict}}. The method assumes a
#' process with
#' homegeneous variance, but the estimated interpolating line is quite robust to
#' heteroskedasticity, anyway, also scale parameter \code{sigma2} is estimated as well as
#' SEs based on it.
#'
#' \code{predict}
#' @param M series of observed integrals
#' @param t_0 time at beginning of the first interval
#' @param t_unit time span of each interval
#' @param t_end time at the end of the last interval. If \code{t_end} is specified,
#' \code{t_unit} is ignored
#' @return A list of class \code{interv_fit}, with the following attributes:
#' \itemize{
#'   \item \code{$data}: series M
#'   \item \code{$knots}: estimated value for the Wiener process in the points between
#'   intervals
#'   \item \code{$sigma2}: estiamated scale parameter
#'   \item \code{$se}: standard errors for the knots values. They depend on sigma2 and on
#'   homoskedasticity assumption
#'   \item \code{$covariances}: covariances of consecutive knots estimates.
#' }
#'
#' @examples M= c(0, 1, -1, 2, 4)
#' mod= interv_fit(M)
#' plot(mod)
#'
#' mod= interv_fit(M, t_0= 5, t_unit= 2)
#' plot(mod)
#' @export
interv_fit= function(M, t_0= 0, t_unit= 1, t_end= NULL) {

	stopifnot(is.numeric(t_0), is.numeric(t_unit), length(t_0) == 1, length(t_unit) == 1)
	if (!is.vector(M) || !is.numeric(M))
		stop('M must be a numeric vector')

	n= length(M)

	if (!is.null(t_end))
		t_unit= (t_end-t_0)/n

	# setting tridiagonal system to find knots values
	B= c(3*M[1], 6*M[-c(1, n)]) + c(6*M[-c(1, n)], 3*M[n])
	diag= c(7, rep(8, n-3), 7)
	coef= rep(2, n-2)
	knots= as.vector(limSolve::Solve.tridiag(coef, diag, coef, B))

	# computations for estimating sigma2
	m_squares= 3*sum(M[c(1, n)]^2) + 12*sum(M[-c(1, n)]^2)
	dot_product= sum(B*knots)
	sigma2= (m_squares - dot_product)/(n-1)

	# Usmani's algorithm for inverting matrix A/2
	m= min(15, n) # if n > 15 inverting a 15x15 matrix is sufficient

	diag= diag/2

	theta= vector('numeric', m)
	theta[1:2]= c(1, 3.5)
	for (i in 2:(m-1))
		theta[i+1]= diag[i]*theta[i]-theta[i-1]

	phi= theta[m:1]
	denominator= theta[m]

	SEs= sqrt(sigma2*theta[-m]*phi[-1]/denominator/2)
	covariances= -sigma2*theta[-c(m-1, m)]*phi[-(1:2)]/denominator/2

	# extending the result for n > 15
	k= n-15
	if (k > 0) {
		SEs= c(SEs[1:7], rep(SEs[7], k), SEs[8:14])
		covariances= c(covariances[1:6], rep(covariances[7], k-1), covariances[8:13])
	}

	out= list(data= M, t_0= t_0, t_unit= t_unit, knots= knots, sigma2= sigma2, se= SEs,
						covariances= covariances)
	class(out)= "interv_fit"
	return(out)
}


#'@export
print.interv_fit= function(x, ...) {
	cat("Object of class interv_fit.\n")
	cat("Number of intervals evaluated:", length(x$data))

	invisible(x)
}


#' Extracting point estimates of the interpolated process
#'
#' For new values of t provided whitin the scope of the model, the corresponding Maximum
#' LIkelihood estimated value of the process is returned
#'
#' @param object an object of class \code{interv_fit}
#' @param t a numeric vector; values must lay between 0 and N, where N corresponds to the
#' lenght of the data vector on which the model was estimated
#' @param se.fit boolean value; if TRUE standard errors of the estimates are returned along
#' the point estimates
#' @param ... further arguments that will not be processed
#' @return Object returned depends on the parameter \code{se.fit}. If set to FALSE, a
#' numerical vector is returned, elsewise a data.frame is returned with a column of estimates
#' and a column with their SEs
#'
#' @examples M= c(0, 1, -1, 2, 4)
#' mod= interv_fit(M, t_0= 0, t_end= 1)
#' predict(mod, c(0, .25, .5, .75, 1))
#' @export
predict.interv_fit= function(object, t, se.fit= FALSE, ...) {

	if (!is.vector(t) || !is.numeric(t))
		stop('A numeric vector is required for t')
	# rescaling t to unit intervals starting from 0
	t= (t-object$t_0)/object$t_unit
	if (any(t < 0 | t > length(object$data)))
		stop('t out of bounds')

	i= ceiling(t) # interval of each point
	i[i == 0]= 1
	t= t+1-i

	coefficients= .parabolas(object)[i, ]
	variables= cbind(1, t, t*t)
	y= rowSums(coefficients*variables)

	if (!se.fit)
		return(y)

	# evaluating SEs requires different methods for extreme intervals
	# for semplicity values from first interval are moved to last one, reverted
	n= length(object$data)
	t[i == 1]= 1-t[i == 1]
	i[i == 1]= n
	extremes= i == n

	# initailizing vector of variances of estimates
	var_y= vector('numeric', length(t))

	# middle intervals
	SE_a= object$se[i[!extremes]-1]
	SE_b= object$se[i[!extremes]]
	cov= object$covariances[i[!extremes]-1]

	t1= t[!extremes]
	coef_a= 1-4*t1+3*t1*t1
	coef_b= -2*t1+3*t1*t1
	var_y[!extremes]= coef_a^2*SE_a^2 + coef_b^2*SE_b^2 + 2*coef_a*coef_b*cov +
													object$sigma2*(t1 -4*t1^2 +6*t1^3 -3*t1^4)

	# extreme intervals
	t2= t[extremes]
	coef_a= 1-3*t2+1.5*t2*t2
	SE_a= object$se[n-1]
	var_y[extremes]= coef_a^2*SE_a^2 + object$sigma2*(t2 - 3*t2^2 + 3*t2^3 - .75*t2^4)

	out= data.frame(prediction= y, SE= sqrt(var_y))
	return(out)
}


#' Plotting a interv_fit object
#'
#' Plots the interpolating curve along the interval data. By default, also plots confidence
#' bands
#'
#' @param x an object of class interv_fit
#' @param from,to limits of the plotting window of time. By default, it plots the window
#' corresponding to first 20 time intervals, or the whole series if it is shorter. A warning
#' will be raised if a part of the series is omitted. By setting these parameters, it is
#' possible to plot a window of any lenght, but windows too wide will make the graph too
#' dense
#' @param points integer, the number of points in which to evaluate the curve
#' @param plot_bands logical; if TRUE, confidence bands are plotted. Default is TRUE
#' @param sigma2 it is possible to fix scale parameter for confidence bands. If
#' \code{sigma2} is fixed, then normal quantiles are used; if is not, t quantiles are used
#' @param col_means color for interval observed integrals
#' @param col_line color for interpolating curve
#' @param col_bands color for the region inside confidence bands. Default uses transparency
#' @param ... any graphical parameter to be used for plotting
#'
#' @examples M= c(0, 1, -1, 2, 4)
#' mod= interv_fit(M)
#' plot(mod)
#' @export
#' @importFrom grDevices rgb
plot.interv_fit= function(x, ..., from= NULL, to= NULL, points= 300,
													plot_bands= TRUE, sigma2= NULL,
													col_means= 1, col_line= 4, col_bands= rgb(0, 0, 1, .25))
	{
	n= length(x$data)

	#TODO: tutti i controlli sui parametrini

	t_0= x$t_0
	t_unit= x$t_unit
	t_end= t_0+n*t_unit
	if (is.null(from))
		from= t_0
	if (is.null(to)) {
		if (t_end-from > 20*t_unit)
			warning('Series is too long, only first 20 intervals were plotted. Manually
							fix plotting window to avoid this limit')
		to= min(from+20*t_unit, t_end)
	}
	t= seq(from, to, length.out= points)

	# computing bands
	ylim= NULL
	if (plot_bands) {
		if (is.null(sigma2)){
			sigma2= x$sigma2
			q= stats::qt(.975, df= n-1)
		}
		else
			q= stats::qnorm(.975)

		p= predict.interv_fit(x, t, se.fit= T)
		y= p$prediction
		upper= y + p$SE*q
		lower= y - p$SE*q
		ylim= c(min(lower), max(upper))
	}

	# setting plotting parameters, the ones in first list are overwritten by those passed to
	# this function, the ones in the second list would rather overwrite them
	pars= list(ylim= ylim, xlab= 't', ylab= 'y', pch= 20)
	pars= utils::modifyList(pars, list(...))
	pars= utils::modifyList(pars, list(x= t, y= y, col= col_line, type= 'l'))
	do.call(graphics::plot, pars)

	if (plot_bands)
		graphics::polygon(c(t, t[points:1]), c(upper, lower[points:1]), border= NA,
											col= col_bands)

	from= floor((from-t_0)/t_unit)
	to= ceiling((to-t_0)/t_unit)
	M= x$data[(from+1):to]
	x= t_0+(from:to)*t_unit
	graphics::segments(x0= x[-length(x)], x1= x[-1], y0= M, y1= M, col= col_means)
	graphics::points(x[-1], M, pch= pars$pch, col= col_means)
}


