#' @name pcp
#' @aliases pcp
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com}
#' @title Precipitation Concentration Period
#' @description Calculates the Precipitation Concentration Period (PCP) on a 
#' daily or monthly precipitation serie.
#' @usage pcp(object)
#' @param object a daily or monthly precipitation serie.
#' @return A data.frame containing the following variables:
#' \itemize{
#' \item \code{year} is the year.
#' \item \code{pcp} is the precipitation concentration period, in degree, corresponding to a year.
#' Results correspond to a month like below when using the `azimuth`  default values:
#' 0 = January, 30 = February, 60 = March, \dots, 300 = November, and 330 = December.
#' }
#' @seealso 
#' \code{\link{pplot.pcp}}
#' \code{\link{read.data}}
#' \code{\link{as.daily}}
#' \code{\link{as.monthly}}
#' @examples 
#' ##
#' # Loading the monthly precipitation serie.
#' data(monthly)
#' 
#' ## 
#' # Performing the Precipitation Concentration Degree analysis
#' pcd(monthly)
#' @references Zhang L.J., Qian Y.F. (2003) Annual distribution features of precipitation in China and their interannual variations. J Acta Meteorological Sinica 17:146-163
#' @keywords precipitation concentration degree PCD
#' @export
# pcp <- function(object) {
	
# 	if (is.element("precintcon.daily", class(object)))
# 		object <- as.precintcon.monthly(object)
	
# 	if (!is.element("precintcon.monthly", class(object)))
# 		stop("Invalid data. Please, check your input object.")
	
# 	azimuth <- 360 * object$month / 12
	
# 	rx <- aggregate(object$precipitation * sin(azimuth), by = list(object$year), FUN = sum)[2]
	
# 	ry <- aggregate(object$precipitation * cos(azimuth), by = list(object$year), FUN = sum)[2]
	
# 	pcp = atan(rx / ry) / 0.0174532925
	
# 	r <- data.frame(year = unique(object$year), pcp = pcp)
	
# 	colnames(r) <- c("year", "pcp")
	
#   return(r)
# }

pcp <- function(object) {

  if (inherits(object, "precintcon.daily")) {
    object <- as.precintcon.monthly(object)
  }

  if (!inherits(object, "precintcon.monthly")) {
    stop("Invalid data. Please, check your input object.")
  }

  # 月份转换为弧度：1月为0，12月为11π/6
  theta <- 2 * pi * (object$month - 1) / 12

  # 按年合计向量分量及降水总量；缺失值保留
  components <- aggregate(
    data.frame(
      rx = object$precipitation * sin(theta),
      ry = object$precipitation * cos(theta),
      total = object$precipitation
    ),
    by = list(year = object$year),
    FUN = sum
  )

  # 合成向量长度
  resultant <- sqrt(components$rx^2 + components$ry^2)

  # 无降水、含缺失值或合成向量接近零时，集中期不定义
  valid <- is.finite(components$total) &
    is.finite(resultant) &
    components$total > 0 &
    resultant > 1e-12 * components$total

  pcp <- rep(NA_real_, nrow(components))

  # 按论文坐标约定：rx 为正弦分量，ry 为余弦分量
  pcp[valid] <- (
    atan2(components$rx[valid], components$ry[valid]) *
      180 / pi
  ) %% 360

  data.frame(
    year = components$year,
    pcp = pcp
  )
}