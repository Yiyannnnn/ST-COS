
#' Generate spatial coordinates for spots
#' 
#' @param x_len Length of x dimension
#' @param y_len Length of y dimension
#' @param pattern Pattern type: "grid", "hex", or "random"
#' @param spot_distance Distance between spots (for grid/hex patterns)
#' @param n_spot Number of spots (for random pattern)
#' @param seed Random seed
#' @return data.table with x and y coordinates
#' @export
generate_coordinates <- function(x_len, y_len, pattern, spot_distance=NULL, n_spot=NULL, seed=123) {
  set.seed(seed)
  if(pattern == "grid"){
    x = seq(0, x_len, spot_distance)
    y = seq(0, y_len, spot_distance)
    return(data.table(x_coord = rep(x, times=length(y)),
                     y_coord = rep(y, each=length(x))))
  } else if(pattern == "hex"){
    x = seq(0, x_len, spot_distance)
    y = seq(0, y_len, spot_distance*sqrt(3)/2)
    if(length(y)%%2==0){
      x_coord = rep(c(x,x+spot_distance/2), times=length(y)/2)
    } else {
      x_coord = c(rep(c(x,x+spot_distance/2), times=length(y)%/%2), x)
    }
    return(data.table(x_coord = x_coord,
                     y_coord = rep(y, each=length(x))) %>% filter(x_coord<=x_len))
  } else if(pattern == "random"){
    x <- runif(n_spot, min = 0, max = x_len)
    y <- runif(n_spot, min = 0, max = y_len)
    return(data.table(x_coord = rep(x, times=length(y)),
                     y_coord = rep(y, each=length(x)))[sample(seq(1,n_spot^2), n_spot),])
  } else {
    return('Not implemented')
  }
}