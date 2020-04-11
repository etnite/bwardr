#' Create a Color Bar
#' 
#' @param pal Input color palette
#' @param min Minimum value for colorbar axis
#' @param max Maximum value for colorbar axis
#' @param width Value specifying width of the color bar rectangle - must be
#'   between 0 and 10
#' @param nticks Integer specifying number of tick marks to label color bar axis
#'   with
#' @param ticks Numeric vector specifying position of tick marks, for instance as
#'   produced by seq()
#' @param title Title for the color bar plot
#' @details All credit goes to Stackoverflow user John Colby:
#' https://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette
#' @export
color_bar <- function(pal, min, max = -min, width = 1, 
                      nticks = 11, ticks = seq(min, max, len = nticks), 
                      title = "") {
  
  scale = (length(pal) - 1)/(max - min)
  
  grDevices::dev.new(width = 1.75, height = 5)
  graphics::plot(c(0, 10), c(min, max), type = "n", bty = "n", 
                 xaxt = "n", xlab = "", 
                 yaxt = "n", ylab = "", main = title)
  graphics::axis(2, ticks, las = 1)
  for (i in 1:(length(pal) - 1)) {
    y = (i - 1)/scale + min
    graphics::rect(0, y, width, y + 1/scale, col = pal[i], border = NA)
  }
}