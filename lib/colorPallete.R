# DEV - WIP
# List of color function.
#

col.base.content <- "black"

# red
col.enhanced.risk <- "red"

# green
col.suppressed.risk <- "#5fd38d"

# A function to add transparency to #?????? color
add.alpha.channel <- function(hex.color, alpha) {
  # alpha is between 0 and 1 (like percentage)

  rgb.channel <- col2rgb(hex.color)

  rgb.channel <- rbind(rgb.channel, alpha * 255)

  hex.color <- apply(rgb.channel, 2, function(i) do.call(rgb, as.list(i / 255)))

  return(hex.color)
}
