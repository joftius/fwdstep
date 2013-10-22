
library(xtable)

tex_table = function(texfile, matrix, caption) {
  f = file(texfile)
  textable = xtable(matrix, caption = caption)
  writeLines(print(textable), f)
  close(f)
}
