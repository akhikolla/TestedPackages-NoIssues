## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval = FALSE)

## -----------------------------------------------------------------------------
#  library(fplot)
#  setFplot_page(page = "a4", margins = 1)

## ---- eval = FALSE------------------------------------------------------------
#  pdf_fit("first_export.pdf")
#  plot(1, 1, type = "n", ann = FALSE)
#  text(1, 1, "This text will be displayed in 10pt.")
#  fit.off()

## ---- eval = FALSE------------------------------------------------------------
#  pdf_fit("second_export.pdf", pt = 12)
#  plot(1, 1, type = "n", ann = FALSE)
#  text(1, 1, "This text will be displayed in 12pt.")
#  fit.off()

## ---- eval = FALSE------------------------------------------------------------
#  pdf_fit("third_export.pdf", pt = 12, height = 1)
#  plot(1, 1, type = "n", ann = FALSE)
#  text(1, 1, "This text will be displayed in 12pt.")
#  fit.off()

## ---- eval = FALSE------------------------------------------------------------
#  # You can also set the point size globally
#  setFplot_page(pt = 12)
#  pdf_fit("fourth_export.pdf", sideways = TRUE)
#  plot(1, 1, type = "n", ann = FALSE)
#  text(1, 1, "This text will be displayed in 12pt.")
#  fit.off()

