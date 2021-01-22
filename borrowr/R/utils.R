
#  borrowr: estimate population average treatment effects with borrowing between data sources.
#  Copyright (C) 2019  Jeffrey A. Verdoliva Boatman
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.

notallzero <- function(x) any(x != 0)

dropzeros <- function(x) {
  if (notallzero(x)) return(x)
  NULL
}

quiet <- function(f) {
  function(...) {
    capture.output(out <- f(...))
    out
  }
}

# none <- Negate(any)

`%!in%` <- Negate(`%in%`)

stdize_matrix <- function(x, m, s) {
  m <- t(array(m, dim = rev(dim(x))))
  s <- t(array(s, dim = rev(dim(x))))
  (x - m) / s
}
