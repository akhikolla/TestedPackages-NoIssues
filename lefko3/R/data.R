#' Demographic Dataset of \emph{Cypripedium candidum} Population, in Horizontal
#' Format
#' 
#' A dataset containing the states and fates of \emph{Cypripedium candidum} 
#' (white lady's slipper orchids), family Orchidaceae, from a population in 
#' Illinois, USA, resulting from monitoring that occurred annually between 2004 
#' and 2009.
#' 
#' @docType data
#' 
#' @usage data(cypdata)
#' 
#' @format A data frame with 77 individuals and 27 variables. Each row 
#' corresponds to an unique individual, and each variable from \code{size.04} 
#' on refers to the state of the individual in a particular year.
#' 
#' \describe{
#'   \item{plantid}{A numberic variable giving a unique number to each 
#'   individual.}
#'   \item{patch}{A variable refering to patch within the population.}
#'   \item{censor}{A variable coding for whether the data point is valid. An
#'   entry of 1 means that it is so.}
#'   \item{Inf2.04}{Number of double inflorescences in 2004.}
#'   \item{Inf.04}{Number of inflorescences in 2004.}
#'   \item{Veg.04}{Number of stems without inflorescences in 2004.}
#'   \item{Pod.04}{Number of fruits in 2004.}
#'   \item{Inf2.05}{Number of double inflorescences in 2005.}
#'   \item{Inf.05}{Number of inflorescences in 2005.}
#'   \item{Veg.05}{Number of stems without inflorescences in 2005.}
#'   \item{Pod.05}{Number of fruits in 2005.}
#'   \item{Inf2.06}{Number of double inflorescences in 2006.}
#'   \item{Inf.06}{Number of inflorescences in 2006.}
#'   \item{Veg.06}{Number of stems without inflorescences in 2006.}
#'   \item{Pod.06}{Number of fruits in 2006.}
#'   \item{Inf2.07}{Number of double inflorescences in 2007.}
#'   \item{Inf.07}{Number of inflorescences in 2007.}
#'   \item{Veg.07}{Number of stems without inflorescences in 2007.}
#'   \item{Pod.07}{Number of fruits in 2007.}
#'   \item{Inf2.08}{Number of double inflorescences in 2008.}
#'   \item{Inf.08}{Number of inflorescences in 2008.}
#'   \item{Veg.08}{Number of stems without inflorescences in 2008.}
#'   \item{Pod.08}{Number of fruits in 2008.}
#'   \item{Inf2.09}{Number of double inflorescences in 2009.}
#'   \item{Inf.09}{Number of inflorescences in 2009.}
#'   \item{Veg.09}{Number of stems without inflorescences in 2009.}
#'   \item{Pod.09}{Number of fruits in 2009.}
#' }
#' 
#' @source Shefferson, R.P., R. Mizuta, and M.J. Hutchings. 2017. Predicting
#' evolution in response to climate change: the example of sprouting probability
#' in three dormancy-prone orchid species. \emph{Royal Society Open Science} 
#' 4(1):160647.
#' 
#' @examples 
#' data(cypdata)
#' 
#' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
#' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", 
#'                  "Sm", "Md", "Lg", "XLg")
#' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 2.5, 7)
#' 
#' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector, 
#'                           repstatus = repvector, obsstatus = obsvector, 
#'                           matstatus = matvector, propstatus = propvector, 
#'                           immstatus = immvector, indataset = indataset, 
#'                           binhalfwidth = binvec)
#' 
#' cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004, 
#'                           patchidcol = "patch", individcol = "plantid", 
#'                           blocksize = 4, sizeacol = "Inf2.04", sizebcol = "Inf.04", 
#'                           sizeccol = "Veg.04", repstracol = "Inf.04", 
#'                           repstrbcol = "Inf2.04", fecacol = "Pod.04", 
#'                           stageassign = cypframe_raw, stagesize = "sizeadded", 
#'                           NAas0 = TRUE, NRasRep = TRUE)
#' 
#' rep_cyp_raw <- matrix(0, 11, 11)
#' rep_cyp_raw[1:2,7:11] <- 0.5
#' 
#' cypover2r <- overwrite(stage3 = c("SD", "P1", "P2", "P3", "SL", "SL", "D", 
#'                        "XSm", "Sm"), stage2 = c("SD", "SD", "P1", "P2", "P3", 
#'                        "SL", "SL", "SL", "SL"), eststage3 = c(NA, NA, NA, NA, 
#'                        NA, NA, "D", "XSm", "Sm"), eststage2 = c(NA, NA, NA, NA, 
#'                        NA, NA, "XSm", "XSm", "XSm"), givenrate = c(0.1, 0.2, 
#'                        0.2, 0.2, 0.25, 0.4, NA, NA, NA), type = c("S", "S", "S",
#'                        "S", "S", "S", "S", "S", "S"))
#' 
#' cypmatrix2r <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, year = "all", 
#'                        patch = "all", stages = c("stage3", "stage2"),
#'                        size = c("size3added", "size2added"),
#'                        repmatrix = rep_cyp_raw, overwrite = cypover2r,
#'                        yearcol = "year2", patchcol = "patchid",
#'                        indivcol = "individ")
#' cypmatrix2r$A[[1]]
#' 
"cypdata"

#' Demographic Dataset of \emph{Cypripedium candidum} Population, in Vertical
#' Format
#' 
#' A dataset containing the states and fates of \emph{Cypripedium candidum} 
#' (white lady's slipper orchids), family Orchidaceae, from a population in 
#' Illinois, USA, resulting from monitoring that occurred annually between 2004 
#' and 2009. Same dataset as \code{cypdata}, but arranged in an ahistorical
#' vertical format.
#' 
#' @docType data
#' 
#' @usage data(cypvert)
#' 
#' @format A data frame with 77 individuals, 331 rows, and 12 variables. Each
#' row corresponds to a specific two-year transition for a specific individual.
#' Variable codes are similar to those for \code{cypdata}, but use \code{.2} to
#' identify time \emph{t} and \code{.3} to identify time \emph{t}+1.
#' 
#' \describe{
#'   \item{plantid}{A numberic variable giving a unique number to each 
#'   individual.}
#'   \item{patch}{A variable refering to patch within the population.}
#'   \item{censor}{A variable coding for whether the data point is valid. An
#'   entry of 1 means that it is so.}
#'   \item{year2}{Year in time \emph{t}.}
#'   \item{Inf2.2}{Number of double inflorescences in time \emph{t}.}
#'   \item{Inf.2}{Number of inflorescences in time \emph{t}.}
#'   \item{Veg.2}{Number of stems without inflorescences in time \emph{t}.}
#'   \item{Pod.2}{Number of fruits in time \emph{t}.}
#'   \item{Inf2.3}{Number of double inflorescences in time \emph{t}+1.}
#'   \item{Inf.3}{Number of inflorescences in time \emph{t}+1.}
#'   \item{Veg.3}{Number of stems without inflorescences in time \emph{t}+1.}
#'   \item{Pod.3}{Number of fruits in time \emph{t}+1.}
#' }
#' 
#' @source Shefferson, R.P., R. Mizuta, and M.J. Hutchings. 2017. Predicting
#' evolution in response to climate change: the example of sprouting probability
#' in three dormancy-prone orchid species. \emph{Royal Society Open Science} 
#' 4(1):160647.
#' 
#' @examples 
#' data(cypvert)
#' 
#' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
#' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md", "Lg", "XLg")
#' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 2.5, 7)
#' 
#' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector, 
#'                           repstatus = repvector, obsstatus = obsvector,
#'                           matstatus = matvector, propstatus = propvector,
#'                           immstatus = immvector, indataset = indataset,
#'                           binhalfwidth = binvec)
#' 
#' cypframe_raw
#' 
#' cypraw_v2 <- historicalize3(data = cypvert, patchidcol = "patch", individcol = "plantid",
#'                             year2col = "year2", sizea2col = "Inf2.2", sizea3col = "Inf2.3",
#'                             sizeb2col = "Inf.2", sizeb3col = "Inf.3", sizec2col = "Veg.2",
#'                             sizec3col = "Veg.3", repstra2col = "Inf2.2", repstra3col = "Inf2.3",
#'                             repstrb2col = "Inf.2", repstrb3col = "Inf.3", feca2col = "Pod.2",
#'                             feca3col = "Pod.3", repstrrel = 2, stageassign = cypframe_raw,
#'                             stagesize = "sizeadded", censorcol = "censor", censor = FALSE,
#'                             NAas0 = TRUE, NRasRep = TRUE, reduce = TRUE)
#' summary(cypraw_v2)
#' 
"cypvert"

#' Demographic Dataset of \emph{Lathyrus vernus} Population
#' 
#' A dataset containing the states and fates of \emph{Lathyrus vernus} (spring
#' vetch), family Fabaceae, from a population in Sweden monitored annually
#' from 1988 to 1991 in six study plots.
#' 
#' @docType data
#' 
#' @usage data(lathyrus)
#' 
#' @format A data frame with 1119 individuals and 34 variables. Each row
#' corresponds to a unique individual, and each variable from \code{Volume88}
#' on refers to the state of the individual in a given year.
#' 
#' \describe{
#'   \item{SUBPLOT}{A variable refering to patch within the population.}
#'   \item{GENET}{A numberic variable giving a unique number to each 
#'   individual.}
#'   \item{Volume88}{Aboveground volume in cubic mm in 1988.}
#'   \item{lnVol88}{Natural logarithm of \code{Volume88}.}
#'   \item{FCODE88}{Equals 1 if flowering and 0 if not flowering in 1988.}
#'   \item{Flow88}{Number of flowers in 1988.}
#'   \item{Intactseed88}{Number of intact mature seeds produced in 1988.
#'   Not always an integer, as in some cases seed number was estimated via 
#'   linear modeling.}
#'   \item{Dead1988}{Marked as 1 if known to be dead in 1988.}
#'   \item{Dormant1988}{Marked as 1 if known to be alive but vegetatively 
#'   dormant in 1988.}
#'   \item{Missing1988}{Marked as 1 if not found in 1988.}
#'   \item{Seedling1988}{Marked as 1, 2, or 3 if observed as a seedling in year
#'   \emph{t}. Numbers refer to certainty of assignment: 1 = certain that plant
#'   is a seedling in 1988, 2 = likely that plant is a seedling in 1988,
#'   3 = probable that plant is a seedling in 1988.}
#'   \item{Volume89}{Aboveground volume in cubic mm in 1989.}
#'   \item{lnVol89}{Natural logarithm of \code{Volume89}.}
#'   \item{FCODE89}{Equals 1 if flowering and 0 if not flowering in 1989.}
#'   \item{Flow89}{Number of flowers in 1989.}
#'   \item{Intactseed89}{NZumber of intact mature seeds produced in 1989.
#'   Not always an integer, as in some cases seed number was estimated via
#'   linear modeling.}
#'   \item{Dead1989}{Marked as 1 if known to be dead in 1989.}
#'   \item{Dormant1989}{Marked as 1 if known to be alive but vegetatively 
#'   dormant in 1989.}
#'   \item{Missing1989}{Marked as 1 if not found in 1989.}
#'   \item{Seedling1989}{Marked as 1, 2, or 3 if observed as a seedling in
#'   year \emph{t}. Numbers refer to certainty of assignment: 1 = certain 
#'   that plant is a seedling in 1989, 2 = likely that plant is a seedling 
#'   in 1989, 3 = probable that plant is a seedling in 1989.}
#'   \item{Volume90}{Aboveground volume in mm<sup>3</sup> in 1990.}
#'   \item{lnVol90}{Natural logarithm of \code{Volume90}.}
#'   \item{FCODE90}{Equals 1 if flowering and 0 if not flowering in 1990.}
#'   \item{Flow90}{Number of flowers in 1990.}
#'   \item{Intactseed90}{NZumber of intact mature seeds produced in 1990.
#'   Not always an integer, as in some cases seed number was estimated via 
#'   linear modeling.}
#'   \item{Dead1990}{Marked as 1 if known to be dead in 1990.}
#'   \item{Dormant1990}{Marked as 1 if known to be alive but vegetatively 
#'   dormant in 1990.}
#'   \item{Missing1990}{Marked as 1 if not found in 1990.}
#'   \item{Seedling1990}{Marked as 1, 2, or 3 if observed as a seedling in
#'   year \emph{t}. Numbers refer to certainty of assignment: 1 = certain 
#'   that plant is a seedling in 1990, 2 = likely that plant is a seedling
#'   in 1990, 3 = probable that plant is a seedling in 1990.}
#'   \item{Volume91}{Aboveground volume in mm<sup>3</sup> in 1991.}
#'   \item{lnVol91}{Natural logarithm of \code{Volume91}.}
#'   \item{FCODE91}{Equals 1 if flowering and 0 if not flowering in 1991.}
#'   \item{Flow91}{Number of flowers in 1991.}
#'   \item{Intactseed91}{NZumber of intact mature seeds produced in 1991.
#'   Not always an integer, as in some cases seed number was estimated via
#'   linear modeling.}
#'   \item{Dead1991}{Marked as 1 if known to be dead in 1991.}
#'   \item{Dormant1991}{Marked as 1 if known to be alive but vegetatively 
#'   dormant in 1991.}
#'   \item{Missing1991}{Marked as 1 if not found in 1991.}
#'   \item{Seedling1991}{Marked as 1, 2, or 3 if observed as a seedling 
#'   in year \emph{t}. Numbers refer to certainty of assignment: 
#'   1 = certain that plant is a seedling in 1991, 2 = likely that plant 
#'   is a seedling in 1991, 3 = probable that plant is a seedling in 
#'   1991.}
#' }
#' 
#' @source Ehrlen, J. 2000. The dynamics of plant populations: does the 
#' history of individuals matter? \emph{Ecology} 81(6):1675-1684.
#' 
#' @examples
#' data(lathyrus)
#' 
#' sizevector <- c(0, 4.6, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9)
#' stagevector <- c("Sd", "Sdl", "Dorm", "Sz1nr", "Sz2nr", "Sz3nr", "Sz4nr", "Sz5nr",
#'                  "Sz6nr", "Sz7nr", "Sz8nr", "Sz9nr", "Sz1r", "Sz2r", "Sz3r", "Sz4r",
#'                  "Sz5r", "Sz6r", "Sz7r", "Sz8r", "Sz9r")
#' repvector <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 4.6, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
#'             0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
#' 
#' lathframeln <- sf_create(sizes = sizevector, stagenames = stagevector, repstatus = repvector,
#'                          obsstatus = obsvector, matstatus = matvector, immstatus = immvector,
#'                          indataset = indataset, binhalfwidth = binvec, propstatus = propvector)
#' 
#' lathvertln <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, patchidcol = "SUBPLOT",
#'                            individcol = "GENET", blocksize = 9, juvcol = "Seedling1988",
#'                            sizeacol = "lnVol88", repstracol = "FCODE88",
#'                            fecacol = "Intactseed88", deadacol = "Dead1988",
#'                            nonobsacol = "Dormant1988", stageassign = lathframeln,
#'                            stagesize = "sizea", censorcol = "Missing1988",
#'                            censorkeep = NA, NAas0 = TRUE, censor = TRUE)
#' 
#' summary(lathvertln)
#' 
"lathyrus"
