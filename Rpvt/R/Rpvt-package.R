#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
#' @importFrom Rcpp evalCpp
#' @importFrom Rdpack reprompt
## usethis namespace: end
NULL

#' Create a matrix of PVT properties for dry- and wet-gas samples
#'
#' The pvt_gas() generates a table of gas PVT properties at reservoir temperature and pressures from the atmospheric condition up to the initial reservoir pressure. The estimated properties are compressibility factor, formation volume factor, density, compressibility, viscosity, and pseudo-pressure.
#'
#' @param input_unit input unit system for parameters, a character string either 'SI' or 'Field'
#' @param output_unit output unit system for properties, a character string either 'SI' or 'Field'
#' @param fluid fluid type, a character string either 'dry_gas' or 'wet_gas'
#' @param pvt_model PVT model, the character string 'DAK'
#' @param visc_model viscosity model, the character string 'Sutton'
#' @param t temperature, a numeric value either in 'C' or 'F' depending on the 'input_unit'
#' @param p pressure, a numeric value either in 'kPag' or 'Psig' depending on the 'input_unit'
#' @param gas_spgr gas specific gravity (Air = 1.0)
#' @param nhc_composition a vector of mole fractions for nitrogen, hydrogen sulfide, and carbon dioxide, respectively
#' @param cgr condensate to gas ratio, a numeric value in 'm3/m3' or 'STB/MMSCF' depending on the 'input_unit'
#' @param cond_api condensate API
#' @param warning a charater string either 'yes' or 'no'
#' @useDynLib Rpvt, .registration = TRUE
#'
#' @export
#'
#' @examples
#' pvt_gas_results_1 <- pvt_gas(input_unit = "Field", output_unit = "Field",
#' fluid = "dry_gas", pvt_model = "DAK", visc_model = "Sutton",
#' t = 400, p = 20000, gas_spgr = 0.65, nhc_composition = c(0.05,0.02,0.04),
#' cgr = 0.0, cond_api = NULL, warning = "yes")
#'
#' head(pvt_gas_results_1)
#'
#' pvt_gas_results_2 <- pvt_gas(input_unit = "Field", output_unit = "Field",
#' fluid = "wet_gas", pvt_model = "DAK", visc_model = "Sutton",
#' t = 300, p = 20000, gas_spgr = 0.75, nhc_composition = c(0.02,0.05,0.08),
#' cgr = 10.0, cond_api = 42.4, warning = "yes")
#'
#' head(pvt_gas_results_2)
#'
#' @references
#'
#' \insertRef{Sutton2007}{Rpvt}
#'
#' \insertRef{Wichert1972}{Rpvt}



pvt_gas <- function(input_unit = "Field", output_unit = "Field", fluid = "wet_gas", pvt_model = "DAK", visc_model = "Sutton", t = 300, p = 5000, gas_spgr = 0.69, nhc_composition = c(0.0, 0.0, 0.0), cgr = 3.0, cond_api = 42.0, warning = "yes") {

   kpa_to_psi <- 0.14503773800722
   if (!is.character(input_unit)) stop("'input_unit' must be a character string either 'SI' or 'Field'.")
   if (!is.character(output_unit)) stop("'output_unit' must be a character string either 'SI' or 'Field'.")
   if (!is.character(fluid)) stop("'fluid' must be a character string either 'dry_gas' or 'wet_gas'.")
   if (!is.character(pvt_model)) stop("'pvt_model' must be the character string 'DAK'.")
   if (!is.character(visc_model)) stop("'visc_model' must be the character string 'Sutton'.")
   if (!(warning %in% c("no","yes"))) stop("'warning' must be a character string either 'yes' or 'no'.")

   if (fluid == "dry_gas") {
      if (pvt_model != "DAK") stop("Change the gas PVT model to 'DAK'.")
      if (visc_model != "Sutton") stop("Change the gas viscosity model to 'Sutton'.")
      if (!is.numeric(t)) stop("'t' must be a numeric.")
      if (!is.numeric(p)) stop("'p' must be a numeric.")
      if (!is.numeric(gas_spgr)) stop("'gas_spgr' must be a numeric.")
      if (!is.numeric(nhc_composition)) stop("'nhc_composition' must be a numeric vector.")
      if (cgr != 0) stop("'cgr' must be zero for dry gas samples.")
      if (input_unit == "SI") {
         p_ <- p * kpa_to_psi
         t_ <- t * 1.8 + 32
      } else {
         p_ <- p
         t_ <- t
      }
      if (warning == "yes") {
         if (nhc_composition[1] > 0.217) warning('N2 composition is greater than dry gas PVT correlation upper limit: 0.217 mol fraction')
         if (nhc_composition[2] > 0.10) warning('H2S composition is greater than dry gas PVT correlation upper limit: 0.10 mol fraction')
         if (nhc_composition[3] > 0.558) warning('CO2 composition is greater than dry gas PVT correlation upper limit: 0.558 mol fraction')
         if ((gas_spgr < 0.554) | (gas_spgr > 1.862)) warning('gas specific gravity is outside of the dry gas PVT correlation range: [0.554, 1.861]')
         if (p_ > 10000) warning('pressure is greater than dry gas PVT correlation upper limit: 10000 psig')
         if ((t_ < 32) | (t_ > 460)) warning('temperature is outside of the dry gas PVT correlation range: [32 F, 460 F]')
         if (nhc_composition[1] > 0.052) warning('N2 composition is greater than gas viscosity correlation upper limit: 0.052 mol fraction')
         if (nhc_composition[2] > 0.017) warning('H2S composition is greater than gas viscosity correlation upper limit: 0.017 mol fraction')
         if (nhc_composition[3] > 0.089) warning('CO2 composition is greater than gas viscosity correlation upper limit: 0.089 mol fraction')
         if ((gas_spgr < 0.554) | (gas_spgr > 1.861)) warning('gas specific gravity is outside of the gas viscosity correlation range: [0.554, 1.861]')
         if (p_ > 20000) warning('pressure is greater than gas viscosity correlation upper limit: 20000 psig')
         if (t_ > 652) warning('temperature is greater than gas viscosity correlation upper limit: 652 F')
      }
      gas_lst <- list(input_unit = input_unit, output_unit = output_unit, fluid = fluid, pvt_model = pvt_model, visc_model = visc_model, t = t, p = p, gas_spgr = gas_spgr, nhc_composition = nhc_composition, cgr = 0.0, cond_api = 40.0, warning = warning)
      mat_results <- PVT_GAS_PROPERTIES(gas_lst)
      return(mat_results)
   }
   if (fluid == "wet_gas") {
      if (pvt_model != "DAK") stop("Change gas PVT model to 'DAK'.")
      if (visc_model != "Sutton") stop("Change gas viscosity model to 'Sutton'.")
      if (!is.numeric(t)) stop("'t' must be a numeric.")
      if (!is.numeric(p)) stop("'p' must be a numeric.")
      if (!is.numeric(gas_spgr)) stop("'gas_spgr' must be a numeric.")
      if (!is.numeric(nhc_composition)) stop("'nhc_composition' must be a numeric vector.")
      if (cgr <= 0) stop("'cgr' must be greater than zero for wet gas samples.")
      if (cgr > 0) {
         if (is.null(cond_api)) stop("'cond_api' value is required for wet_gas samples.")
      }
      if (input_unit == "SI") {
         p_ <- p * kpa_to_psi
         t_ <- t * 1.8 + 32
      } else {
         p_ <- p
         t_ <- t
      }
      if (warning == "yes") {
         if (nhc_composition[1] > 0.333) warning('N2 composition is greater than wet gas PVT correlation upper limit: 0.333 mol fraction')
         if (nhc_composition[2] > 0.90) warning('H2S composition is greater than wet gas PVT correlation upper limit: 0.90 mol fraction')
         if (nhc_composition[3] > 0.899) warning('CO2 composition is greater than wet gas PVT correlation upper limit: 0.899 mol fraction')
         if ((gas_spgr < 0.554) | (gas_spgr > 2.667)) warning('gas specific gravity is outside of the wet gas correlation range: [0.554, 2.667]')
         if (p_ > 17000) warning('pressure is greater than wet gas PVT correlation upper limit: 17000 psig')
         if ((t_ < 32) | (t_ > 460)) warning('temperature is outside of the wet gas PVT correlation range: [32 F, 460 F]')
         if (nhc_composition[1] > 0.052) warning('N2 composition is greater than gas viscosity correlation upper limit: 0.052 mol fraction')
         if (nhc_composition[2] > 0.017) warning('H2S composition is greater than gas viscosity correlation upper limit: 0.017 mol fraction')
         if (nhc_composition[3] > 0.089) warning('CO2 composition is greater than gas viscosity correlation upper limit: 0.089 mol fraction')
         if ((gas_spgr < 0.554) | (gas_spgr > 1.861)) warning('gas specific gravity is outside of the gas viscosity correlation range: [0.554, 1.861]')
         if (p_ > 20000) warning('pressure is greater than gas viscosity correlation upper limit: 20000 psig')
         if (t_ > 652) warning('temperature is greater than gas viscosity correlation upper limit: 652 F')
      }
      gas_lst <- list(input_unit = input_unit, output_unit = output_unit, fluid = fluid, pvt_model = pvt_model, visc_model = visc_model,
                      t = t, p = p, gas_spgr = gas_spgr, nhc_composition = nhc_composition, cgr = cgr, cond_api = cond_api, warning = warning)
      mat_results <- PVT_GAS_PROPERTIES(gas_lst)
      return(mat_results)
   }
}




#' Create a matrix of PVT properties for black oil samples
#'
#' The pvt_oil() generates a table of oil and gas PVT properties at reservoir temperature and pressures from the atmospheric condition up to the initial reservoir pressure. The estimated oil properties are solution gas-oil ratio, formation volume factor, density, compressibility, and viscosity. Estimated PVT properties for the associated gas are compressibility factor, formation volume factor, density, compressibility, viscosity, and pseudo-pressure.
#'
#' @param input_unit input unit system for parameters, a character string either 'SI' or 'Field'
#' @param output_unit output unit system for properties, a character string either 'SI' or 'Field'
#' @param fluid fluid type, the character string 'black_oil'
#' @param pvt_model PVT model, a character string. 'Standing', 'Vasquez_Beggs', 'Farshad_Petrosky', 'Al_Marhoun', and 'Glaso' models are currently available
#' @param visc_model viscosity model, a character string. 'Beggs_Robinson', and 'Al_Marhoun' models are currently available
#' @param t temperature, a numeric value either in 'C' or 'F' depending on the 'input_unit'
#' @param p pressure, a numeric value either in 'kPag' or 'Psig' depending on the 'input_unit'
#' @param oil_api API gravity of oil
#' @param gas_spgr gas specific gravity (Air = 1.0)
#' @param nhc_composition a vector of mole fractions for nitrogen, hydrogen sulfide, and carbon dioxide, respectively
#' @param rsi initial solution gas oil ratio in  'm3/m3' or 'SCF/STB' depending on the 'input_unit'. It is either NULL or a numeric value. If 'rsi' is NULL, then a numeric value must be assigned to 'pb'
#' @param pb bubble point pressure, a numeric value either in 'kPag' or 'Psig' depending on the 'input_unit'. it is either NULL or a numeric value. If 'pb' is NULL, then a numeric value must be assigned to 'rsi'
#' @param warning a charater string either 'yes' or 'no'
#' @useDynLib Rpvt, .registration = TRUE
#'
#' @export
#'
#' @examples
#' pvt_oil_results_1 <- pvt_oil(input_unit = "Field", output_unit = "Field", fluid = "black_oil",
#' pvt_model = "Standing", visc_model = "Beggs_Robinson",
#' t = 200, p = 3000, oil_api = 35, gas_spgr = 0.8,
#' nhc_composition = c(0.05,0.02,0.04), rsi = 650, pb = NULL, warning = "no")
#'
#' head(pvt_oil_results_1)
#'
#' pvt_oil_results_2 <- pvt_oil(input_unit = "SI", output_unit = "Field", fluid = "black_oil",
#' pvt_model = "Vasquez_Beggs", visc_model = "Al_Marhoun",
#' t = 80, p = 20000, oil_api = 34.4, gas_spgr = 0.65,
#' nhc_composition = c(0.0,0.0,0.0), rsi = NULL, pb = 11350, warning = "yes")
#'
#' head(pvt_oil_results_2)
#'
#' @references
#'
#' \insertRef{Standing1947}{Rpvt}
#'
#' \insertRef{Vasquez1980}{Rpvt}
#'
#' \insertRef{Glaso1980}{Rpvt}
#'
#' \insertRef{PetroskyJr.1998}{Rpvt}
#'
#' \insertRef{Al-Marhoun1988}{Rpvt}
#'
#' \insertRef{Beggs1975}{Rpvt}
#'
#' \insertRef{Al-Marhoun2004}{Rpvt}
#'
#' \insertRef{Spivey2007}{Rpvt}
#'
#' \insertRef{Sutton2007}{Rpvt}



pvt_oil <- function(input_unit = "SI", output_unit = "SI", fluid = "black_oil", pvt_model = "Standing", visc_model = "Beggs_Robinson", t = 85.4, p = 35000, oil_api = 38, gas_spgr = 0.67, nhc_composition = c(00, 0.0, 0.0), rsi = NULL, pb = 29500 , warning = "yes") {

   kpa_to_psi <- 0.14503773800722
   if (input_unit == "SI") {
      p_ <- p * kpa_to_psi
      t_ <- t * 1.8 + 32
   } else {
      p_ <- p
      t_ <- t
   }
   if (!is.character(input_unit)) stop("'input_unit' must be a character string either 'SI' or 'Field'.")
   if (!is.character(output_unit)) stop("'output_unit' must be a character string either 'SI' or 'Field'.")
   if (!is.character(fluid)) stop("'fluid' must be the character string 'black_oil'.")
   if (!is.character(pvt_model)) stop("'pvt_model' must be a character string.")
   if (!is.character(visc_model)) stop("'visc_model' must be the character string.")
   if (!(warning %in% c("no","yes"))) stop("'warning' must be a character string either 'yes' or 'no'.")
   if (!(pvt_model %in% c("Standing","Vasquez_Beggs","Farshad_Petrosky","Al_Marhoun"))) stop("Selected oil PVT model must be one of the 'Standing', 'Vasquez_Beggs', 'Farshad_Petrosky', or 'Al_Marhoun' correlations.")
   if (!(visc_model %in% c("Beggs_Robinson","Al_Marhoun"))) stop("Selected oil viscosity model must be either 'Beggs_Robinson' or 'Al_Marhoun'.")
   if (!is.numeric(t)) stop("'t' must be a numeric.")
   if (!is.numeric(p)) stop("'p' must be a numeric.")
   if (!is.numeric(oil_api)) stop("'oil_api' must be a numeric.")
   if (!is.numeric(gas_spgr)) stop("'gas_spgr' must be a numeric.")
   if (!is.numeric(nhc_composition)) stop("'nhc_composition' must be a numeric vector of mole fractions of nitrogen, hydrogen sulfide, and carbon dioxide, respectively in the surface gas.")
   if (is.null(rsi) & is.null(pb)) stop("Either 'rsi' or 'pb' must be provided as a numeric.")
   # Priority is given to Pb
   if (is.null(rsi)) {
      if (pb <= 0) stop("'pb' must be greater than zero.")
      if (warning == "yes") {
         if (nhc_composition[1] > 0.217) warning('N2 composition is greater than dry gas PVT correlation upper limit: 0.217 mol fraction')
         if (nhc_composition[2] > 0.10) warning('H2S composition is greater than dry gas PVT correlation upper limit: 0.10 mol fraction')
         if (nhc_composition[3] > 0.558) warning('CO2 composition is greater than dry gas PVT correlation upper limit: 0.558 mol fraction')
         if ((oil_api < 16) | (oil_api > 58)) warning('oil API gravity is outside of the black_oil PVT correlation range: [16 API ,58 API]')
         if ((gas_spgr < 0.554) | (gas_spgr > 1.862)) warning('gas specific gravity is outside of the dry gas PVT correlation range: [0.554, 1.861]')
         if (p_ > 10000) warning('pressure is greater than dry gas PVT correlation upper limit: 10000 psig')
         if ((t_ < 32) | (t_ > 460)) warning('temperature is outside of the dry gas PVT correlation range: [32 F, 460 F]')
         if ((gas_spgr < 0.56) | (gas_spgr > 1.18)) warning('gas specific gravity is outside of the black_oil PVT correlation range: [0.56, 1.18]')
         if (p_ > 5250) warning("pressure is greater than black oil PVT correlation upper limit: 5250 psig",call. = FALSE)
         if ((t_ < 70) | (t_ > 295)) warning('temperature is outside of the black_oil PVT correlation range: [70 F, 295 F]')
         if (nhc_composition[1] > 0.052) warning('N2 composition is greater than gas viscosity correlation upper limit: 0.052 mol fraction')
         if (nhc_composition[2] > 0.017) warning('H2S composition is greater than gas viscosity correlation upper limit: 0.017 mol fraction')
         if (nhc_composition[3] > 0.089) warning('CO2 composition is greater than gas viscosity correlation upper limit: 0.089 mol fraction')
         if ((gas_spgr < 0.554) | (gas_spgr > 1.861)) warning('gas specific gravity is outside of the gas viscosity correlation range: [0.554, 1.861]')
         if (p_ > 20000) warning('pressure is greater than gas viscosity correlation upper limit: 20000 psig')
         if (t_ > 652) warning('temperature is greater than gas viscosity correlation upper limit: 652 F')
      }
      oil_lst <- list(input_unit = input_unit, output_unit = output_unit, fluid = fluid, pvt_model = pvt_model, visc_model = visc_model, t = t, p = p, oil_api = oil_api, gas_spgr = gas_spgr, nhc_composition = nhc_composition, sat_cond = pb, flag_sat = 1, warning = warning)
      mat_results <- PVT_OIL_PROPERTIES(oil_lst)
   } else {
      if (rsi <= 0) stop("'rsi' must be greater than zero.")
      if (warning == "yes") {
         if (nhc_composition[1] > 0.217) warning('N2 composition is greater than dry gas PVT correlation upper limit: 0.217 mol fraction')
         if (nhc_composition[2] > 0.10) warning('H2S composition is greater than dry gas PVT correlation upper limit: 0.10 mol fraction')
         if (nhc_composition[3] > 0.558) warning('CO2 composition is greater than dry gas PVT correlation upper limit: 0.558 mol fraction')
         if ((oil_api < 16) | (oil_api > 58)) warning('oil API gravity is outside of the black_oil PVT correlation range: [16 API ,58 API]')
         if ((gas_spgr < 0.554) | (gas_spgr > 1.862)) warning('gas specific gravity is outside of the dry gas PVT correlation range: [0.554, 1.861]')
         if (p_ > 10000) warning('pressure is greater than dry gas PVT correlation upper limit: 10000 psig')
         if ((t_ < 32) | (t_ > 460)) warning('temperature is outside of the dry gas PVT correlation range: [32 F, 460 F]')
         if ((gas_spgr < 0.56) | (gas_spgr > 1.18)) warning('gas specific gravity is outside of the black_oil PVT correlation range: [0.56, 1.18]')
         if (p_ > 5250) warning('pressure is greater than black oil PVT correlation upper limit: 5250 psig')
         if ((t_ < 70) | (t_ > 295)) warning('temperature is outside of the black_oil PVT correlation range: [70 F, 295 F]')
         if (nhc_composition[1] > 0.052) warning('N2 composition is greater than gas viscosity correlation upper limit: 0.052 mol fraction')
         if (nhc_composition[2] > 0.017) warning('H2S composition is greater than gas viscosity correlation upper limit: 0.017 mol fraction')
         if (nhc_composition[3] > 0.089) warning('CO2 composition is greater than gas viscosity correlation upper limit: 0.089 mol fraction')
         if ((gas_spgr < 0.554) | (gas_spgr > 1.861)) warning('gas specific gravity is outside of the gas viscosity correlation range: [0.554, 1.861]')
         if (p_ > 20000) warning('pressure is greater than gas viscosity correlation upper limit: 20000 psig')
         if (t_ > 652) warning('temperature is greater than gas viscosity correlation upper limit: 652 F')
      }
      oil_lst <- list(input_unit = input_unit, output_unit = output_unit, fluid = fluid, pvt_model = pvt_model, visc_model = visc_model, t = t, p = p, oil_api = oil_api, gas_spgr = gas_spgr, nhc_composition = nhc_composition, sat_cond = rsi, flag_sat = 0, warning = warning)
      mat_results <- PVT_OIL_PROPERTIES(oil_lst)
   }
   return(mat_results)
}





#' Create a matrix of PVT properties for reservoir water samples
#'
#' The pvt_water() generates a table of water PVT properties at reservoir temperature and pressures from the atmospheric condition up to the initial reservoir pressure. The estimated water properties are solution gas-water ratio, formation volume factor, density, compressibility, and viscosity.
#'
#' @param input_unit input unit system for parameters, a character string either 'SI' or 'Field'.
#' @param output_unit output unit system for properties, a character string
#' @param fluid fluid type, the character string 'water'
#' @param pvt_model PVT model, a character string. 'Spivey', 'Meehan', and 'McCain' models are currently available
#' @param visc_model viscosity model, a character string. 'Spivey', 'Meehan', and 'McCain' models are currently available
#' @param t temperature, a numeric value either in 'C' or 'F' depending on the 'input_unit'
#' @param p pressure, a numeric value either in 'kPag' or 'Psig' depending on the 'input_unit'
#' @param gas_saturated a charater string either 'yes' or 'no'
#' @param salinity water salinity in weight percent TDS
#' @param warning a charater string either 'yes' or 'no'
#' @useDynLib Rpvt, .registration = TRUE
#'
#' @export
#'
#' @examples
#' pvt_water_results_1 <- pvt_water(input_unit = "Field", output_unit = "Field",
#' fluid = "water", pvt_model = "McCain", visc_model = "McCain",
#' t = 300, p = 5000, salinity = 10, gas_saturated = "yes", warning = "no")
#'
#' head(pvt_water_results_1)
#'
#' pvt_water_results_2 <- pvt_water(input_unit = "SI", output_unit = "SI",
#' fluid = "water", pvt_model = "Spivey", visc_model = "Spivey",
#' t = 100, p = 15000, salinity = 0.0, gas_saturated = "no", warning = "no")
#'
#' head(pvt_water_results_2)
#'
#' @references
#' \insertRef{Meehan1980a}{Rpvt}
#'
#' \insertRef{Meehan1980b}{Rpvt}
#'
#' \insertRef{McCainJr.1991}{Rpvt}
#'
#' \insertRef{Spivey2004}{Rpvt}
#'
#' \insertRef{McCainJr.2011}{Rpvt}

pvt_water <- function(input_unit = "Field", output_unit = "Filed", fluid = "water", pvt_model = "Spivey", visc_model = "Spivey", t = 220, p = 6000, salinity = 10.0, gas_saturated = "yes", warning = "yes") {

   kpa_to_psi <- 0.14503773800722
   if (!is.character(input_unit)) stop("'input_unit' must be a character string either 'SI' or 'Field'.")
   if (!is.character(output_unit)) stop("'output_unit' must be a character string either 'SI' or 'Field'.")
   if (!is.character(fluid)) stop("'fluid' must be the character string 'water'.")
   if (!is.character(pvt_model)) stop("'pvt_model' must be one of 'Spivey', 'Meehan' or 'McCain' correlations.")
   if (!(pvt_model %in% c("Meehan","McCain","Spivey"))) stop("'pvt_model' must be one of 'Spivey', 'Meehan' or 'McCain' correlations.")
   if (!is.character(visc_model)) stop("'visc_model' must be one of 'Spivey', 'Meehan' or 'McCain' correlations.")
   if (!(visc_model %in% c("Meehan","McCain","Spivey"))) stop("'visc_model' must be one of 'Spivey', 'Meehan' or 'McCain' correlations.")
   if (!(warning %in% c("no","yes"))) stop("'warning' must be a character string either 'yes' or 'no'.")
   if (!is.numeric(t)) stop("'t' must be a numeric.")
   if (t < 0) stop("'tempearture' must be a positive value.")
   if (!is.numeric(p)) stop("'p' must be a numeric.")
   if (p < 0) stop("'pressure' must be a positive value.")
   if (!is.character(gas_saturated)) stop("'gas_saturated' must be a character string either 'yes' or 'no'.")
   if (!(gas_saturated %in% c("no","yes"))) stop("'gas_saturated' must be assigned 'yes' or 'no'.")
   if (!is.numeric(salinity)) stop("'salinity' must be a numeric.")
   if (salinity < 0) stop("'salinity' must be equal or greater than zero for water samples.")
   if (input_unit == "SI") {
      p_ <- p * kpa_to_psi
      t_ <- t * 1.8 + 32
   } else {
      p_ <- p
      t_ <- t
   }
   if (warning == "yes") {
      if (pvt_model == "Spivey") {
         warning("The pressure lower limit in 'Spivey' correlation is 115 psig (790 kPag)")
         if (p_ > 29000) warning('pressure is greater than PVT correlation upper limit: 29000 psig')
         if (p_ < 130) warning('pressure is less than PVT correlation lower limit: 115 psig')
         if (t_ > 527) warning('temperature is greater than PVT correlation upper limit: 527 F')
         if (t_ < 68) warning('temperature is less than PVT correlation lower limit: 60 F')
         if (salinity > 25) warning('salinity is greater than PVT correlation upper limit: 25 wt%')
      }
      if (pvt_model == "McCain") {
         if (p_ > 5000) warning('pressure is greater than PVT correlation upper limit: 5000 psig')
         if (p_ < 1000) warning('pressure is less than PVT correlation lower limit for gas solubility estimation: 1000 psig')
         if (t_ > 260) warning('temperature is greater than PVT correlation upper limit: 260 F')
         if (salinity > 25) warning('salinity is greater than PVT correlation upper limit: 25 wt%')
      }
      if (pvt_model == "Meehan") {
         if (p_ > 6000) warning('pressure is greater than PVT correlation upper limit: 6000 psig')
         if (p_ < 1000) warning('pressure is less than PVT correlation lower limit: 1000 psig')
         if (t_ > 250) warning('temperature is greater than PVT correlation upper limit: 250 F')
         if (t_ < 80) warning('temperature is less than PVT correlation lower limit: 80 F')
         if (salinity > 25) warning('salinity is greater than PVT correlation upper limit: 25 wt%')
      }
   }
   water_lst <- list(input_unit = input_unit, output_unit = output_unit, fluid = fluid, pvt_model = pvt_model, visc_model = visc_model, t = t, p = p, salinity = salinity, gas_saturated = gas_saturated, warning = warning)
   mat_results <- PVT_WATER_PROPERTIES(water_lst)
   return(mat_results)
}

