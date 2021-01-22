context("test-check")

test_that("error if gas pvt_model is not 'DAK'", {

   expect_error(pvt_gas(input_unit = "Field", output_unit = "Field",
                        fluid = "dry_gas", pvt_model = "Sutton", visc_model = "Sutton",
                        t = 400, p = 20000, gas_spgr = 0.65, nhc_composition = c(0.05,0.02,0.04),
                        cgr = 0, cond_api = NULL, warning = "yes"))

})


test_that("error if gas visc_model is not 'Sutton'", {

   expect_error(pvt_gas(input_unit = "Field", output_unit = "Field",
                        fluid = "dry_gas", pvt_model = "DAK", visc_model = "DAK",
                        t = 400, p = 20000, gas_spgr = 0.65, nhc_composition = c(0.05,0.02,0.04),
                        cgr = 0, cond_api = NULL, warning = "yes"))

})


test_that("error if a different function is used for gas PVT results", {

   expect_error(pvt_gass(input_unit = "Field", output_unit = "Field",
                        fluid = "dry_gas", pvt_model = "DAK", visc_model = "Sutton",
                        t = 400, p = 20000, gas_spgr = 0.65, nhc_composition = c(0.05,0.02,0.04),
                        cgr = 0, cond_api = NULL, warning = "yes"))

})


test_that("warning if pressure is out of the range", {

   expect_warning(pvt_gas(input_unit = "Field", output_unit = "Field",
                        fluid = "dry_gas", pvt_model = "DAK", visc_model = "Sutton",
                        t = 400, p = 30000, gas_spgr = 0.65, nhc_composition = c(0.05,0.02,0.04),
                        cgr = 0, cond_api = NULL, warning = "yes"))

})
