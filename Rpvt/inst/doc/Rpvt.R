## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ---- eval = FALSE, warning = FALSE-------------------------------------------
#  install.packages("Rpvt")

## ---- fig.width = 6, fig.height= 8, fig.align = "center", warning = TRUE------
library(Rpvt)
library(ggplot2)
library(magrittr)
library(ggpubr)

pvt_gas_1 <- pvt_gas(input_unit = "Field", output_unit = "Field", 
fluid = "dry_gas", pvt_model = "DAK", visc_model = "Sutton", 
t = 400, p = 30000, gas_spgr = 0.65, nhc_composition = c(0.03,0.012,0.018),
cgr = 0, cond_api = NULL, warning = "yes")

attributes(pvt_gas_1)

pvt_gas_1 <- as.data.frame(pvt_gas_1)

head(pvt_gas_1,10)

colnames(pvt_gas_1) <- c("T(F)", "P(Psig)", "Z-FACTOR", "Bg(rb/scf)", 
                         "Density(lb/ft3)", "Cg(1/Psi)", "MUg(cp)", 
                         "m(p)(Psia^2/cp)")
Z_plot <- pvt_gas_1 %>% ggplot(aes(x = `P(Psig)`, y = `Z-FACTOR`)) +
  geom_point(color = "blue") +
  xlab(label = "P (Psig)") +
  ylab(label = "Z-Factor") +
  theme_bw()
Density_plot <- pvt_gas_1 %>% ggplot(aes(x = `P(Psig)`, y = `Density(lb/ft3)`)) +
  geom_point(color = "blue") +
  xlab(label = "P (Psig)") +
  ylab(label = "Density (lb/ft3)") +
  theme_bw()
Bg_plot <- pvt_gas_1 %>% ggplot(aes(x = `P(Psig)`, y = `Bg(rb/scf)`)) +
  geom_point(color = "blue") +
  scale_y_log10() +
  xlab(label = "P (Psig)") +
  ylab(label = "Gas FVF (rb/scf)") +
  theme_bw()
Cg_plot <- pvt_gas_1 %>% ggplot(aes(x = `P(Psig)`, y = `Cg(1/Psi)`)) +
  geom_point(color = "blue") +
  scale_y_log10() +
  xlab(label = "P (Psig)") +
  ylab(label = "Compressibility (1/Psi)") +
  theme_bw()
MUg_plot <- pvt_gas_1 %>% ggplot(aes(x = `P(Psig)`, y = `MUg(cp)`)) +
  geom_point(color = "blue") +
  xlab(label = "P (Psig)") +
  ylab(label = "Viscosity (cp)") +
  theme_bw()
mp_plot <- pvt_gas_1 %>% ggplot(aes(x = `P(Psig)`, y = `m(p)(Psia^2/cp)`)) +
  geom_point(color = "blue") +
  xlab(label = "P (Psig)") +
  ylab(label = "Pseudopressure (Psia^2/cp)") +
  theme_bw()

gas_pvt_plots <- ggarrange(Z_plot, Density_plot, Bg_plot, Cg_plot, MUg_plot, mp_plot,
                      ncol = 2, nrow = 3, align = "v")


gas_pvt_plots 

## -----------------------------------------------------------------------------

pvt_gas_2 <- pvt_gas(input_unit = "Field", output_unit = "SI", 
fluid = "wet_gas", pvt_model = "DAK", visc_model = "Sutton", 
t = 300, p = 5000, gas_spgr = 0.75, nhc_composition = c(0.05,0.01,0.04),
cgr = 5, cond_api = 42.3, warning = "no")

tail(as.data.frame(pvt_gas_2), 10)


## ---- fig.width = 6, fig.height= 12, fig.align = "center", warning = TRUE-----
library(Rpvt)
library(ggplot2)
library(magrittr)
library(ggpubr)

pvt_oil_1 <- pvt_oil(input_unit = "Field", output_unit = "Field", 
                                   fluid = "black_oil", pvt_model = "Standing",
                                   visc_model = "Beggs_Robinson", t = 200, p = 3000, 
                                   oil_api = 35, gas_spgr = 0.8,
                                   nhc_composition = c(0.05,0.02,0.04), 
                                   rsi = 650, pb = NULL, warning = "yes")

attributes(pvt_oil_1)

pvt_oil_1 <- as.data.frame(pvt_oil_1)

tail(pvt_oil_1, 10) 

Rs_plot <- pvt_oil_1 %>% ggplot(aes(x = `P_(Psig)`, y = `Rso_(scf/stb)`)) +
  geom_point(color = "blue", size = 2) +
  xlab(label = "P (Psig)") +
  ylab(label = "Rs (SCF/STB)") +
  theme_bw()
Bo_plot <- pvt_oil_1 %>% ggplot(aes(x = `P_(Psig)`, y = `Bo_(rb/stb)`)) +
  geom_point(color = "blue", size = 2) +
  xlab(label = "P (Psig)") +
  ylab(label = "Oil FVF (rb/STB)") +
  theme_bw()
RHOo_plot <- pvt_oil_1 %>% ggplot(aes(x = `P_(Psig)`, y = `Oil_Density_(lb/ft3)`)) +
  geom_point(color = "blue", size = 2) +
  xlab(label = "P (Psig)") +
  ylab(label = "Oil Density (lb/ft3)") +
  theme_bw()
Co_plot <- pvt_oil_1 %>% ggplot(aes(x = `P_(Psig)`, y = `Co_(1/Psia)`)) +
  geom_point(color = "blue", size = 2) +
  scale_y_log10() +
  xlab(label = "P (Psig)") +
  ylab(label = "Oil Compressibility (1/Psi)") +
  theme_bw()
MUo_plot <- pvt_oil_1 %>% ggplot(aes(x = `P_(Psig)`, y = `Oil_Viscosity_(cp)`)) +
  geom_point(color = "blue", size = 2) +
  scale_y_log10() +
  xlab(label = "P (Psig)") +
  ylab(label = "Oil Viscosity (cp)") +
  theme_bw()
Bg_plot <- pvt_oil_1 %>% ggplot(aes(x = `P_(Psig)`, y = `Bg_(rb/scf)`)) +
  geom_point(color = "blue", size = 2) +
  scale_y_log10() +
  xlab(label = "P (Psig)") +
  ylab(label = "Gas FVF (rb/SCF)") +
  theme_bw()
Cg_plot <- pvt_oil_1 %>% ggplot(aes(x = `P_(Psig)`, y = `Cg_(1/Psia)`)) +
  geom_point(color = "blue", size = 2) +
  scale_y_log10() +
  xlab(label = "P (Psig)") +
  ylab(label = "Gas Compressibility (1/Psi)") +
  theme_bw()
MUg_plot <- pvt_oil_1 %>% ggplot(aes(x = `P_(Psig)`, y = `Gas_Viscosity_(cp)`)) +
  geom_point(color = "blue", size = 2) +
  xlab(label = "P (Psig)") +
  ylab(label = "Gas Viscosity (cp)") +
  theme_bw()
oil_pvt_plots <- ggarrange(Rs_plot, RHOo_plot, Bo_plot, Bg_plot, Co_plot, Cg_plot, MUo_plot, MUg_plot,
                      ncol = 2, nrow = 4, align = "v")

oil_pvt_plots

## -----------------------------------------------------------------------------

pvt_oil_2 <- pvt_oil(input_unit = "SI", output_unit = "SI", fluid = "black_oil",
                     pvt_model = "Vasquez_Beggs", visc_model = "Al_Marhoun", 
                     t = 100, p = 20000, oil_api = 40, gas_spgr = 0.75, 
                     nhc_composition = c(0.05,0.02,0.04), 
                     rsi = NULL, pb = 13000, warning = "yes")

head(as.data.frame(pvt_oil_2), 10)


## -----------------------------------------------------------------------------

pvt_oil_3 <- pvt_oil(input_unit = "Field", output_unit = "SI", fluid = "black_oil",
                     pvt_model = "Farshad_Petrosky", visc_model = "Al_Marhoun", 
                     t = 260, p = 4000, oil_api = 38, gas_spgr = 0.68, 
                     nhc_composition = c(0.03,0.07,0.08), 
                     rsi = NULL, pb = 2500, warning = "yes")

head(as.data.frame(pvt_oil_3), 10)


## ---- fig.width = 6, fig.height= 8, fig.align = "center", warning = TRUE------
library(Rpvt)
library(ggplot2)
library(magrittr)
library(ggpubr)

pvt_water_1 <- as.data.frame(pvt_water(input_unit = "SI", output_unit = "SI", 
fluid = "water", pvt_model = "Spivey", visc_model = "Spivey", t = 150, p = 200000, 
salinity = 5, gas_saturated = "yes", warning = "yes"))

head(pvt_water_1,10)

colnames(pvt_water_1) <- c("T(C)", "P(kPag)", "Rsw(rm3/sm3)", "Bw(rm3/sm3)", 
                         "Density(kg/m3)", "Cw(1/kPa)", "MUw(mPa.s)")
Solubility_plot <- pvt_water_1 %>% ggplot(aes(x = `P(kPag)`, y = `Rsw(rm3/sm3)`)) +
  geom_point(color = "blue") +
  xlab(label = "P (kPag)") +
  ylab(label = "Gas solubility in water (rm3/sm3)") +
  theme_bw()
Density_plot <- pvt_water_1 %>% ggplot(aes(x = `P(kPag)`, y = `Density(kg/m3)`)) +
  geom_point(color = "blue") +
  xlab(label = "P (kPag)") +
  ylab(label = "Density (kg/m3)") +
  theme_bw()
Bw_plot <- pvt_water_1 %>% ggplot(aes(x = `P(kPag)`, y = `Bw(rm3/sm3)`)) +
  geom_point(color = "blue") +
  scale_y_log10() +
  xlab(label = "P (kPag)") +
  ylab(label = "Water FVF (rm3/sm3)") +
  theme_bw()
Cw_plot <- pvt_water_1 %>% ggplot(aes(x = `P(kPag)`, y = `Cw(1/kPa)`)) +
  geom_point(color = "blue") +
  scale_y_log10() +
  xlab(label = "P (kPag)") +
  ylab(label = "Compressibility (1/kPa)") +
  theme_bw()
MUw_plot <- pvt_water_1 %>% ggplot(aes(x = `P(kPag)`, y = `MUw(mPa.s)`)) +
  geom_point(color = "blue") +
  xlab(label = "P (kPag)") +
  ylab(label = "Viscosity (mPa.s)") +
  theme_bw()
water_pvt_plots <- ggarrange(Solubility_plot, Density_plot, Bw_plot, Cw_plot, MUw_plot,
                      ncol = 2, nrow = 3, align = "v")

water_pvt_plots 

## -----------------------------------------------------------------------------

pvt_water_2 <- as.data.frame(pvt_water(input_unit = "Field", output_unit = "SI", 
fluid = "water", pvt_model = "McCain", visc_model = "McCain", t = 200, p = 4000, 
salinity = 15, gas_saturated = "yes", warning = "yes"))

head(pvt_water_2, 10)


## -----------------------------------------------------------------------------

pvt_water_3 <- pvt_water(input_unit = "Field", output_unit = "Field",
fluid = "water", pvt_model = "Meehan", visc_model = "Meehan",
t = 300, p = 3000,  salinity = 10, gas_saturated = "no", warning = "yes")

tail(as.data.frame(pvt_water_3), 10)


