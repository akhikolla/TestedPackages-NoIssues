
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Rpvt

`Rpvt` is a correlation-based PVT (Pressure-Volume-Temperature) package
for dry gas, wet gas, black oil, and water samples. It generates PVT
properties of hydrocarbons and water samples in a tabular format at a
constant temperature from atmospheric pressure up to the pressure of
interest.

## Installation

You can install the released version of Rpvt from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("Rpvt")
```

## Example

### Dry Gas Example

``` r
library(Rpvt)

pvt_gas_results <- pvt_gas(input_unit = "Field", output_unit = "Field", 
                           fluid = "dry_gas", pvt_model = "DAK", visc_model = "Sutton", 
                           t = 400, p = 30000, gas_spgr = 0.65, 
                           nhc_composition = c(0.03,0.012,0.018), 
                           cgr = 0, cond_api = NULL, warning = "yes")
#> Warning in pvt_gas(input_unit = "Field", output_unit = "Field", fluid =
#> "dry_gas", : pressure is greater than dry gas PVT correlation upper limit: 10000
#> psig
#> Warning in pvt_gas(input_unit = "Field", output_unit = "Field", fluid =
#> "dry_gas", : pressure is greater than gas viscosity correlation upper limit:
#> 20000 psig

attributes(pvt_gas_results)
#> $dim
#> [1] 3001    8
#> 
#> $dimnames
#> $dimnames[[1]]
#> NULL
#> 
#> $dimnames[[2]]
#> [1] "T_(F)"             "P_(Psig)"          "Z-Factor"         
#> [4] "Bg_(rb/scf)"       "Density_(lb/cuft)" "Cg_(1/Psia)"      
#> [7] "Viscosity_(cp)"    "m(p)_(Psia^2/cp)" 
#> 
#> 
#> $`gas pseudocritical temperature (F)`
#> [1] -102.2183
#> 
#> $`gas pseudocritical pressure (Psia)`
#> [1] 648.5108

pvt_gas_results <- as.data.frame(pvt_gas_results)

head(pvt_gas_results,10)
#>    T_(F) P_(Psig)  Z-Factor Bg_(rb/scf) Density_(lb/cuft) Cg_(1/Psia)
#> 1    400        0 0.9996026  0.29449758        0.03000184 0.068072651
#> 2    400       10 0.9993343  0.17520145        0.05043034 0.040519143
#> 3    400       20 0.9990678  0.12467204        0.07086968 0.028848352
#> 4    400       30 0.9988031  0.09675307        0.09131977 0.022399782
#> 5    400       40 0.9985402  0.07904301        0.11178051 0.018309117
#> 6    400       50 0.9982790  0.06680793        0.13225179 0.015482979
#> 7    400       60 0.9980196  0.05784892        0.15273351 0.013413498
#> 8    400       70 0.9977620  0.05100557        0.17322557 0.011832661
#> 9    400       80 0.9975062  0.04560763        0.19372786 0.010585663
#> 10   400       90 0.9972522  0.04124093        0.21424029 0.009576844
#>    Viscosity_(cp) m(p)_(Psia^2/cp)
#> 1      0.01652119             0.00
#> 2      0.01652162         23856.42
#> 3      0.01652218         59833.30
#> 4      0.01652286        107935.97
#> 5      0.01652364        168169.44
#> 6      0.01652452        240538.49
#> 7      0.01652550        325047.64
#> 8      0.01652656        421701.16
#> 9      0.01652770        530503.09
#> 10     0.01652893        651457.25
```

### Black Oil Example

``` r
library(Rpvt)

pvt_oil_results <- pvt_oil(input_unit = "Field", output_unit = "Field", 
                                   fluid = "black_oil", pvt_model = "Standing",
                                   visc_model = "Beggs_Robinson", t = 200, p = 3000, 
                                   oil_api = 35, gas_spgr = 0.8,
                                   nhc_composition = c(0.05,0.02,0.04), 
                                   rsi = 650, pb = NULL, warning = "yes")
#> Warning in pvt_oil(input_unit = "Field", output_unit = "Field", fluid =
#> "black_oil", : H2S composition is greater than gas viscosity correlation upper
#> limit: 0.017 mol fraction

attributes(pvt_oil_results)
#> $dim
#> [1] 301  13
#> 
#> $dimnames
#> $dimnames[[1]]
#> NULL
#> 
#> $dimnames[[2]]
#>  [1] "T_(F)"                "P_(Psig)"             "Rso_(scf/stb)"       
#>  [4] "Bo_(rb/stb)"          "Oil_Density_(lb/ft3)" "Co_(1/Psia)"         
#>  [7] "Oil_Viscosity_(cp)"   "Z-Factor"             "Bg_(rb/scf)"         
#> [10] "Gas_Density_(lb/ft3)" "Cg_(1/Psia)"          "Gas_Viscosity_(cp)"  
#> [13] "m(p)_(Psia^2/cp)"    
#> 
#> 
#> $`gas pseudocritical temperature (F)`
#> [1] -61.59983
#> 
#> $`gas pseudocritical pressure (Psia)`
#> [1] 647.0921

pvt_oil_results <- as.data.frame(pvt_oil_results)

tail(pvt_oil_results, 10) 
#>     T_(F) P_(Psig) Rso_(scf/stb) Bo_(rb/stb) Oil_Density_(lb/ft3)  Co_(1/Psia)
#> 292   200     2910           650    1.380185             43.52833 1.319029e-05
#> 293   200     2920           650    1.380014             43.53371 1.316513e-05
#> 294   200     2930           650    1.379844             43.53907 1.314019e-05
#> 295   200     2940           650    1.379675             43.54441 1.311547e-05
#> 296   200     2950           650    1.379506             43.54974 1.309097e-05
#> 297   200     2960           650    1.379338             43.55504 1.306667e-05
#> 298   200     2970           650    1.379171             43.56034 1.304258e-05
#> 299   200     2980           650    1.379004             43.56561 1.301870e-05
#> 300   200     2990           650    1.378837             43.57087 1.299502e-05
#> 301   200     3000           650    1.378671             43.57611 1.297154e-05
#>     Oil_Viscosity_(cp)  Z-Factor  Bg_(rb/scf) Gas_Density_(lb/ft3)  Cg_(1/Psia)
#> 292          0.5075536 0.8571989 0.0009737542             11.16752 0.0003017984
#> 293          0.5080544 0.8575451 0.0009708281             11.20118 0.0003001107
#> 294          0.5085569 0.8578959 0.0009679270             11.23475 0.0002984348
#> 295          0.5090612 0.8582513 0.0009650507             11.26824 0.0002967708
#> 296          0.5095672 0.8586112 0.0009621989             11.30164 0.0002951185
#> 297          0.5100750 0.8589757 0.0009593714             11.33495 0.0002934779
#> 298          0.5105844 0.8593446 0.0009565678             11.36817 0.0002918488
#> 299          0.5110956 0.8597181 0.0009537878             11.40130 0.0002902311
#> 300          0.5116085 0.8600959 0.0009510313             11.43435 0.0002886249
#> 301          0.5121231 0.8604782 0.0009482979             11.46731 0.0002870300
#>     Gas_Viscosity_(cp) m(p)_(Psia^2/cp)
#> 292         0.01894741        613123202
#> 293         0.01898046        616726950
#> 294         0.01901354        620335238
#> 295         0.01904666        623948017
#> 296         0.01907980        627565235
#> 297         0.01911298        631186844
#> 298         0.01914620        634812793
#> 299         0.01917944        638443035
#> 300         0.01921272        642077519
#> 301         0.01924602        645716199
```

### Water Example

``` r
library(Rpvt)

pvt_water_results <- pvt_water(input_unit = "SI", output_unit = "SI", 
                               fluid = "water", pvt_model = "Spivey", 
                               visc_model = "Spivey", t = 150, p = 200000, 
                               salinity = 5, gas_saturated = "yes", warning = "yes")
#> Warning in pvt_water(input_unit = "SI", output_unit = "SI", fluid = "water", :
#> The pressure lower limit in 'Spivey' correlation is 115 psig (790 kPag)
#> Warning in pvt_water(input_unit = "SI", output_unit = "SI", fluid = "water", :
#> pressure is greater than PVT correlation upper limit: 29000 psig

head(pvt_water_results,10)
#>       T_(C)  P_(kPag) Rsw_(rm3/sm3) Bw_(rm3/sm3) Density_(kg/m3)  Cw_(1/kPaa)
#>  [1,]   150   0.00000           NaN          NaN             NaN          NaN
#>  [2,]   150  68.94757           NaN          NaN             NaN          NaN
#>  [3,]   150 137.89515           NaN          NaN             NaN          NaN
#>  [4,]   150 206.84272           NaN          NaN             NaN          NaN
#>  [5,]   150 275.79029           NaN          NaN             NaN          NaN
#>  [6,]   150 344.73786           NaN          NaN             NaN          NaN
#>  [7,]   150 413.68544   0.001877106     1.086347        952.9482 2.160339e-05
#>  [8,]   150 482.63301   0.009433052     1.086319        952.9771 3.105915e-05
#>  [9,]   150 551.58058   0.019506617     1.086296        953.0035 3.324147e-05
#> [10,]   150 620.52816   0.031017232     1.086276        953.0285 3.317271e-05
#>       Viscosity_(mPa.s)
#>  [1,]         0.2066466
#>  [2,]         0.2066671
#>  [3,]         0.2066877
#>  [4,]         0.2067082
#>  [5,]         0.2067288
#>  [6,]         0.2067493
#>  [7,]         0.2067698
#>  [8,]         0.2067904
#>  [9,]         0.2068109
#> [10,]         0.2068314
```
