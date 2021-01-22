#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Main Formula Creation for \code{modelsearch()}
//'
//' Function \code{stovokor()} creates formulae to be used as input in the global
//' model calls used in function \code{\link{modelsearch}()}.
//'
//' @param surv A vector of strings indicating the names of the variables coding 
//' survival.
//' @param obs A vector of strings indicating the names of the variables coding 
//' observation status.
//' @param size A vector of strings indicating the names of the variables coding 
//' size.
//' @param repst A vector of strings indicating the names of the variables coding 
//' reproductive status.
//' @param fec A vector of strings indicating the names of the variables coding 
//' fecundity.
//' @param vitalrates A vector of strings indicating which vital rates will be
//' estimated.
//' @param historical A logical value indicating whether to create global models
//' with historical effects.
//' @param suite A string indicating the scope of independent factors included in
//' the global models. Options include \code{"full"}, \code{"main"}, \code{"size"}, \code{"rep"},
//' and \code{"const"}.
//' @param approach A string indicating whether to use mixed model encoding 
//' (\code{"mixed"}) or GLM encoding (\code{"glm"}).
//' @param sizedist A string variable indicating the distribution to use to
//' model size.
//' @param fecdist A string variable indicating the distribution to use to
//' model fecundity.
//' @param nojuvs A logical value indicating that juvenile rates should be
//' estimated (\code{FALSE}) or not (\code{TRUE}).
//' @param age A string indicating the name of the variable coding age.
//' @param indcova A vector of strings indicating the names in times \emph{t}+1, \emph{t},
//' and \emph{t}-1 of a specific individual covariate used in the dataset.
//' @param indcovb A vector of strings indicating the names in times \emph{t}+1, \emph{t},
//' and \emph{t}-1 of a specific individual covariate used in the dataset.
//' @param indcovc A vector of strings indicating the names in times \emph{t}+1, \emph{t},
//' and \emph{t}-1 of a specific individual covariate used in the dataset.
//' @param indiv A string indicating the name of the variable coding individual
//' identity.
//' @param patch A string indicating the name of the variable coding patch identity.
//' @param year A string indicating the name of the variable coding time \emph{t}.
//' @param pasrand A logical value indicating whether to treat patch as a random
//' variable within mixed models.
//' @param yasrand A logical value indicating whether to treat year as a random
//' variable within mixed models.
//' @param fectime An integer indicating whether to use reproductive output in time 
//' \emph{t} (2) or time \emph{t}+1 (3) as the response for fecundity.
//' @param juvsize A logical value indicating whether to include size terms in
//' juvenile models.
//' @param size0 A logical value indicating whether size distribution should
//' be zero-inflated. Only applies to Poisson and negative binomial distributions.
//' @param fec0 A logical value indicating whether size distribution should
//' be zero-inflated. Only applies to Poisson and negative binomial distributions.
//' 
//' @return Vector of 9 strings, each a formula to be used as input in function.
//' \code{modelsearch()}.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List stovokor(StringVector surv, StringVector obs, StringVector size, StringVector repst, 
              StringVector fec, StringVector vitalrates, bool historical, String suite, 
              String approach, String sizedist, String fecdist, bool nojuvs, String age, 
              StringVector indcova, StringVector indcovb, StringVector indcovc,
              String indiv, String patch, String year, bool pasrand, bool yasrand, 
              int fectime, bool juvsize, bool size0, bool fec0) {
  
  if (nojuvs) juvsize = FALSE;
  
  int nvitalrates = vitalrates.length();
  bool survcheck = 0;
  bool obscheck = 0;
  bool sizecheck = 0;
  bool repstcheck = 0;
  bool feccheck = 0;
  
  int sizel = size.length();
  int repstl = repst.length();
  
  String fullsurvmodel;
  String fullobsmodel;
  String fullsizemodel;
  String fullrepstmodel;
  String fullfecmodel;
  String juvsurvmodel;
  String juvobsmodel;
  String juvsizemodel;
  String juvrepstmodel;
  String randomtackonp;
  String randomtackony;
  String randomtackoni;
  String randomtackon;
  String fixedtackonp;
  String fixedtackony;
  String fixedtackon;
  String sizesuffix;
  String jsizesuffix;
  String fecsuffix;
  
  int covcount {0};
  if (indcova(1) != "none") covcount += 1;
  if (indcovb(1) != "none") covcount += 1;
  if (indcovc(1) != "none") covcount += 1;
  
  // This section determines which vital rates need global model formulae
  for (int i = 0; i < nvitalrates; i++) {
    if (vitalrates(i) == "surv") survcheck = 1;
    if (vitalrates(i) == "obs") obscheck = 1;
    if (vitalrates(i) == "size") sizecheck = 1;
    if (vitalrates(i) == "repst") repstcheck = 1;
    if (vitalrates(i) == "fec") feccheck = 1;
  }
  
  // This section tests to see if the inputs are appropriate for the suite
  if (suite == "full" || suite == "main") {
    if (historical) {
      if (sizel != 3 && repstl != 3) {
        if (sizel == 2 && repstl == 2) {
          historical = FALSE;
        } else if (repstl == 3) {
          suite = "rep";
        } else if (sizel == 3) {
          suite = "size";
        }
      }
    } else {
      if (sizel != 2 && repstl != 2) {
        if (sizel != 3 && repstl != 3) {
          suite = "const";
        } else if (repstl > 1 && sizel < 2) {
          suite = "rep";
        } else if (sizel > 1 && repstl < 2) {
          suite = "size";
        }
      }
    }
  } else if (suite == "rep") {
    if (historical) {
      if (repstl != 3) {
        if (repstl == 2) {
          historical = FALSE;
        } else {
          suite = "const";
        }
      }
    } else {
      if (repstl != 2) {
        if (repstl != 3) {
          suite = "const";
        }
      }
    }
  } else if (suite == "size") {
    if (historical) {
      if (sizel != 3) {
        if (sizel == 2) {
          historical = FALSE;
        } else {
          suite = "const";
        }
      }
    } else {
      if (sizel != 2) {
        if (sizel != 3) {
          suite = "const";
        }
      }
    }
  }
  
  // Here we determine the nature of the potentially random variables
  randomtackonp = "";
  randomtackony = "";
  randomtackoni = "";
  fixedtackonp = "";
  fixedtackony = "";
  
  if (approach != "mixed") {
    yasrand = FALSE;
    pasrand = FALSE;
  }
  
  if (year!= "none") {
    if (yasrand) {
      randomtackony += " + ";
      randomtackony += "(1 | ";
      randomtackony += year;
      randomtackony += ")";
    } else {
      fixedtackony += " + as.factor(";
      fixedtackony += year;
      fixedtackony += ")";
    }
  }
  
  if (patch!= "none") {
    if (pasrand) {
      randomtackonp += " + ";
      randomtackonp += "(1 | ";
      randomtackonp += patch;
      randomtackonp += ")";
    } else {
      fixedtackonp += " + as.factor(";
      fixedtackonp += patch;
      fixedtackonp += ")";
    }
  }
  
  if (indiv != "none" && approach == "mixed") {
    randomtackoni += " + (1 | ";
    randomtackoni += indiv;
    randomtackoni += ")";
  }
  
  randomtackon = randomtackony;
  randomtackon += randomtackonp;
  randomtackon += randomtackoni;
  fixedtackon = fixedtackony;
  fixedtackon += fixedtackonp;
  
  // Now we will build the global models
  if (survcheck) {
    fullsurvmodel = surv(0);
    fullsurvmodel += " ~ ";
    
    if (!nojuvs) {
      juvsurvmodel = surv(0);
      juvsurvmodel += " ~ ";
      
      if (suite == "full" || suite == "main" || suite == "size") {
        if (juvsize) {
          juvsurvmodel += size(1);
        } else juvsurvmodel += "1";
      } else juvsurvmodel += "1";
      
      juvsurvmodel += fixedtackon;
      juvsurvmodel += randomtackon;
    }
    
    if (age != "none") {
      fullsurvmodel += age;
    }
    if (indcova(1) != "none") {
      if (age != "none") fullsurvmodel += " + ";
      fullsurvmodel += indcova(1);
      
      if (historical) {
        fullsurvmodel += " + ";
        fullsurvmodel += indcova(2);
      }
    }
    if (indcovb(1) != "none") {
      if (age != "none" || indcova(1) != "none") fullsurvmodel += " + ";
      fullsurvmodel += indcovb(1);
      
      if (historical) {
        fullsurvmodel += " + ";
        fullsurvmodel += indcovb(2);
      }
    }
    if (indcovc(1) != "none") {
      if (age != "none" || covcount > 1) fullsurvmodel += " + ";
      fullsurvmodel += indcovc(1);
      
      if (historical) {
        fullsurvmodel += " + ";
        fullsurvmodel += indcovc(2);
      }
    }
    
    if (suite == "full") {
      if (age != "none" || covcount > 0) fullsurvmodel += " + ";
      fullsurvmodel += size(1);
      fullsurvmodel += " + ";
      fullsurvmodel += repst(1);
      
      if (historical) {
        fullsurvmodel += " + ";
        fullsurvmodel += size(2);
        fullsurvmodel += " + ";
        fullsurvmodel += repst(2);
      }
      
      fullsurvmodel += " + ";
      fullsurvmodel += size(1);
      fullsurvmodel += ":";
      fullsurvmodel += repst(1);
      
      if (historical) {
        fullsurvmodel += " + ";
        fullsurvmodel += size(1);
        fullsurvmodel += ":";
        fullsurvmodel += size(2);
        
        fullsurvmodel += " + ";
        fullsurvmodel += repst(1);
        fullsurvmodel += ":";
        fullsurvmodel += repst(2);
        
        fullsurvmodel += " + ";
        fullsurvmodel += size(1);
        fullsurvmodel += ":";
        fullsurvmodel += repst(2);
        
        fullsurvmodel += " + ";
        fullsurvmodel += repst(1);
        fullsurvmodel += ":";
        fullsurvmodel += size(2);
      }
      
      if (age != "none") {
        fullsurvmodel += " + ";
        fullsurvmodel += age;
        fullsurvmodel += ":";
        fullsurvmodel += size(1);
        
        fullsurvmodel += " + ";
        fullsurvmodel += age;
        fullsurvmodel += ":";
        fullsurvmodel += repst(1);
        
        if (historical) {
          fullsurvmodel += " + ";
          fullsurvmodel += age;
          fullsurvmodel += ":";
          fullsurvmodel += size(2);
          
          fullsurvmodel += " + ";
          fullsurvmodel += age;
          fullsurvmodel += ":";
          fullsurvmodel += repst(2);
        }
      }
      
      if (indcova(1) != "none") {
        fullsurvmodel += " + ";
        fullsurvmodel += indcova(1);
        fullsurvmodel += ":";
        fullsurvmodel += size(1);
        
        fullsurvmodel += " + ";
        fullsurvmodel += indcova(1);
        fullsurvmodel += ":";
        fullsurvmodel += repst(1);
        
        if (historical && indcova(2) != "none") {
          fullsurvmodel += " + ";
          fullsurvmodel += indcova(2);
          fullsurvmodel += ":";
          fullsurvmodel += size(2);
          
          fullsurvmodel += " + ";
          fullsurvmodel += indcova(2);
          fullsurvmodel += ":";
          fullsurvmodel += repst(2);
        }
      }
      if (indcovb(1) != "none") {
        fullsurvmodel += " + ";
        fullsurvmodel += indcovb(1);
        fullsurvmodel += ":";
        fullsurvmodel += size(1);
        
        fullsurvmodel += " + ";
        fullsurvmodel += indcovb(1);
        fullsurvmodel += ":";
        fullsurvmodel += repst(1);
        
        if (historical && indcovb(2) != "none") {
          fullsurvmodel += " + ";
          fullsurvmodel += indcovb(2);
          fullsurvmodel += ":";
          fullsurvmodel += size(2);
          
          fullsurvmodel += " + ";
          fullsurvmodel += indcovb(2);
          fullsurvmodel += ":";
          fullsurvmodel += repst(2);
        }
        
        if (indcova(1) != "none") {
          fullsurvmodel += " + ";
          fullsurvmodel += indcova(1);
          fullsurvmodel += ":";
          fullsurvmodel += indcovb(1);
          
          if (historical) {
            if (indcova(2) != "none" && indcovb(2) != "none") {
              fullsurvmodel += " + ";
              fullsurvmodel += indcova(2);
              fullsurvmodel += ":";
              fullsurvmodel += indcovb(2);
              
              fullsurvmodel += " + ";
              fullsurvmodel += indcova(1);
              fullsurvmodel += ":";
              fullsurvmodel += indcovb(2);
              
              fullsurvmodel += " + ";
              fullsurvmodel += indcova(2);
              fullsurvmodel += ":";
              fullsurvmodel += indcovb(1);
            } else if (indcova(2) != "none") {
              fullsurvmodel += " + ";
              fullsurvmodel += indcova(2);
              fullsurvmodel += ":";
              fullsurvmodel += indcovb(1);
            } else if (indcovb(2) != "none") {
              fullsurvmodel += " + ";
              fullsurvmodel += indcova(1);
              fullsurvmodel += ":";
              fullsurvmodel += indcovb(2);
            }
          } 
        }
      }
      
      if (indcovc(1) != "none") {
        fullsurvmodel += " + ";
        fullsurvmodel += indcovc(1);
        fullsurvmodel += ":";
        fullsurvmodel += size(1);
        
        fullsurvmodel += " + ";
        fullsurvmodel += indcovc(1);
        fullsurvmodel += ":";
        fullsurvmodel += repst(1);
        
        if (historical && indcovc(2) != "none") {
          fullsurvmodel += " + ";
          fullsurvmodel += indcovc(2);
          fullsurvmodel += ":";
          fullsurvmodel += size(2);
          
          fullsurvmodel += " + ";
          fullsurvmodel += indcovc(2);
          fullsurvmodel += ":";
          fullsurvmodel += repst(2);
        }
        
        if (indcova(1) != "none") {
          fullsurvmodel += " + ";
          fullsurvmodel += indcova(1);
          fullsurvmodel += ":";
          fullsurvmodel += indcovc(1);
          
          if (historical) {
            if (indcova(2) != "none" && indcovc(2) != "none") {
              fullsurvmodel += " + ";
              fullsurvmodel += indcova(2);
              fullsurvmodel += ":";
              fullsurvmodel += indcovc(2);
              
              fullsurvmodel += " + ";
              fullsurvmodel += indcova(1);
              fullsurvmodel += ":";
              fullsurvmodel += indcovc(2);
              
              fullsurvmodel += " + ";
              fullsurvmodel += indcova(2);
              fullsurvmodel += ":";
              fullsurvmodel += indcovc(1);
            } else if (indcova(2) != "none") {
              fullsurvmodel += " + ";
              fullsurvmodel += indcova(2);
              fullsurvmodel += ":";
              fullsurvmodel += indcovc(1);
            } else if (indcovc(2) != "none") {
              fullsurvmodel += " + ";
              fullsurvmodel += indcova(1);
              fullsurvmodel += ":";
              fullsurvmodel += indcovc(2);
            }
          } 
        }
        if (indcovb(1) != "none") {
          fullsurvmodel += " + ";
          fullsurvmodel += indcovb(1);
          fullsurvmodel += ":";
          fullsurvmodel += indcovc(1);
          
          if (historical) {
            if (indcovb(2) != "none" && indcovc(2) != "none") {
              fullsurvmodel += " + ";
              fullsurvmodel += indcovb(2);
              fullsurvmodel += ":";
              fullsurvmodel += indcovc(2);
              
              fullsurvmodel += " + ";
              fullsurvmodel += indcovb(1);
              fullsurvmodel += ":";
              fullsurvmodel += indcovc(2);
              
              fullsurvmodel += " + ";
              fullsurvmodel += indcovb(2);
              fullsurvmodel += ":";
              fullsurvmodel += indcovc(1);
            } else if (indcovb(2) != "none") {
              fullsurvmodel += " + ";
              fullsurvmodel += indcovb(2);
              fullsurvmodel += ":";
              fullsurvmodel += indcovc(1);
            } else if (indcovc(2) != "none") {
              fullsurvmodel += " + ";
              fullsurvmodel += indcovb(1);
              fullsurvmodel += ":";
              fullsurvmodel += indcovc(2);
            }
          } 
        }
      }
    } else if (suite == "main") {
      if (age != "none" || covcount > 0) fullsurvmodel += " + ";
      fullsurvmodel += size(1);
      
      if (historical && sizecheck) {
        fullsurvmodel += " + ";
        fullsurvmodel += size(2);
      }
      
      fullsurvmodel += " + ";
      fullsurvmodel += repst(1);
      
      if (historical && repstcheck) {
        fullsurvmodel += " + ";
        fullsurvmodel += repst(2);
      }
      
    } else if (suite == "size") {
      if (age != "none" || covcount > 0) fullsurvmodel += " + ";
      fullsurvmodel += size(1);
      
      if (historical && sizecheck) {
        fullsurvmodel += " + ";
        fullsurvmodel += size(2);
        
        fullsurvmodel += " + ";
        fullsurvmodel += size(1);
        fullsurvmodel += ":";
        fullsurvmodel += size(2);
      }
      
    } else if (suite == "rep") {
      if (age != "none" || covcount > 0) fullsurvmodel += " + ";
      fullsurvmodel += repst(1);
      
      if (historical && repstcheck) {
        fullsurvmodel += " + ";
        fullsurvmodel += repst(2);
        
        fullsurvmodel += " + ";
        fullsurvmodel += repst(1);
        fullsurvmodel += ":";
        fullsurvmodel += repst(2);
      }
      
    } else {
      if (age == "none" && covcount == 0) fullsurvmodel += "1";
    }
    
    fullsurvmodel += fixedtackon;
    fullsurvmodel += randomtackon;
    
  } else {
    fullsurvmodel = "none";
  }
  
  if (obscheck) {
    fullobsmodel = obs(0);
    fullobsmodel += " ~ ";
    
    if (!nojuvs) {
      juvobsmodel = obs(0);
      juvobsmodel += " ~ ";
      
      if (suite == "full" || suite == "main" || suite == "size") {
        if (juvsize) {
          juvobsmodel += size(1);
        } else juvobsmodel += "1";
      } else juvobsmodel += "1";
      
      juvobsmodel += fixedtackon;
      juvobsmodel += randomtackon;
    }
    
    if (age != "none") {
      fullobsmodel += age;
    }
    if (indcova(1) != "none") {
      if (age != "none") fullobsmodel += " + ";
      fullobsmodel += indcova(1);
      
      if (historical) {
        fullobsmodel += " + ";
        fullobsmodel += indcova(2);
      }
    }
    if (indcovb(1) != "none") {
      if (age != "none" || indcova(1) != "none") fullobsmodel += " + ";
      fullobsmodel += indcovb(1);
      
      if (historical) {
        fullobsmodel += " + ";
        fullobsmodel += indcovb(2);
      }
    }
    if (indcovc(1) != "none") {
      if (age != "none" || covcount > 1) fullobsmodel += " + ";
      fullobsmodel += indcovc(1);
      
      if (historical) {
        fullobsmodel += " + ";
        fullobsmodel += indcovc(2);
      }
    }
    
    if (suite == "full") {
      if (age != "none" || covcount > 0) fullobsmodel += " + ";
      fullobsmodel += size(1);
      fullobsmodel += " + ";
      fullobsmodel += repst(1);
      
      if (historical) {
        fullobsmodel += " + ";
        fullobsmodel += size(2);
        fullobsmodel += " + ";
        fullobsmodel += repst(2);
      }
      
      fullobsmodel += " + ";
      fullobsmodel += size(1);
      fullobsmodel += ":";
      fullobsmodel += repst(1);
      
      if (historical) {
        fullobsmodel += " + ";
        fullobsmodel += size(1);
        fullobsmodel += ":";
        fullobsmodel += size(2);
        
        fullobsmodel += " + ";
        fullobsmodel += repst(1);
        fullobsmodel += ":";
        fullobsmodel += repst(2);
        
        fullobsmodel += " + ";
        fullobsmodel += size(1);
        fullobsmodel += ":";
        fullobsmodel += repst(2);
        
        fullobsmodel += " + ";
        fullobsmodel += repst(1);
        fullobsmodel += ":";
        fullobsmodel += size(2);
      }
      
      if (age != "none") {
        fullobsmodel += " + ";
        fullobsmodel += age;
        fullobsmodel += ":";
        fullobsmodel += size(1);
        
        fullobsmodel += " + ";
        fullobsmodel += age;
        fullobsmodel += ":";
        fullobsmodel += repst(1);
        
        if (historical) {
          fullobsmodel += " + ";
          fullobsmodel += age;
          fullobsmodel += ":";
          fullobsmodel += size(2);
          
          fullobsmodel += " + ";
          fullobsmodel += age;
          fullobsmodel += ":";
          fullobsmodel += repst(2);
        }
      }
      
      if (indcova(1) != "none") {
        fullobsmodel += " + ";
        fullobsmodel += indcova(1);
        fullobsmodel += ":";
        fullobsmodel += size(1);
        
        fullobsmodel += " + ";
        fullobsmodel += indcova(1);
        fullobsmodel += ":";
        fullobsmodel += repst(1);
        
        if (historical && indcova(2) != "none") {
          fullobsmodel += " + ";
          fullobsmodel += indcova(2);
          fullobsmodel += ":";
          fullobsmodel += size(2);
          
          fullobsmodel += " + ";
          fullobsmodel += indcova(2);
          fullobsmodel += ":";
          fullobsmodel += repst(2);
        }
      }
      if (indcovb(1) != "none") {
        fullobsmodel += " + ";
        fullobsmodel += indcovb(1);
        fullobsmodel += ":";
        fullobsmodel += size(1);
        
        fullobsmodel += " + ";
        fullobsmodel += indcovb(1);
        fullobsmodel += ":";
        fullobsmodel += repst(1);
        
        if (historical && indcovb(2) != "none") {
          fullobsmodel += " + ";
          fullobsmodel += indcovb(2);
          fullobsmodel += ":";
          fullobsmodel += size(2);
          
          fullobsmodel += " + ";
          fullobsmodel += indcovb(2);
          fullobsmodel += ":";
          fullobsmodel += repst(2);
        }
        
        if (indcova(1) != "none") {
          fullobsmodel += " + ";
          fullobsmodel += indcova(1);
          fullobsmodel += ":";
          fullobsmodel += indcovb(1);
          
          if (historical) {
            if (indcova(2) != "none" && indcovb(2) != "none") {
              fullobsmodel += " + ";
              fullobsmodel += indcova(2);
              fullobsmodel += ":";
              fullobsmodel += indcovb(2);
              
              fullobsmodel += " + ";
              fullobsmodel += indcova(1);
              fullobsmodel += ":";
              fullobsmodel += indcovb(2);
              
              fullobsmodel += " + ";
              fullobsmodel += indcova(2);
              fullobsmodel += ":";
              fullobsmodel += indcovb(1);
            } else if (indcova(2) != "none") {
              fullobsmodel += " + ";
              fullobsmodel += indcova(2);
              fullobsmodel += ":";
              fullobsmodel += indcovb(1);
            } else if (indcovb(2) != "none") {
              fullobsmodel += " + ";
              fullobsmodel += indcova(1);
              fullobsmodel += ":";
              fullobsmodel += indcovb(2);
            }
          } 
        }
      }
      
      if (indcovc(1) != "none") {
        fullobsmodel += " + ";
        fullobsmodel += indcovc(1);
        fullobsmodel += ":";
        fullobsmodel += size(1);
        
        fullobsmodel += " + ";
        fullobsmodel += indcovc(1);
        fullobsmodel += ":";
        fullobsmodel += repst(1);
        
        if (historical && indcovc(2) != "none") {
          fullobsmodel += " + ";
          fullobsmodel += indcovc(2);
          fullobsmodel += ":";
          fullobsmodel += size(2);
          
          fullobsmodel += " + ";
          fullobsmodel += indcovc(2);
          fullobsmodel += ":";
          fullobsmodel += repst(2);
        }
        if (indcova(1) != "none") {
          fullobsmodel += " + ";
          fullobsmodel += indcova(1);
          fullobsmodel += ":";
          fullobsmodel += indcovc(1);
          
          if (historical) {
            if (indcova(2) != "none" && indcovc(2) != "none") {
              fullobsmodel += " + ";
              fullobsmodel += indcova(2);
              fullobsmodel += ":";
              fullobsmodel += indcovc(2);
              
              fullobsmodel += " + ";
              fullobsmodel += indcova(1);
              fullobsmodel += ":";
              fullobsmodel += indcovc(2);
              
              fullobsmodel += " + ";
              fullobsmodel += indcova(2);
              fullobsmodel += ":";
              fullobsmodel += indcovc(1);
            } else if (indcova(2) != "none") {
              fullobsmodel += " + ";
              fullobsmodel += indcova(2);
              fullobsmodel += ":";
              fullobsmodel += indcovc(1);
            } else if (indcovc(2) != "none") {
              fullobsmodel += " + ";
              fullobsmodel += indcova(1);
              fullobsmodel += ":";
              fullobsmodel += indcovc(2);
            }
          } 
        }
        if (indcovb(1) != "none") {
          fullobsmodel += " + ";
          fullobsmodel += indcovb(1);
          fullobsmodel += ":";
          fullobsmodel += indcovc(1);
          
          if (historical) {
            if (indcovb(2) != "none" && indcovc(2) != "none") {
              fullobsmodel += " + ";
              fullobsmodel += indcovb(2);
              fullobsmodel += ":";
              fullobsmodel += indcovc(2);
              
              fullobsmodel += " + ";
              fullobsmodel += indcovb(1);
              fullobsmodel += ":";
              fullobsmodel += indcovc(2);
              
              fullobsmodel += " + ";
              fullobsmodel += indcovb(2);
              fullobsmodel += ":";
              fullobsmodel += indcovc(1);
            } else if (indcovb(2) != "none") {
              fullobsmodel += " + ";
              fullobsmodel += indcovb(2);
              fullobsmodel += ":";
              fullobsmodel += indcovc(1);
            } else if (indcovc(2) != "none") {
              fullobsmodel += " + ";
              fullobsmodel += indcovb(1);
              fullobsmodel += ":";
              fullobsmodel += indcovc(2);
            }
          } 
        }
      }
    } else if (suite == "main") {
      if (age != "none" || covcount > 0) fullobsmodel += " + ";
      fullobsmodel += size(1);
      
      if (historical && sizecheck) {
        fullobsmodel += " + ";
        fullobsmodel += size(2);
      }
      
      fullobsmodel += " + ";
      fullobsmodel += repst(1);
      
      if (historical && repstcheck) {
        fullobsmodel += " + ";
        fullobsmodel += repst(2);
      }
      
    } else if (suite == "size") {
      if (age != "none" || covcount > 0) fullobsmodel += " + ";
      fullobsmodel += size(1);
      
      if (historical && sizecheck) {
        fullobsmodel += " + ";
        fullobsmodel += size(2);
        
        fullobsmodel += " + ";
        fullobsmodel += size(1);
        fullobsmodel += ":";
        fullobsmodel += size(2);
      }
      
    } else if (suite == "rep") {
      if (age != "none" || covcount > 0) fullobsmodel += " + ";
      fullobsmodel += repst(1);
      
      if (historical && repstcheck) {
        fullobsmodel += " + ";
        fullobsmodel += repst(2);
        
        fullobsmodel += " + ";
        fullobsmodel += repst(1);
        fullobsmodel += ":";
        fullobsmodel += repst(2);
      }
      
    } else {
      if (age == "none" && covcount == 0) fullobsmodel += "1";
    }
    
    fullobsmodel += fixedtackon;
    fullobsmodel += randomtackon;
    
  } else {
    fullobsmodel = "none";
  }
  
  if (sizecheck) {
    fullsizemodel = size(0);
    fullsizemodel += " ~ ";
    sizesuffix = "";
    
    if (!nojuvs) {
      juvsizemodel = size(0);
      juvsizemodel += " ~ ";
      
      if (suite == "full" || suite == "main" || suite == "size") {
        
        if (juvsize) {
          juvsizemodel += size(1);
          jsizesuffix = size(1);
        } else {
          juvsizemodel += "1";
        }
        
      } else if (approach == "mixed") {
        juvsizemodel += "1"; // This should be altered, since it will probably result in 1 for glm with year or other factor
      }
      
      juvsizemodel += fixedtackon;
      juvsizemodel += randomtackon;
      jsizesuffix += fixedtackon;
      
      if (approach == "glm" && size0 && sizedist == "poisson") {
        juvsizemodel += " | ";
        juvsizemodel += jsizesuffix;
      }
    }
    
    if (age != "none") {
      fullsizemodel += age;
      sizesuffix += age;
    }
    if (indcova(1) != "none") {
      if (age != "none") {
        fullsizemodel += " + ";
        sizesuffix += " + ";
      }
      fullsizemodel += indcova(1);
      sizesuffix += indcova(1);
      
      if (historical) {
        fullsizemodel += " + ";
        fullsizemodel += indcova(2);
        
        sizesuffix += " + ";
        sizesuffix += indcova(2);
      }
    }
    if (indcovb(1) != "none") {
      if (age != "none" || indcova(1) != "none") {
        fullsizemodel += " + ";
        sizesuffix += " + ";
      }
      fullsizemodel += indcovb(1);
      sizesuffix += indcovb(1);
      
      if (historical) {
        fullsizemodel += " + ";
        fullsizemodel += indcovb(2);
        
        sizesuffix += " + ";
        sizesuffix += indcovb(2);
      }
    }
    if (indcovc(1) != "none") {
      if (age != "none" || covcount > 1) {
        fullsizemodel += " + ";
        sizesuffix += " + ";
      }
      fullsizemodel += indcovc(1);
      sizesuffix += indcovc(1);
      
      if (historical) {
        fullsizemodel += " + ";
        fullsizemodel += indcovc(2);
        
        sizesuffix += " + ";
        sizesuffix += indcovc(2);
      }
    }
    
    if (suite == "full") {
      if (age != "none" || covcount > 0) {
        fullsizemodel += " + ";
        sizesuffix += " + ";
      }
      fullsizemodel += size(1);
      fullsizemodel += " + ";
      fullsizemodel += repst(1);
      
      sizesuffix += size(1);
      sizesuffix += " + ";
      sizesuffix += repst(1);
      
      if (historical) {
        fullsizemodel += " + ";
        fullsizemodel += size(2);
        fullsizemodel += " + ";
        fullsizemodel += repst(2);
        
        sizesuffix += " + ";
        sizesuffix += size(2);
        sizesuffix += " + ";
        sizesuffix += repst(2);
      }
      
      fullsizemodel += " + ";
      fullsizemodel += size(1);
      fullsizemodel += ":";
      fullsizemodel += repst(1);
      
      sizesuffix += " + ";
      sizesuffix += size(1);
      sizesuffix += ":";
      sizesuffix += repst(1);
      
      if (historical) {
        fullsizemodel += " + ";
        fullsizemodel += size(1);
        fullsizemodel += ":";
        fullsizemodel += size(2);
        
        fullsizemodel += " + ";
        fullsizemodel += repst(1);
        fullsizemodel += ":";
        fullsizemodel += repst(2);
        
        fullsizemodel += " + ";
        fullsizemodel += size(1);
        fullsizemodel += ":";
        fullsizemodel += repst(2);
        
        fullsizemodel += " + ";
        fullsizemodel += repst(1);
        fullsizemodel += ":";
        fullsizemodel += size(2);
        
        sizesuffix += " + ";
        sizesuffix += size(1);
        sizesuffix += ":";
        sizesuffix += size(2);
        
        sizesuffix += " + ";
        sizesuffix += repst(1);
        sizesuffix += ":";
        sizesuffix += repst(2);
        
        sizesuffix += " + ";
        sizesuffix += size(1);
        sizesuffix += ":";
        sizesuffix += repst(2);
        
        sizesuffix += " + ";
        sizesuffix += repst(1);
        sizesuffix += ":";
        sizesuffix += size(2);
      }
      
      if (age != "none") {
        fullsizemodel += " + ";
        fullsizemodel += age;
        fullsizemodel += ":";
        fullsizemodel += size(1);
        
        fullsizemodel += " + ";
        fullsizemodel += age;
        fullsizemodel += ":";
        fullsizemodel += repst(1);
        
        sizesuffix += " + ";
        sizesuffix += age;
        sizesuffix += ":";
        sizesuffix += size(1);
        
        sizesuffix += " + ";
        sizesuffix += age;
        sizesuffix += ":";
        sizesuffix += repst(1);
        
        if (historical) {
          fullsizemodel += " + ";
          fullsizemodel += age;
          fullsizemodel += ":";
          fullsizemodel += size(2);
          
          fullsizemodel += " + ";
          fullsizemodel += age;
          fullsizemodel += ":";
          fullsizemodel += repst(2);
          
          sizesuffix += " + ";
          sizesuffix += age;
          sizesuffix += ":";
          sizesuffix += size(2);
          
          sizesuffix += " + ";
          sizesuffix += age;
          sizesuffix += ":";
          sizesuffix += repst(2);
        }
      }
      
      if (indcova(1) != "none") {
        fullsizemodel += " + ";
        fullsizemodel += indcova(1);
        fullsizemodel += ":";
        fullsizemodel += size(1);
        
        fullsizemodel += " + ";
        fullsizemodel += indcova(1);
        fullsizemodel += ":";
        fullsizemodel += repst(1);
        
        sizesuffix += " + ";
        sizesuffix += indcova(1);
        sizesuffix += ":";
        sizesuffix += size(1);
        
        sizesuffix += " + ";
        sizesuffix += indcova(1);
        sizesuffix += ":";
        sizesuffix += repst(1);
        
        if (historical && indcova(2) != "none") {
          fullsizemodel += " + ";
          fullsizemodel += indcova(2);
          fullsizemodel += ":";
          fullsizemodel += size(2);
          
          fullsizemodel += " + ";
          fullsizemodel += indcova(2);
          fullsizemodel += ":";
          fullsizemodel += repst(2);
          
          sizesuffix += " + ";
          sizesuffix += indcova(2);
          sizesuffix += ":";
          sizesuffix += size(2);
          
          sizesuffix += " + ";
          sizesuffix += indcova(2);
          sizesuffix += ":";
          sizesuffix += repst(2);
        }
      }
      if (indcovb(1) != "none") {
        fullsizemodel += " + ";
        fullsizemodel += indcovb(1);
        fullsizemodel += ":";
        fullsizemodel += size(1);
        
        fullsizemodel += " + ";
        fullsizemodel += indcovb(1);
        fullsizemodel += ":";
        fullsizemodel += repst(1);
        
        sizesuffix += " + ";
        sizesuffix += indcovb(1);
        sizesuffix += ":";
        sizesuffix += size(1);
        
        sizesuffix += " + ";
        sizesuffix += indcovb(1);
        sizesuffix += ":";
        sizesuffix += repst(1);
        
        if (historical && indcovb(2) != "none") {
          fullsizemodel += " + ";
          fullsizemodel += indcovb(2);
          fullsizemodel += ":";
          fullsizemodel += size(2);
          
          fullsizemodel += " + ";
          fullsizemodel += indcovb(2);
          fullsizemodel += ":";
          fullsizemodel += repst(2);
          
          sizesuffix += " + ";
          sizesuffix += indcovb(2);
          sizesuffix += ":";
          sizesuffix += size(2);
          
          sizesuffix += " + ";
          sizesuffix += indcovb(2);
          sizesuffix += ":";
          sizesuffix += repst(2);
        }
        
        if (indcova(1) != "none") {
          fullsizemodel += " + ";
          fullsizemodel += indcova(1);
          fullsizemodel += ":";
          fullsizemodel += indcovb(1);
          
          sizesuffix += " + ";
          sizesuffix += indcova(1);
          sizesuffix += ":";
          sizesuffix += indcovb(1);
          
          if (historical) {
            if (indcova(2) != "none" && indcovb(2) != "none") {
              fullsizemodel += " + ";
              fullsizemodel += indcova(2);
              fullsizemodel += ":";
              fullsizemodel += indcovb(2);
              
              fullsizemodel += " + ";
              fullsizemodel += indcova(1);
              fullsizemodel += ":";
              fullsizemodel += indcovb(2);
              
              fullsizemodel += " + ";
              fullsizemodel += indcova(2);
              fullsizemodel += ":";
              fullsizemodel += indcovb(1);
              
              sizesuffix += " + ";
              sizesuffix += indcova(2);
              sizesuffix += ":";
              sizesuffix += indcovb(2);
              
              sizesuffix += " + ";
              sizesuffix += indcova(1);
              sizesuffix += ":";
              sizesuffix += indcovb(2);
              
              sizesuffix += " + ";
              sizesuffix += indcova(2);
              sizesuffix += ":";
              sizesuffix += indcovb(1);
            } else if (indcova(2) != "none") {
              fullsizemodel += " + ";
              fullsizemodel += indcova(2);
              fullsizemodel += ":";
              fullsizemodel += indcovb(1);
              
              sizesuffix += " + ";
              sizesuffix += indcova(2);
              sizesuffix += ":";
              sizesuffix += indcovb(1);
            } else if (indcovb(2) != "none") {
              fullsizemodel += " + ";
              fullsizemodel += indcova(1);
              fullsizemodel += ":";
              fullsizemodel += indcovb(2);
              
              sizesuffix += " + ";
              sizesuffix += indcova(1);
              sizesuffix += ":";
              sizesuffix += indcovb(2);
            }
          } 
        }
      }
      
      if (indcovc(1) != "none") {
        fullsizemodel += " + ";
        fullsizemodel += indcovc(1);
        fullsizemodel += ":";
        fullsizemodel += size(1);
        
        fullsizemodel += " + ";
        fullsizemodel += indcovc(1);
        fullsizemodel += ":";
        fullsizemodel += repst(1);
        
        sizesuffix += " + ";
        sizesuffix += indcovc(1);
        sizesuffix += ":";
        sizesuffix += size(1);
        
        sizesuffix += " + ";
        sizesuffix += indcovc(1);
        sizesuffix += ":";
        sizesuffix += repst(1);
        
        if (historical && indcovc(2) != "none") {
          fullsizemodel += " + ";
          fullsizemodel += indcovc(2);
          fullsizemodel += ":";
          fullsizemodel += size(2);
          
          fullsizemodel += " + ";
          fullsizemodel += indcovc(2);
          fullsizemodel += ":";
          fullsizemodel += repst(2);
          
          sizesuffix += " + ";
          sizesuffix += indcovc(2);
          sizesuffix += ":";
          sizesuffix += size(2);
          
          sizesuffix += " + ";
          sizesuffix += indcovc(2);
          sizesuffix += ":";
          sizesuffix += repst(2);
        }
        if (indcova(1) != "none") {
          fullsizemodel += " + ";
          fullsizemodel += indcova(1);
          fullsizemodel += ":";
          fullsizemodel += indcovc(1);
          
          sizesuffix += " + ";
          sizesuffix += indcova(1);
          sizesuffix += ":";
          sizesuffix += indcovc(1);
          
          if (historical) {
            if (indcova(2) != "none" && indcovc(2) != "none") {
              fullsizemodel += " + ";
              fullsizemodel += indcova(2);
              fullsizemodel += ":";
              fullsizemodel += indcovc(2);
              
              fullsizemodel += " + ";
              fullsizemodel += indcova(1);
              fullsizemodel += ":";
              fullsizemodel += indcovc(2);
              
              fullsizemodel += " + ";
              fullsizemodel += indcova(2);
              fullsizemodel += ":";
              fullsizemodel += indcovc(1);
              
              sizesuffix += " + ";
              sizesuffix += indcova(2);
              sizesuffix += ":";
              sizesuffix += indcovc(2);
              
              sizesuffix += " + ";
              sizesuffix += indcova(1);
              sizesuffix += ":";
              sizesuffix += indcovc(2);
              
              sizesuffix += " + ";
              sizesuffix += indcova(2);
              sizesuffix += ":";
              sizesuffix += indcovc(1);
            } else if (indcova(2) != "none") {
              fullsizemodel += " + ";
              fullsizemodel += indcova(2);
              fullsizemodel += ":";
              fullsizemodel += indcovc(1);
              
              sizesuffix += " + ";
              sizesuffix += indcova(2);
              sizesuffix += ":";
              sizesuffix += indcovc(1);
            } else if (indcovc(2) != "none") {
              fullsizemodel += " + ";
              fullsizemodel += indcova(1);
              fullsizemodel += ":";
              fullsizemodel += indcovc(2);
              
              sizesuffix += " + ";
              sizesuffix += indcova(1);
              sizesuffix += ":";
              sizesuffix += indcovc(2);
            }
          } 
        }
        if (indcovb(1) != "none") {
          fullsizemodel += " + ";
          fullsizemodel += indcovb(1);
          fullsizemodel += ":";
          fullsizemodel += indcovc(1);
          
          sizesuffix += " + ";
          sizesuffix += indcovb(1);
          sizesuffix += ":";
          sizesuffix += indcovc(1);
          
          if (historical) {
            if (indcovb(2) != "none" && indcovc(2) != "none") {
              fullsizemodel += " + ";
              fullsizemodel += indcovb(2);
              fullsizemodel += ":";
              fullsizemodel += indcovc(2);
              
              fullsizemodel += " + ";
              fullsizemodel += indcovb(1);
              fullsizemodel += ":";
              fullsizemodel += indcovc(2);
              
              fullsizemodel += " + ";
              fullsizemodel += indcovb(2);
              fullsizemodel += ":";
              fullsizemodel += indcovc(1);
              
              sizesuffix += " + ";
              sizesuffix += indcovb(2);
              sizesuffix += ":";
              sizesuffix += indcovc(2);
              
              sizesuffix += " + ";
              sizesuffix += indcovb(1);
              sizesuffix += ":";
              sizesuffix += indcovc(2);
              
              sizesuffix += " + ";
              sizesuffix += indcovb(2);
              sizesuffix += ":";
              sizesuffix += indcovc(1);
            } else if (indcovb(2) != "none") {
              fullsizemodel += " + ";
              fullsizemodel += indcovb(2);
              fullsizemodel += ":";
              fullsizemodel += indcovc(1);
              
              sizesuffix += " + ";
              sizesuffix += indcovb(2);
              sizesuffix += ":";
              sizesuffix += indcovc(1);
            } else if (indcovc(2) != "none") {
              fullsizemodel += " + ";
              fullsizemodel += indcovb(1);
              fullsizemodel += ":";
              fullsizemodel += indcovc(2);
              
              sizesuffix += " + ";
              sizesuffix += indcovb(1);
              sizesuffix += ":";
              sizesuffix += indcovc(2);
            }
          } 
        }
      }
    } else if (suite == "main") {
      if (age != "none" || covcount > 0) {
        fullsizemodel += " + ";
        sizesuffix += " + ";
      }
      
      fullsizemodel += size(1);
      sizesuffix += size(1);
      
      if (historical && sizecheck) {
        fullsizemodel += " + ";
        fullsizemodel += size(2);
        
        sizesuffix += " + ";
        sizesuffix += size(2);
      }
      
      fullsizemodel += " + ";
      fullsizemodel += repst(1);
      
      sizesuffix += " + ";
      sizesuffix += repst(1);
      
      if (historical && repstcheck) {
        fullsizemodel += " + ";
        fullsizemodel += repst(2);
        
        sizesuffix += " + ";
        sizesuffix += repst(2);
      }
      
    } else if (suite == "size") {
      if (age != "none" || covcount > 0) {
        fullsizemodel += " + ";
        sizesuffix += " + ";
      }
      
      fullsizemodel += size(1);
      sizesuffix += size(1);
      
      if (historical && sizecheck) {
        fullsizemodel += " + ";
        fullsizemodel += size(2);
        
        sizesuffix += " + ";
        sizesuffix += size(2);
        
        fullsizemodel += " + ";
        fullsizemodel += size(1);
        fullsizemodel += ":";
        fullsizemodel += size(2);
        
        sizesuffix += " + ";
        sizesuffix += size(1);
        sizesuffix += ":";
        sizesuffix += size(2);
      }
      
    } else if (suite == "rep") {
      if (age != "none" || covcount > 0) {
        fullsizemodel += " + ";
        sizesuffix += " + ";
      }
      
      fullsizemodel += repst(1);
      sizesuffix += repst(1);
      
      if (historical && repstcheck) {
        fullsizemodel += " + ";
        fullsizemodel += repst(2);
        
        sizesuffix += " + ";
        sizesuffix += repst(2);
        
        fullsizemodel += " + ";
        fullsizemodel += repst(1);
        fullsizemodel += ":";
        fullsizemodel += repst(2);
        
        sizesuffix += " + ";
        sizesuffix += repst(1);
        sizesuffix += ":";
        sizesuffix += repst(2);
      }
      
    } else {
      if (age == "none" && covcount == 1) fullsizemodel += "1";
    }
    
    fullsizemodel += fixedtackon;
    sizesuffix += fixedtackon;
    fullsizemodel += randomtackon;
    
    if (size0 && suite != "const" && approach == "glm") {
      if (sizedist == "poisson") {
        fullsizemodel += " | ";
        fullsizemodel += sizesuffix;
      }
    }
    
  } else {
    fullsizemodel = "none";
  }
  
  if (repstcheck) {
    fullrepstmodel = repst(0);
    fullrepstmodel += " ~ ";
    
    if (!nojuvs) {
      juvrepstmodel = repst(0);
      juvrepstmodel += " ~ ";
      
      if (suite == "full" || suite == "main" || suite == "size") {
        if (juvsize) {
          juvrepstmodel += size(1);
        } else juvrepstmodel += "1";
      } else juvrepstmodel += "1";
      
      juvrepstmodel += fixedtackon;
      juvrepstmodel += randomtackon;
    }
    
    if (age != "none") {
      fullrepstmodel += age;
    }
    if (indcova(1) != "none") {
      if (age != "none") fullrepstmodel += " + ";
      fullrepstmodel += indcova(1);
      
      if (historical) {
        fullrepstmodel += " + ";
        fullrepstmodel += indcova(2);
      }
    }
    if (indcovb(1) != "none") {
      if (age != "none" || indcova(1) != "none") fullrepstmodel += " + ";
      fullrepstmodel += indcovb(1);
      
      if (historical) {
        fullrepstmodel += " + ";
        fullrepstmodel += indcovb(2);
      }
    }
    if (indcovc(1) != "none") {
      if (age != "none" || covcount > 1) fullrepstmodel += " + ";
      fullrepstmodel += indcovc(1);
      
      if (historical) {
        fullrepstmodel += " + ";
        fullrepstmodel += indcovc(2);
      }
    }
    
    if (suite == "full") {
      if (age != "none" || covcount > 0) fullrepstmodel += " + ";
      fullrepstmodel += size(1);
      fullrepstmodel += " + ";
      fullrepstmodel += repst(1);
      
      if (historical) {
        fullrepstmodel += " + ";
        fullrepstmodel += size(2);
        fullrepstmodel += " + ";
        fullrepstmodel += repst(2);
      }
      
      fullrepstmodel += " + ";
      fullrepstmodel += size(1);
      fullrepstmodel += ":";
      fullrepstmodel += repst(1);
      
      if (historical) {
        fullrepstmodel += " + ";
        fullrepstmodel += size(1);
        fullrepstmodel += ":";
        fullrepstmodel += size(2);
        
        fullrepstmodel += " + ";
        fullrepstmodel += repst(1);
        fullrepstmodel += ":";
        fullrepstmodel += repst(2);
        
        fullrepstmodel += " + ";
        fullrepstmodel += size(1);
        fullrepstmodel += ":";
        fullrepstmodel += repst(2);
        
        fullrepstmodel += " + ";
        fullrepstmodel += repst(1);
        fullrepstmodel += ":";
        fullrepstmodel += size(2);
      }
      
      if (age != "none") {
        fullrepstmodel += " + ";
        fullrepstmodel += age;
        fullrepstmodel += ":";
        fullrepstmodel += size(1);
        
        fullrepstmodel += " + ";
        fullrepstmodel += age;
        fullrepstmodel += ":";
        fullrepstmodel += repst(1);
        
        if (historical) {
          fullrepstmodel += " + ";
          fullrepstmodel += age;
          fullrepstmodel += ":";
          fullrepstmodel += size(2);
          
          fullrepstmodel += " + ";
          fullrepstmodel += age;
          fullrepstmodel += ":";
          fullrepstmodel += repst(2);
        }
      }
      
      if (indcova(1) != "none") {
        fullrepstmodel += " + ";
        fullrepstmodel += indcova(1);
        fullrepstmodel += ":";
        fullrepstmodel += size(1);
        
        fullrepstmodel += " + ";
        fullrepstmodel += indcova(1);
        fullrepstmodel += ":";
        fullrepstmodel += repst(1);
        
        if (historical && indcova(2) != "none") {
          fullrepstmodel += " + ";
          fullrepstmodel += indcova(2);
          fullrepstmodel += ":";
          fullrepstmodel += size(2);
          
          fullrepstmodel += " + ";
          fullrepstmodel += indcova(2);
          fullrepstmodel += ":";
          fullrepstmodel += repst(2);
        }
      }
      if (indcovb(1) != "none") {
        fullrepstmodel += " + ";
        fullrepstmodel += indcovb(1);
        fullrepstmodel += ":";
        fullrepstmodel += size(1);
        
        fullrepstmodel += " + ";
        fullrepstmodel += indcovb(1);
        fullrepstmodel += ":";
        fullrepstmodel += repst(1);
        
        if (historical && indcovb(2) != "none") {
          fullrepstmodel += " + ";
          fullrepstmodel += indcovb(2);
          fullrepstmodel += ":";
          fullrepstmodel += size(2);
          
          fullrepstmodel += " + ";
          fullrepstmodel += indcovb(2);
          fullrepstmodel += ":";
          fullrepstmodel += repst(2);
        }
        
        if (indcova(1) != "none") {
          fullrepstmodel += " + ";
          fullrepstmodel += indcova(1);
          fullrepstmodel += ":";
          fullrepstmodel += indcovb(1);
          
          if (historical) {
            if (indcova(2) != "none" && indcovb(2) != "none") {
              fullrepstmodel += " + ";
              fullrepstmodel += indcova(2);
              fullrepstmodel += ":";
              fullrepstmodel += indcovb(2);
              
              fullrepstmodel += " + ";
              fullrepstmodel += indcova(1);
              fullrepstmodel += ":";
              fullrepstmodel += indcovb(2);
              
              fullrepstmodel += " + ";
              fullrepstmodel += indcova(2);
              fullrepstmodel += ":";
              fullrepstmodel += indcovb(1);
            } else if (indcova(2) != "none") {
              fullrepstmodel += " + ";
              fullrepstmodel += indcova(2);
              fullrepstmodel += ":";
              fullrepstmodel += indcovb(1);
            } else if (indcovb(2) != "none") {
              fullrepstmodel += " + ";
              fullrepstmodel += indcova(1);
              fullrepstmodel += ":";
              fullrepstmodel += indcovb(2);
            }
          } 
        }
      }
      
      if (indcovc(1) != "none") {
        fullrepstmodel += " + ";
        fullrepstmodel += indcovc(1);
        fullrepstmodel += ":";
        fullrepstmodel += size(1);
        
        fullrepstmodel += " + ";
        fullrepstmodel += indcovc(1);
        fullrepstmodel += ":";
        fullrepstmodel += repst(1);
        
        if (historical && indcovc(2) != "none") {
          fullrepstmodel += " + ";
          fullrepstmodel += indcovc(2);
          fullrepstmodel += ":";
          fullrepstmodel += size(2);
          
          fullrepstmodel += " + ";
          fullrepstmodel += indcovc(2);
          fullrepstmodel += ":";
          fullrepstmodel += repst(2);
        }
        if (indcova(1) != "none") {
          fullrepstmodel += " + ";
          fullrepstmodel += indcova(1);
          fullrepstmodel += ":";
          fullrepstmodel += indcovc(1);
          
          if (historical) {
            if (indcova(2) != "none" && indcovc(2) != "none") {
              fullrepstmodel += " + ";
              fullrepstmodel += indcova(2);
              fullrepstmodel += ":";
              fullrepstmodel += indcovc(2);
              
              fullrepstmodel += " + ";
              fullrepstmodel += indcova(1);
              fullrepstmodel += ":";
              fullrepstmodel += indcovc(2);
              
              fullrepstmodel += " + ";
              fullrepstmodel += indcova(2);
              fullrepstmodel += ":";
              fullrepstmodel += indcovc(1);
            } else if (indcova(2) != "none") {
              fullrepstmodel += " + ";
              fullrepstmodel += indcova(2);
              fullrepstmodel += ":";
              fullrepstmodel += indcovc(1);
            } else if (indcovc(2) != "none") {
              fullrepstmodel += " + ";
              fullrepstmodel += indcova(1);
              fullrepstmodel += ":";
              fullrepstmodel += indcovc(2);
            }
          } 
        }
        if (indcovb(1) != "none") {
          fullrepstmodel += " + ";
          fullrepstmodel += indcovb(1);
          fullrepstmodel += ":";
          fullrepstmodel += indcovc(1);
          
          if (historical) {
            if (indcovb(2) != "none" && indcovc(2) != "none") {
              fullrepstmodel += " + ";
              fullrepstmodel += indcovb(2);
              fullrepstmodel += ":";
              fullrepstmodel += indcovc(2);
              
              fullrepstmodel += " + ";
              fullrepstmodel += indcovb(1);
              fullrepstmodel += ":";
              fullrepstmodel += indcovc(2);
              
              fullrepstmodel += " + ";
              fullrepstmodel += indcovb(2);
              fullrepstmodel += ":";
              fullrepstmodel += indcovc(1);
            } else if (indcovb(2) != "none") {
              fullrepstmodel += " + ";
              fullrepstmodel += indcovb(2);
              fullrepstmodel += ":";
              fullrepstmodel += indcovc(1);
            } else if (indcovc(2) != "none") {
              fullrepstmodel += " + ";
              fullrepstmodel += indcovb(1);
              fullrepstmodel += ":";
              fullrepstmodel += indcovc(2);
            }
          } 
        }
        
      }
    } else if (suite == "main") {
      if (age != "none" || covcount > 0) fullrepstmodel += " + ";
      fullrepstmodel += size(1);
      
      if (historical && sizecheck) {
        fullrepstmodel += " + ";
        fullrepstmodel += size(2);
      }
      
      fullrepstmodel += " + ";
      fullrepstmodel += repst(1);
      
      if (historical && repstcheck) {
        fullrepstmodel += " + ";
        fullrepstmodel += repst(2);
      }
      
    } else if (suite == "size") {
      if (age != "none" || covcount > 0) fullrepstmodel += " + ";
      fullrepstmodel += size(1);
      
      if (historical && sizecheck) {
        fullrepstmodel += " + ";
        fullrepstmodel += size(2);
        
        fullrepstmodel += " + ";
        fullrepstmodel += size(1);
        fullrepstmodel += ":";
        fullrepstmodel += size(2);
      }
      
    } else if (suite == "rep") {
      if (age != "none" || covcount > 0) fullrepstmodel += " + ";
      fullrepstmodel += repst(1);
      
      if (historical && repstcheck) {
        fullrepstmodel += " + ";
        fullrepstmodel += repst(2);
        
        fullrepstmodel += " + ";
        fullrepstmodel += repst(1);
        fullrepstmodel += ":";
        fullrepstmodel += repst(2);
      }
      
    } else {
      if (age == "none" && covcount == 0) fullrepstmodel += "1";
    }
    
    fullrepstmodel += fixedtackon;
    fullrepstmodel += randomtackon;
    
  } else {
    fullrepstmodel = "none";
  }
  
  if (feccheck) {
    if (fectime == 3) {
      fullfecmodel = fec(0);
    } else {
      fullfecmodel = fec(1);
    }
    fullfecmodel += " ~ ";
    fecsuffix = "";
    
    if (age != "none") {
      fullfecmodel += age;
      fecsuffix += age;
    }
    if (indcova(1) != "none") {
      if (age != "none") {
        
        fullfecmodel += " + ";
        fecsuffix += " + ";
      }
      fullfecmodel += indcova(1);
      fecsuffix += indcova(1);
      
      if (historical) {
        fullfecmodel += " + ";
        fullfecmodel += indcova(2);
        
        fecsuffix += " + ";
        fecsuffix += indcova(2);
      }
    }
    if (indcovb(1) != "none") {
      if (age != "none" || indcova(1) != "none") {
        
        fullfecmodel += " + ";
        fecsuffix += " + ";
      }
      fullfecmodel += indcovb(1);
      fecsuffix += indcovb(1);
      
      if (historical) {
        fullfecmodel += " + ";
        fullfecmodel += indcovb(2);
        
        fecsuffix += " + ";
        fecsuffix += indcovb(2);
      }
    }
    if (indcovc(1) != "none") {
      if (age != "none" || covcount > 1) {
        
        fullfecmodel += " + ";
        fecsuffix += " + ";
      }
      fullfecmodel += indcovc(1);
      fecsuffix += indcovc(1);
      
      if (historical) {
        fullfecmodel += " + ";
        fullfecmodel += indcovc(2);
        
        fecsuffix += " + ";
        fecsuffix += indcovc(2);
      }
    }
    
    if (suite == "full") {
      if (age != "none" || covcount > 0) {
        
        fullfecmodel += " + ";
        fecsuffix += " + ";
      }
      
      fullfecmodel += size(1);
      fecsuffix += size(1);
      
      if (fectime == 3) {
        fullfecmodel += " + ";
        fullfecmodel += repst(1);
        
        fecsuffix += " + ";
        fecsuffix += repst(1);
      }
      
      if (historical) {
        fullfecmodel += " + ";
        fullfecmodel += size(2);
        fullfecmodel += " + ";
        fullfecmodel += repst(2);
        
        fecsuffix += " + ";
        fecsuffix += size(2);
        fecsuffix += " + ";
        fecsuffix += repst(2);
      }
      
      if (fectime == 3) {
        fullfecmodel += " + ";
        fullfecmodel += size(1);
        fullfecmodel += ":";
        fullfecmodel += repst(1);
        
        fecsuffix += " + ";
        fecsuffix += size(1);
        fecsuffix += ":";
        fecsuffix += repst(1);
      }
      
      if (historical) {
        fullfecmodel += " + ";
        fullfecmodel += size(1);
        fullfecmodel += ":";
        fullfecmodel += size(2);
        
        fecsuffix += " + ";
        fecsuffix += size(1);
        fecsuffix += ":";
        fecsuffix += size(2);
        
        if (fectime == 3) {
          fullfecmodel += " + ";
          fullfecmodel += repst(1);
          fullfecmodel += ":";
          fullfecmodel += repst(2);
          
          fecsuffix += " + ";
          fecsuffix += repst(1);
          fecsuffix += ":";
          fecsuffix += repst(2);
        }
        
        fullfecmodel += " + ";
        fullfecmodel += size(1);
        fullfecmodel += ":";
        fullfecmodel += repst(2);
        
        fecsuffix += " + ";
        fecsuffix += size(1);
        fecsuffix += ":";
        fecsuffix += repst(2);
        
        if (fectime == 3) {
          fullfecmodel += " + ";
          fullfecmodel += size(2);
          fullfecmodel += ":";
          fullfecmodel += repst(1);
          
          fecsuffix += " + ";
          fecsuffix += size(2);
          fecsuffix += ":";
          fecsuffix += repst(1);
        }
      }
      
      if (age != "none") {
        fullfecmodel += " + ";
        fullfecmodel += age;
        fullfecmodel += ":";
        fullfecmodel += size(1);
        
        fecsuffix += " + ";
        fecsuffix += age;
        fecsuffix += ":";
        fecsuffix += size(1);
        
        if (fectime == 3) {
          fullfecmodel += " + ";
          fullfecmodel += age;
          fullfecmodel += ":";
          fullfecmodel += repst(1);
          
          fecsuffix += " + ";
          fecsuffix += age;
          fecsuffix += ":";
          fecsuffix += repst(1);
        }
        
        if (historical) {
          fullfecmodel += " + ";
          fullfecmodel += age;
          fullfecmodel += ":";
          fullfecmodel += size(2);
          
          fullfecmodel += " + ";
          fullfecmodel += age;
          fullfecmodel += ":";
          fullfecmodel += repst(2);
          
          fecsuffix += " + ";
          fecsuffix += age;
          fecsuffix += ":";
          fecsuffix += size(2);
          
          fecsuffix += " + ";
          fecsuffix += age;
          fecsuffix += ":";
          fecsuffix += repst(2);
        }
      }
      
      if (indcova(1) != "none") {
        fullfecmodel += " + ";
        fullfecmodel += indcova(1);
        fullfecmodel += ":";
        fullfecmodel += size(1);
        
        fecsuffix += " + ";
        fecsuffix += indcova(1);
        fecsuffix += ":";
        fecsuffix += size(1);
        
        if (fectime == 3) {
          fullfecmodel += " + ";
          fullfecmodel += indcova(1);
          fullfecmodel += ":";
          fullfecmodel += repst(1);
          
          fecsuffix += " + ";
          fecsuffix += indcova(1);
          fecsuffix += ":";
          fecsuffix += repst(1);
        }
        
        if (historical && indcova(2) != "none") {
          fullfecmodel += " + ";
          fullfecmodel += indcova(2);
          fullfecmodel += ":";
          fullfecmodel += size(2);
          
          fullfecmodel += " + ";
          fullfecmodel += indcova(2);
          fullfecmodel += ":";
          fullfecmodel += repst(2);
          
          fecsuffix += " + ";
          fecsuffix += indcova(2);
          fecsuffix += ":";
          fecsuffix += size(2);
          
          fecsuffix += " + ";
          fecsuffix += indcova(2);
          fecsuffix += ":";
          fecsuffix += repst(2);
        }
      }
      
      if (indcovb(1) != "none") {
        fullfecmodel += " + ";
        fullfecmodel += indcovb(1);
        fullfecmodel += ":";
        fullfecmodel += size(1);
        
        fecsuffix += " + ";
        fecsuffix += indcovb(1);
        fecsuffix += ":";
        fecsuffix += size(1);
        
        if (fectime == 3) {
          fullfecmodel += " + ";
          fullfecmodel += indcovb(1);
          fullfecmodel += ":";
          fullfecmodel += repst(1);
          
          fecsuffix += " + ";
          fecsuffix += indcovb(1);
          fecsuffix += ":";
          fecsuffix += repst(1);
        }
        
        if (historical && indcovb(2) != "none") {
          fullfecmodel += " + ";
          fullfecmodel += indcovb(2);
          fullfecmodel += ":";
          fullfecmodel += size(2);
          
          fullfecmodel += " + ";
          fullfecmodel += indcovb(2);
          fullfecmodel += ":";
          fullfecmodel += repst(2);
          
          fecsuffix += " + ";
          fecsuffix += indcovb(2);
          fecsuffix += ":";
          fecsuffix += size(2);
          
          fecsuffix += " + ";
          fecsuffix += indcovb(2);
          fecsuffix += ":";
          fecsuffix += repst(2);
        }
        if (indcova(1) != "none") {
          fullfecmodel += " + ";
          fullfecmodel += indcova(1);
          fullfecmodel += ":";
          fullfecmodel += indcovb(1);
          
          fecsuffix += " + ";
          fecsuffix += indcova(1);
          fecsuffix += ":";
          fecsuffix += indcovb(1);
          
          if (historical) {
            if (indcova(2) != "none" && indcovb(2) != "none") {
              fullfecmodel += " + ";
              fullfecmodel += indcova(2);
              fullfecmodel += ":";
              fullfecmodel += indcovb(2);
              
              fullfecmodel += " + ";
              fullfecmodel += indcova(1);
              fullfecmodel += ":";
              fullfecmodel += indcovb(2);
              
              fullfecmodel += " + ";
              fullfecmodel += indcova(2);
              fullfecmodel += ":";
              fullfecmodel += indcovb(1);
              
              fecsuffix += " + ";
              fecsuffix += indcova(2);
              fecsuffix += ":";
              fecsuffix += indcovb(2);
              
              fecsuffix += " + ";
              fecsuffix += indcova(1);
              fecsuffix += ":";
              fecsuffix += indcovb(2);
              
              fecsuffix += " + ";
              fecsuffix += indcova(2);
              fecsuffix += ":";
              fecsuffix += indcovb(1);
            } else if (indcova(2) != "none") {
              fullfecmodel += " + ";
              fullfecmodel += indcova(2);
              fullfecmodel += ":";
              fullfecmodel += indcovb(1);
              
              fecsuffix += " + ";
              fecsuffix += indcova(2);
              fecsuffix += ":";
              fecsuffix += indcovb(1);
            } else if (indcovb(2) != "none") {
              fullfecmodel += " + ";
              fullfecmodel += indcova(1);
              fullfecmodel += ":";
              fullfecmodel += indcovb(2);
              
              fecsuffix += " + ";
              fecsuffix += indcova(1);
              fecsuffix += ":";
              fecsuffix += indcovb(2);
            }
          } 
        }
      }
      
      if (indcovc(1) != "none") {
        fullfecmodel += " + ";
        fullfecmodel += indcovc(1);
        fullfecmodel += ":";
        fullfecmodel += size(1);
        
        fecsuffix += " + ";
        fecsuffix += indcovc(1);
        fecsuffix += ":";
        fecsuffix += size(1);
        
        if (fectime == 3) {
          fullfecmodel += " + ";
          fullfecmodel += indcovc(1);
          fullfecmodel += ":";
          fullfecmodel += repst(1);
          
          fecsuffix += " + ";
          fecsuffix += indcovc(1);
          fecsuffix += ":";
          fecsuffix += repst(1);
        }
        
        if (historical && indcovc(2) != "none") {
          fullfecmodel += " + ";
          fullfecmodel += indcovc(2);
          fullfecmodel += ":";
          fullfecmodel += size(2);
          
          fullfecmodel += " + ";
          fullfecmodel += indcovc(2);
          fullfecmodel += ":";
          fullfecmodel += repst(2);
          
          fecsuffix += " + ";
          fecsuffix += indcovc(2);
          fecsuffix += ":";
          fecsuffix += size(2);
          
          fecsuffix += " + ";
          fecsuffix += indcovc(2);
          fecsuffix += ":";
          fecsuffix += repst(2);
        }
        if (indcova(1) != "none") {
          fullfecmodel += " + ";
          fullfecmodel += indcova(1);
          fullfecmodel += ":";
          fullfecmodel += indcovc(1);
          
          fecsuffix += " + ";
          fecsuffix += indcova(1);
          fecsuffix += ":";
          fecsuffix += indcovc(1);
          
          if (historical) {
            if (indcova(2) != "none" && indcovc(2) != "none") {
              fullfecmodel += " + ";
              fullfecmodel += indcova(2);
              fullfecmodel += ":";
              fullfecmodel += indcovc(2);
              
              fullfecmodel += " + ";
              fullfecmodel += indcova(1);
              fullfecmodel += ":";
              fullfecmodel += indcovc(2);
              
              fullfecmodel += " + ";
              fullfecmodel += indcova(2);
              fullfecmodel += ":";
              fullfecmodel += indcovc(1);
              
              fecsuffix += " + ";
              fecsuffix += indcova(2);
              fecsuffix += ":";
              fecsuffix += indcovc(2);
              
              fecsuffix += " + ";
              fecsuffix += indcova(1);
              fecsuffix += ":";
              fecsuffix += indcovc(2);
              
              fecsuffix += " + ";
              fecsuffix += indcova(2);
              fecsuffix += ":";
              fecsuffix += indcovc(1);
            } else if (indcova(2) != "none") {
              fullfecmodel += " + ";
              fullfecmodel += indcova(2);
              fullfecmodel += ":";
              fullfecmodel += indcovc(1);
              
              fecsuffix += " + ";
              fecsuffix += indcova(2);
              fecsuffix += ":";
              fecsuffix += indcovc(1);
            } else if (indcovc(2) != "none") {
              fullfecmodel += " + ";
              fullfecmodel += indcova(1);
              fullfecmodel += ":";
              fullfecmodel += indcovc(2);
              
              fecsuffix += " + ";
              fecsuffix += indcova(1);
              fecsuffix += ":";
              fecsuffix += indcovc(2);
            }
          } 
        }
        if (indcovb(1) != "none") {
          fullfecmodel += " + ";
          fullfecmodel += indcovb(1);
          fullfecmodel += ":";
          fullfecmodel += indcovc(1);
          
          fecsuffix += " + ";
          fecsuffix += indcovb(1);
          fecsuffix += ":";
          fecsuffix += indcovc(1);
          if (historical) {
            if (indcovb(2) != "none" && indcovc(2) != "none") {
              fullfecmodel += " + ";
              fullfecmodel += indcovb(2);
              fullfecmodel += ":";
              fullfecmodel += indcovc(2);
              
              fullfecmodel += " + ";
              fullfecmodel += indcovb(1);
              fullfecmodel += ":";
              fullfecmodel += indcovc(2);
              
              fullfecmodel += " + ";
              fullfecmodel += indcovb(2);
              fullfecmodel += ":";
              fullfecmodel += indcovc(1);
              
              fecsuffix += " + ";
              fecsuffix += indcovb(2);
              fecsuffix += ":";
              fecsuffix += indcovc(2);
              
              fecsuffix += " + ";
              fecsuffix += indcovb(1);
              fecsuffix += ":";
              fecsuffix += indcovc(2);
              
              fecsuffix += " + ";
              fecsuffix += indcovb(2);
              fecsuffix += ":";
              fecsuffix += indcovc(1);
            } else if (indcovb(2) != "none") {
              fullfecmodel += " + ";
              fullfecmodel += indcovb(2);
              fullfecmodel += ":";
              fullfecmodel += indcovc(1);
              
              fecsuffix += " + ";
              fecsuffix += indcovb(2);
              fecsuffix += ":";
              fecsuffix += indcovc(1);
            } else if (indcovc(2) != "none") {
              fullfecmodel += " + ";
              fullfecmodel += indcovb(1);
              fullfecmodel += ":";
              fullfecmodel += indcovc(2);
              
              fecsuffix += " + ";
              fecsuffix += indcovb(1);
              fecsuffix += ":";
              fecsuffix += indcovc(2);
            }
          } 
        }
        
      }
    } else if (suite == "main") {
      
      if (age != "none" || covcount > 0) {
        fullfecmodel += " + ";
        fecsuffix += " + ";
      }
      
      fullfecmodel += size(1);
      fecsuffix += size(1);
      
      if (historical && sizecheck) {
        
        fullfecmodel += " + ";
        fullfecmodel += size(2);
        
        fecsuffix += " + ";
        fecsuffix += size(2);
      }
      
      if (fectime == 3) {
        fullfecmodel += " + ";
        fullfecmodel += repst(1);
        
        fecsuffix += " + ";
        fecsuffix += repst(1);
      }
      
      if (historical && repstcheck) {
        fullfecmodel += " + ";
        fullfecmodel += repst(2);
        
        fecsuffix += " + ";
        fecsuffix += repst(2);
      }
      
    } else if (suite == "size") {
      
      if (age != "none" || covcount > 0) {
        fullfecmodel += " + ";
        fecsuffix += " + ";
      }
      
      fullfecmodel += size(1);
      fecsuffix += size(1);
      
      if (historical && sizecheck) {
        fullfecmodel += " + ";
        fullfecmodel += size(2);
        
        fullfecmodel += " + ";
        fullfecmodel += size(1);
        fullfecmodel += ":";
        fullfecmodel += size(2);
        
        fecsuffix += " + ";
        fecsuffix += size(2);
        
        fecsuffix += " + ";
        fecsuffix += size(1);
        fecsuffix += ":";
        fecsuffix += size(2);
      }
      
    } else if (suite == "rep") {
      
      if (age != "none" || covcount > 0) {
        fullfecmodel += " + ";
        fecsuffix += " + ";
      }
      
      if (fectime == 3) {
        fullfecmodel += repst(1);
        fecsuffix += repst(1);
      } else if (!historical) {
        fullfecmodel += "1";
      }
      
      if (historical && repstcheck) {
        fullfecmodel += " + ";
        fullfecmodel += repst(2);
        
        fullfecmodel += " + ";
        fullfecmodel += repst(1);
        fullfecmodel += ":";
        fullfecmodel += repst(2);
        
        fecsuffix += " + ";
        fecsuffix += repst(2);
        
        fecsuffix += " + ";
        fecsuffix += repst(1);
        fecsuffix += ":";
        fecsuffix += repst(2);
      }
      
    } else {
      if (age == "none" && covcount == 0) fullfecmodel += "1";
    }
    
    fullfecmodel += fixedtackon;
    fecsuffix += fixedtackon;
    fullfecmodel += randomtackon;
    
    if (fec0 && suite != "const" && approach == "glm") {
      if (fecdist == "poisson") {
        fullfecmodel += " | ";
        fullfecmodel += fecsuffix;
      }
    }
    
  } else {
    fullfecmodel = "none";
  }
  
  StringVector fullnames {"time t", "individual", "patch", "alive in time t+1", "observed in time t+1", 
                          "size in time t+1", "reproductive status in time t+1", "fecundity in time t+1",
                          "fecundity in time t", "size in time t", "size in time t-1", 
                          "reproductive status in time t", "reprodutive status in time t-1", "age in time t",
                          "individual covariate a in time t", "individual covariate a in time t-1", 
                          "individual covariate b in time t", "individual covariate b in time t-1",
                          "individual covariate c in time t", "individual covariate c in time t-1",};
  StringVector mainparams {"year2", "individ", "patch", "surv3", "obs3", "size3", "repst3", "fec3", "fec2", 
                           "size2", "size1", "repst2", "repst1", "age", "indcova2", "indcova1", "indcovb2",
                           "indcovb1", "indcovc2", "indcovc1"};
  
  StringVector modelparams (20);
  modelparams(0) = year;
  modelparams(1) = indiv;
  modelparams(2) = patch;
  modelparams(3) = surv(0);
  modelparams(4) = obs(0);
  modelparams(5) = size(0);
  modelparams(6) = repst(0);
  modelparams(7) = fec(0);
  modelparams(8) = fec(1);
  modelparams(9) = size(1);
  if (historical) {modelparams(10) = size(2);} else {modelparams(10) = "none";}
  modelparams(11) = repst(1);
  if (historical) {modelparams(12) = repst(2);} else {modelparams(12) = "none";}
  modelparams(13) = age;
  modelparams(14) = indcova(1);
  if (historical) {modelparams(15) = indcova(2);} else {modelparams(15) = "none";}
  modelparams(16) = indcovb(1);
  if (historical) {modelparams(17) = indcovb(2);} else {modelparams(17) = "none";}
  modelparams(18) = indcovc(1);
  if (historical) {modelparams(19) = indcovc(2);} else {modelparams(19) = "none";}
  
  Rcpp::DataFrame paramnames = DataFrame::create(Named("parameter_names") = fullnames, _["mainparams"] = mainparams,
                                                 _["modelparams"] = modelparams);
  
  Rcpp::List output = List::create(Named("full.surv.model") = fullsurvmodel, _["full.obs.model"] = fullobsmodel,  
                                   _["full.size.model"] = fullsizemodel, _["full.repst.model"] = fullrepstmodel,  
                                   _["full.fec.model"] = fullfecmodel,  _["juv.surv.model"] = juvsurvmodel,
                                   _["juv.obs.model"] = juvobsmodel, _["juv.size.model"] = juvsizemodel, 
                                   _["juv.repst.model"] = juvrepstmodel, _["paramnames"] = paramnames);
  
  if (fullsurvmodel == "none") {
    output["full.surv.model"] = 1;
    
    if (!nojuvs) {output["juv.surv.model"] = 1;}
  }
  if (nojuvs) {output["juv.surv.model"] = 0;}
  
  if (fullobsmodel == "none") {
    output["full.obs.model"] = 1;
    if (!nojuvs) {output["juv.obs.model"] = 1;}
  }
  if (nojuvs) {output["juv.obs.model"] = 0;}
  
  if (fullsizemodel == "none") {
    if (!nojuvs) {output["juv.size.model"] = 1;}
    output["full.size.model"] = 1;
  }
  if (nojuvs) {output["juv.size.model"] = 0;}
  
  if (fullrepstmodel == "none") {
    output["full.repst.model"] = 1;
    if (!nojuvs) {output["juv.repst.model"] = 1;}
  }
  if (nojuvs) {output["juv.repst.model"] = 0;}
  
  if (fullfecmodel == "none") output["full.fec.model"] = 1;
  
  return output;
}
