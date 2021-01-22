.onLoad = function(libname, pkgname) {
  .jpackage(pkgname, lib.loc = libname)
  
  #need to check if proper Java is installed by special request of Prof Brian Ripley
  jv = .jcall("java/lang/System", "S", "getProperty", "java.runtime.version")
  if (substr(jv, 1L, 2L) == "1.") {
	  jvn = as.numeric(paste0(strsplit(jv, "[.]")[[1L]][1:2], collapse = "."))
	  if (jvn < 1.7){
		  warning(paste("Java v", jvn, "detected. Java 7 (at minimum) is needed for this package but is does not seem to be available. This message may be in error; apologies if it is."))
	  }		
  }
}

.onAttach = function(libname, pkgname){
  num_gigs_ram_available = .jcall(.jnew("java/lang/Runtime"), "J", "maxMemory") / 1e9
  packageStartupMessage(paste("Welcome to GreedyExperimentalDesign v", utils::packageVersion("GreedyExperimentalDesign"), ".\n", round(num_gigs_ram_available, 2), "GB memory available for the Java searcher.\n", sep = ""))
}