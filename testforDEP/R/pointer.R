#########################################
# Print Environments
#########################################

printEnv = function(curE){
  print(curE)
  if(!identical(curE, emptyenv())){
    printEnv(parent.env(curE))
  }
}

########################################
# Add environment above GlobalEnv
########################################

pointer = function(obj, curE){

  name = deparse(substitute(obj))
  e = globalenv()
  par = parent.env(e)

#  if(exists(name, envir = parent.env(curE), inherits = TRUE)){
#    stop("Cannot assign pointer twice!!")
#  }
#  else{
    e1 = new.env()
    e1[[name]] = obj
    parent.env(e1) = par
    parent.env(e) = e1
    remove(list = name, envir = curE, inherits = FALSE)
    return (e1)
#  }
}

########################################
# Delete Environment
########################################

delete = function(eName, curE){

  del = function(e, name, curE){

    getChild = function(e, curE){

      if(identical(parent.env(curE), e)){return (curE)}
      if(identical(parent.env(curE), baseenv())){return (NULL)}
      return (getChild(e, parent.env(curE)))

    }

     child = getChild(e, curE)
     parent.env(child) = parent.env(e)
     parent.env(e) = emptyenv()
     remove(list = name, envir = curE, inherits = FALSE)

  }

  del(curE[[eName]], eName, curE)
}

########################################
# Example
########################################



