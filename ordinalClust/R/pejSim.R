pejSim <-
  function(ej,m,mu,p){
    # ----------------------------------------------------------------------------
    # Proba ej (unconditional) with a hierarchical random dichotomic model.
    # input:
    #   ej  interval ej+1 [binf,bsup]
    #   j   dichotomical step
    #   m   nb of modalities
    #   mu  true modality {1,...,m}
    #   p   P(sait) [0,1]
    #   z1tozjm1    values {0,1} of z_1 to z_{j-1} OPTIONAL
    #               1: z is known and is 1
    #               2: z is unknown and is 0
    # output:
    #   pej
    # ----------------------------------------------------------------------------
    # j == 1
    j=m
    if (j == 1){ 
      proba = 1
      return(proba)
    }
    # if ej is not an interval then build [ej,ej]
    if (length(ej) == 1) {ej = c(ej,ej)}
    # create zjm1toz1
    
    z1tozjm1 = matrix(0,1,j-1) # 0 at place j <=> unknown value for z_j
    
    z1tozjm2 = z1tozjm1[1:(length(z1tozjm1)-1)];
    zjm1 = z1tozjm1[length(z1tozjm1)];
    
    
    # version nouvelle plus rapide (exple de gain : 17* pour m=5, 374* pour m=6)
    # ----------------------------
    # idee de la mehode : on coupe les branches de proba nulle, cad ou ejm1 pas
    # inclus dans ej
    # j > 1
    if (zjm1) { # zjm1 is known   ####### ATTENTION : pas sur que ca marche sous R
      proba = 0
      tabint = allej(j-1,m)
      if (!is.matrix(tabint)) {tabint = t(as.matrix(allej(j-1,m)))}
      nbtabint = length(tabint[,1])
      for (i in 1:nbtabint){
        ejm1 = tabint[i,] # ej minus 1
        if ((ej[1] >= ejm1[1]) && (ej[2] <= ejm1[2])){ # pour accelerer, il faut verifier ejm1 est inclus dans ej !
          proba = proba + pejp1zj1_ej(ej,ejm1,mu,p) * pej(ejm1,j-1,m,mu,p,z1tozjm2);
        }
      }
    }
    else # zjm1 is unknown    
    {
      
      proba = 0
      tabint = allej(j-1,m)
      if (!is.matrix(tabint)) {tabint = t(as.matrix(allej(j-1,m)))}
      nbtabint = length(tabint[,1])
      for (i in 1:nbtabint){
        ejm1 = tabint[i,]# ej minus 1
        if ((ej[1] >= ejm1[1]) && (ej[2] <= ejm1[2])){ # pour accelerer, il faut verifier ejm1 est inclus dans ej !
          proba = proba + pejp1_ej(ej,ejm1,mu,p) * pej(ejm1,j-1,m,mu,p,z1tozjm2)
        }
      }
      
    } # zjm1
    return(proba)
  }
