########## R-function: MarronWandParm ##########

# For obtaining the parameters for the Marron & Wand
# family of normal mixture distributions.

# Last changed: 25 OCT 2018

MarronWandParm <- function(densNum)
{
   if (!any(densNum==1:15)) stop("Illegal density number.")
  
   if (densNum==1)
   {
      w <- 1
      mu <- 0
      sigmasq <- 1
   }
   if (densNum==2)
   {
      w <- c(1/5,1/5,3/5)
      mu <- c(0,1/2,13/12)
      sigmasq <- c(1,4/9,25/81)
   }
   if (densNum==3)
   {
      w <- c(1/8,1/8,1/8,1/8,1/8,1/8,1/8,1/8)
      mu <- c(0,-1,-5/3,-19/9,-65/27,-211/81,-665/343,-2059/729)
      mu <- c(0/1,-3/3,-15/9,-57/27,-195/81,-633/243,-1995/729,-6177/2187)
      sigmasq <- c(1,4/9,(4/9)^2,(4/9)^3,(4/9)^4,(4/9)^5,(4/9)^6,(4/9)^7)
   }
   if (densNum==4)
   {
      w <- c(2/3,1/3)
      mu <- c(0,0)
      sigmasq <- c(1,1/100)
   }
   if (densNum==5)
   {               
      w <- c(1/10,9/10)
      mu <- c(0,0)
      sigmasq <- c(1,1/100)
       }
   if (densNum==6)
   {
      w <- c(1/2,1/2)
      mu <- c(-1,1)
      sigmasq <- c(4/9,4/9)
      }
   if (densNum==7)
   {
      w <- c(1/2,1/2)
      mu <- c(-3/2,3/2)
      sigmasq <- c(1/4,1/4)
    }
   if (densNum==8)
   {
      w <- c(3/4,1/4)
      mu <- c(0,3/2)
      sigmasq <- c(1,1/9)
   }
   if (densNum==9)
   {
      w <- c(9/20,9/20,1/10)
      mu <- c(-6/5,6/5,0)
      sigmasq <- c(9/25,9/25,1/16)
   }
   if (densNum==10)
   {
      w <- c(1/2,1/10,1/10,1/10,1/10,1/10)
      mu <- c(0,-1,-1/2,0,1/2,1)
      sigmasq <- c(1,1/100,1/100,1/100,1/100,1/100)  
   }
   if (densNum==11)
   {
      w <- c(0.49,0.49,rep((1/350),7))
      mu <- c(-1,1,((0:6)-3)/2)
      sigmasq <- c(4/9,4/9,rep(0.0001,7))
   }
   if (densNum==12)
   {
      w <- c(0.5,2^(c(3,2,1,0,-1))/31)
      mu <- c(0,((-2:2)+0.5))
      sigmasq <- c(1,1/(100*2^(2*(-2:2))))  
   }
   if (densNum==13)
   {
      w <- c(rep(138,2),rep(1,3),rep(7,3))/300
      mu <- c(-1,1,-0.5,-1,-1.5,0.5,1,1.5)
      sigmasq <- c(rep(4/9,2),rep(1/10000,3),rep(49/10000,3))
   }  
   if (densNum==14)
   {
      ellvec <- 0:5
      w <- (2^(5-ellvec))/63
      mu <- (65 - 96*(0.5^ellvec))/21
      sigmasq <-  (32/(63*2^ellvec))^2
   }  
   if (densNum==15)
   {
      w <- c(rep(6,3),rep(1,3))/21
      mu <- c(-15,-3,9,16,18,20)/7
      sigmasq <- c(rep(36,3),rep(1,3))/441
   }  
   return(list(w=w,mu=mu,sigmasq=sigmasq))
}

############ End of MarronWandParm ############

  
