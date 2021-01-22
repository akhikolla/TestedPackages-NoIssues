/*
 *  borrowr: estimate population average treatment effects with borrowing between data sources.
 *  Copyright (C) 2019  Jeffrey A. Verdoliva Boatman
 *  This code is a modified version from the BART R package from April 2019.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.

 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.

 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "heterbd.h"

bool heterbd(tree& x, xinfo& xi, dinfo& di, pinfo& pi, double *sigma,
	     std::vector<size_t>& nv, std::vector<double>& pv, bool aug, rn& gen)
{
   tree::npv goodbots;  //nodes we could birth at (split on)
   double PBx = getpb(x,xi,pi,goodbots); //prob of a birth at x

   if(gen.uniform() < PBx) { //do birth or death

      //--------------------------------------------------
      //draw proposal
      tree::tree_p nx; //bottom node
      size_t v,c; //variable and cutpoint
      double pr; //part of metropolis ratio from proposal and prior
      bprop(x,xi,pi,goodbots,PBx,nx,v,c,pr,nv,pv,aug,gen);

      //--------------------------------------------------
      //compute sufficient statistics
      size_t nr,nl; //counts in proposed bots
      double bl,br; //sums of weights
      double Ml, Mr; //weighted sum of y in proposed bots
      hetergetsuff(x,nx,v,c,xi,di,nl,bl,Ml,nr,br,Mr,sigma);

      //--------------------------------------------------
      //compute alpha
      double alpha=0.0, lalpha=0.0;
      double lhl, lhr, lht;
      if((nl>=5) && (nr>=5)) { //cludge?
         lhl = heterlh(bl,Ml,sigma[0], pi.gamma);
         lhr = heterlh(br,Mr,sigma[0], pi.gamma);
         lht = heterlh(bl+br,Ml+Mr,sigma[0], pi.gamma);


         alpha=1.0;
         lalpha = log(pr) + (lhl+lhr-lht);
         lalpha = std::min(0.0,lalpha);
      }

      //--------------------------------------------------
      //try metrop
      double mul,mur; //means for new bottom nodes, left and right
      double uu = gen.uniform();
      bool dostep = (alpha > 0) && (log(uu) < lalpha);
      if(dostep) {
         mul = heterdrawnodemu(bl,Ml,sigma[0],pi.gamma,gen);
         mur = heterdrawnodemu(br,Mr,sigma[0],pi.gamma,gen);
         x.birthp(nx,v,c,mul,mur);
	 nv[v]++;
         return true;
      } else {
         return false;
      }
   } else {
      //--------------------------------------------------
      //draw proposal
      double pr;  //part of metropolis ratio from proposal and prior
      tree::tree_p nx; //nog node to death at
      dprop(x,xi,pi,goodbots,PBx,nx,pr,gen);

      //--------------------------------------------------
      //compute sufficient statistics
      double br,bl; //sums of weights
      double Ml, Mr; //weighted sums of y
      hetergetsuff(x, nx->getl(), nx->getr(), xi, di, bl, Ml, br, Mr, sigma);

      //--------------------------------------------------
      //compute alpha
      double lhl, lhr, lht;
      lhl = heterlh(bl,Ml,sigma[0], pi.gamma);
      lhr = heterlh(br,Mr,sigma[0], pi.gamma);
      lht = heterlh(bl+br,Ml+Mr,sigma[0], pi.gamma);


      double lalpha = log(pr) + (lht - lhl - lhr);
      lalpha = std::min(0.0,lalpha);

      //--------------------------------------------------
      //try metrop
      //double a,b,s2,yb;
      double mu;
      if(log(gen.uniform()) < lalpha) {
        mu = heterdrawnodemu(bl+br,Ml+Mr,sigma[0],pi.gamma,gen);
	 nv[nx->getv()]--;
         x.deathp(nx,mu);
         return true;
      } else {
         return false;
      }
   }
}
