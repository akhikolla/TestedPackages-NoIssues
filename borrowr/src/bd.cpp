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

#include "bd.h"

bool bd(tree& x, xinfo& xi, dinfo& di, pinfo& pi, double sigma,
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
      double syl, syr; //sum of y in proposed bots
      getsuff(x,nx,v,c,xi,di,nl,syl,nr,syr);

      //--------------------------------------------------
      //compute alpha
      double alpha=0.0, lalpha=0.0;
      double lhl, lhr, lht;
      if((nl>=5) && (nr>=5)) { //cludge?
         lhl = lh(nl,syl,sigma,pi.gamma);
         lhr = lh(nr,syr,sigma,pi.gamma);
         lht = lh(nl+nr,syl+syr,sigma,pi.gamma);

         alpha=1.0;
         lalpha = log(pr) + (lhl+lhr-lht) + log(sigma);
         lalpha = std::min(0.0,lalpha);
      }

      //--------------------------------------------------
      //try metrop
      double mul,mur; //means for new bottom nodes, left and right
      double uu = gen.uniform();
      bool dostep = (alpha > 0) && (log(uu) < lalpha);
      if(dostep) {
         mul = drawnodemu(nl,syl,pi.gamma,sigma,gen);
         mur = drawnodemu(nr,syr,pi.gamma,sigma,gen);
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
      size_t nr,nl; //counts at bots of nx
      double syl, syr; //sum at bots of nx
      getsuff(x, nx->getl(), nx->getr(), xi, di, nl, syl, nr, syr);

      //--------------------------------------------------
      //compute alpha
      double lhl, lhr, lht;
      lhl = lh(nl,syl,sigma,pi.gamma);
      lhr = lh(nr,syr,sigma,pi.gamma);
      lht = lh(nl+nr,syl+syr,sigma,pi.gamma);

      double lalpha = log(pr) + (lht - lhl - lhr) - log(sigma);
      lalpha = std::min(0.0,lalpha);

      //--------------------------------------------------
      //try metrop
      double mu;
      if(log(gen.uniform()) < lalpha) {
         mu = drawnodemu(nl+nr,syl+syr,pi.gamma,sigma,gen);
	 nv[nx->getv()]--;
         x.deathp(nx,mu);
         return true;
      } else {
         return false;
      }
   }
}
