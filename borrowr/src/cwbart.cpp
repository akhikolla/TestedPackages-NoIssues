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

#include "tree.h"
#include "treefuns.h"
#include "info.h"
#include "bartfuns.h"
#include "bd.h"
#include "bart.h"
#include "heterbart.h"

#ifndef NoRcpp

#define TRDRAW(a, b) trdraw(a, b)
#define TEDRAW(a, b) tedraw(a, b)

RcppExport SEXP cwbart(
   SEXP _in,            //number of observations in training data
   SEXP _ip,		//dimension of x
   SEXP _inp,		//number of observations in test data
   SEXP _ix,		//x, train,  pxn (transposed so rows are contiguous in memory)
   SEXP _iy,		//y, train,  nx1
   SEXP _ixp,		//x, test, pxnp (transposed so rows are contiguous in memory)
   SEXP _im,		//number of trees
   SEXP _inc,		//number of cut points
   SEXP _ind,		//number of kept draws (except for thinnning ..)
   SEXP _iburn,		//number of burn-in draws skipped
   SEXP _ipower,
   SEXP _ibase,
   // SEXP _itau,
   SEXP _inu,
   SEXP _ilambda,
   //mine
   SEXP _igamma,
   //end mine
   SEXP _isigest,
   SEXP _iw,
   SEXP _idart,
   SEXP _itheta,
   SEXP _iomega,
   SEXP _igrp,
   SEXP _ia,
   SEXP _ib,
   SEXP _irho,
   SEXP _iaug,
   SEXP _inkeeptrain,
   SEXP _inkeeptest,
   SEXP _inkeeptestme,
   SEXP _inkeeptreedraws,
   SEXP _inprintevery,
//   SEXP _treesaslists,
   SEXP _Xinfo
)
{

   //--------------------------------------------------
   //process args
   size_t n = Rcpp::as<int>(_in);
   size_t p = Rcpp::as<int>(_ip);
   size_t np = Rcpp::as<int>(_inp);
   Rcpp::NumericVector  xv(_ix);
   double *ix = &xv[0];
   Rcpp::NumericVector  yv(_iy);
   double *iy = &yv[0];
   Rcpp::NumericVector  xpv(_ixp);
   double *ixp = &xpv[0];
   size_t m = Rcpp::as<int>(_im);
   Rcpp::IntegerVector _nc(_inc);
   int *numcut = &_nc[0];
   //size_t nc = Rcpp::as<int>(_inc);
   size_t nd = Rcpp::as<int>(_ind);
   size_t burn = Rcpp::as<int>(_iburn);
   double mybeta = Rcpp::as<double>(_ipower);
   double alpha = Rcpp::as<double>(_ibase);
   //double tau = Rcpp::as<double>(_itau);
   double nu = Rcpp::as<double>(_inu);
   double lambda = Rcpp::as<double>(_ilambda);
   //mine
   double gamma = Rcpp::as<double>(_igamma);
   //end mine
   double sigma=Rcpp::as<double>(_isigest);
   Rcpp::NumericVector  wv(_iw);
   double *iw = &wv[0];
   bool dart;
   if(Rcpp::as<int>(_idart)==1) dart=true;
   else dart=false;
   double a = Rcpp::as<double>(_ia);
   double b = Rcpp::as<double>(_ib);
   double rho = Rcpp::as<double>(_irho);
   bool aug;
   if(Rcpp::as<int>(_iaug)==1) aug=true;
   else aug=false;
   double theta = Rcpp::as<double>(_itheta);
   double omega = Rcpp::as<double>(_iomega);
   Rcpp::IntegerVector _grp(_igrp);
   int *grp = &_grp[0];
   size_t nkeeptrain = Rcpp::as<int>(_inkeeptrain);
   size_t nkeeptest = Rcpp::as<int>(_inkeeptest);
   size_t nkeeptestme = Rcpp::as<int>(_inkeeptestme);
   size_t nkeeptreedraws = Rcpp::as<int>(_inkeeptreedraws);
   size_t printevery = Rcpp::as<int>(_inprintevery);
//   int treesaslists = Rcpp::as<int>(_treesaslists);
   Rcpp::NumericMatrix Xinfo(_Xinfo);
//   Rcpp::IntegerMatrix varcount(nkeeptreedraws, p);

   //return data structures (using Rcpp)
   Rcpp::NumericVector trmean(n); //train
   Rcpp::NumericVector temean(np);
   Rcpp::NumericVector sdraw(nd+burn);
   Rcpp::NumericMatrix trdraw(nkeeptrain,n);
   Rcpp::NumericMatrix tedraw(nkeeptest,np);
//   Rcpp::List list_of_lists(nkeeptreedraws*treesaslists);
   Rcpp::NumericMatrix varprb(nkeeptreedraws,p);
   Rcpp::IntegerMatrix varcnt(nkeeptreedraws,p);

   //random number generation
   arn gen;

   heterbart bm(m);

   if(Xinfo.size()>0) {
     xinfo _xi;
     _xi.resize(p);
     for(size_t i=0;i<p;i++) {
       _xi[i].resize(numcut[i]);
       //Rcpp::IntegerVector cutpts(Xinfo[i]);
       for(size_t j=0;j<numcut[i];j++) _xi[i][j]=Xinfo(i, j);
     }
     bm.setxinfo(_xi);
   }
#else

#define TRDRAW(a, b) trdraw[a][b]
#define TEDRAW(a, b) tedraw[a][b]

void cwbart(
   size_t n,            //number of observations in training data
   size_t p,		//dimension of x
   size_t np,		//number of observations in test data
   double* ix,		//x, train,  pxn (transposed so rows are contiguous in memory)
   double* iy,		//y, train,  nx1
   double* ixp,		//x, test, pxnp (transposed so rows are contiguous in memory)
   size_t m,		//number of trees
   int* numcut,		//number of cut points
   size_t nd,		//number of kept draws (except for thinnning ..)
   size_t burn,		//number of burn-in draws skipped
   double mybeta,
   double alpha,
   //double tau,
   double nu,
   double lambda,
   //mine
   double gamma,
   //end mine
   double sigma,
   double* iw,
   bool dart,
   double theta,
   double omega,
   int *grp,
   double a,
   double b,
   double rho,
   bool aug,
   size_t nkeeptrain,
   size_t nkeeptest,
   size_t nkeeptestme,
   size_t nkeeptreedraws,
   size_t printevery,
//   int treesaslists,
   unsigned int n1, // additional parameters needed to call from C++
   unsigned int n2,
   double* trmean,
   double* temean,
   double* sdraw,
   double* _trdraw,
   double* _tedraw
)
{

   //return data structures (using C++)
   std::vector<double*> trdraw(nkeeptrain);
   std::vector<double*> tedraw(nkeeptest);

   for(size_t i=0; i<nkeeptrain; ++i) trdraw[i]=&_trdraw[i*n];
   for(size_t i=0; i<nkeeptest; ++i) tedraw[i]=&_tedraw[i*np];

   std::vector< std::vector<size_t> > varcnt;
   std::vector< std::vector<double> > varprb;

   //random number generation
   arn gen(n1, n2);

   heterbart bm(m);
#endif

   for(size_t i=0;i<n;i++) trmean[i]=0.0;
   for(size_t i=0;i<np;i++) temean[i]=0.0;

   //printf("*****Into main of wbart\n");
   //-----------------------------------------------------------

   size_t skiptr,skipte,skipteme,skiptreedraws;
   if(nkeeptrain) {skiptr=nd/nkeeptrain;}
   else skiptr = nd+1;
   if(nkeeptest) {skipte=nd/nkeeptest;}
   else skipte=nd+1;
   if(nkeeptestme) {skipteme=nd/nkeeptestme;}
   else skipteme=nd+1;
   if(nkeeptreedraws) {skiptreedraws = nd/nkeeptreedraws;}
   else skiptreedraws=nd+1;

   //--------------------------------------------------
   //print args
   //printf("*****Data:\n");
   //printf("data:n,p,np: %zu, %zu, %zu\n",n,p,np);
   //printf("y1,yn: %lf, %lf\n",iy[0],iy[n-1]);
   //printf("x1,x[n*p]: %lf, %lf\n",ix[0],ix[n*p-1]);
   //if(np) printf("xp1,xp[np*p]: %lf, %lf\n",ixp[0],ixp[np*p-1]);
   //printf("*****Number of Trees: %zu\n",m);
   //printf("*****Number of Cut Points: %d ... %d\n", numcut[0], numcut[p-1]);
   //printf("*****burn and ndpost: %zu, %zu\n",burn,nd);
   //printf("*****Prior:beta,alpha,nu,lambda, gamma: %lf,%lf,%lf,%lf,%.3f\n",
                   // mybeta,alpha,nu,lambda,gamma);
   //printf("*****sigma: %lf\n",sigma);
   //printf("*****w (weights): %lf ... %lf\n",iw[0],iw[n-1]);
   //cout << "*****Dirichlet:sparse,theta,omega,a,b,rho,augment: "
	//<< dart << ',' << theta << ',' << omega << ',' << a << ','
	//<< b << ',' << rho << ',' << aug << endl;
   //printf("*****nkeeptrain,nkeeptest,nkeeptestme,nkeeptreedraws: %zu,%zu,%zu,%zu\n",
               // nkeeptrain,nkeeptest,nkeeptestme,nkeeptreedraws);
   //printf("*****printevery: %zu\n",printevery);
   //printf("*****skiptr,skipte,skipteme,skiptreedraws: %zu,%zu,%zu,%zu\n",skiptr,skipte,skipteme,skiptreedraws);

   //--------------------------------------------------
   //heterbart bm(m);
   //bm.setprior(alpha,mybeta,tau,gamma);
   bm.setprior(alpha,mybeta,gamma);
   bm.setdata(p,n,ix,iy,numcut);
   bm.setdart(a,b,rho,aug,dart,theta,omega);

   //--------------------------------------------------
   //sigma
   //gen.set_df(n+nu);
   double *svec = new double[n];
   for(size_t i=0;i<n;i++) svec[i]=iw[i]*sigma;

   //--------------------------------------------------

   std::stringstream treess;  //string stream to write trees to
   treess.precision(10);
   treess << nkeeptreedraws << " " << m << " " << p << endl;
   // dart iterations
   std::vector<double> ivarprb (p,0.);
   std::vector<size_t> ivarcnt (p,0);

   //--------------------------------------------------
   //temporary storage
   //out of sample fit
   double* fhattest=0; //posterior mean for prediction
   if(np) { fhattest = new double[np]; }
   double restemp=0.0,rss=0.0;


   //--------------------------------------------------
   //mcmc
   //printf("\nMCMC\n");
   //size_t index;
   size_t trcnt=0; //count kept train draws
   size_t tecnt=0; //count kept test draws
   size_t temecnt=0; //count test draws into posterior mean
   size_t treedrawscnt=0; //count kept bart draws
   bool keeptest,keeptestme,keeptreedraw;

   time_t tp;
   int time1 = time(&tp);
   xinfo& xi = bm.getxinfo();

   // mine, trying to get bottom nodes.
   //tree::npv bnv;
   //tree tt;
   //tt << bm.gettree(0);
   //for(int jj = 0; jj < m; jj++) {
    // bm.gettree(jj).getbots(bnv);
     //int ss;
     //ss = bnv.size();
     //printf("number of bottom nodes in tree %u? : %u\n", jj, ss);
   //}

   for(size_t i=0;i<(nd+burn);i++) {
      // if(i%printevery==0) printf("done %zu (out of %lu)\n",i,nd+burn);
      if(i==(burn/2)&&dart) bm.startdart();
      //draw bart
      bm.draw(svec,gen);
      //draw sigma
      rss=0.0;
      //mine, to print sample size
      //printf("n: %u\n", n);
      //mine
      //double m = 0.0;
      for(size_t k=0;k<n;k++) {
        //mine, to print y values
        //printf("iy[%u]: %f\n", k + 1, iy[k]);
        //mine, to compute mean of y values
        //m += iy[k] / n;
        //mine, to print fitted values
        //printf("bm.f(%u): %f\n", k + 1, bm.f(k));
        //mine, to print weights
        //printf("iw[%u]: %f\n", k + 1, iw[k]);
        restemp=(iy[k]-bm.f(k))/(iw[k]); rss += restemp*restemp;
      }
      //mine, to print mean
      //printf("m: %f\n", m);

      //mine, code to compute squared value of terminal node values
      int J = 0; // total number of terminal nodes
      double ssmu = 0.0; //sum of squared terminal node values
      for(int jj = 0; jj < m; jj++) {
        tree::npv bnv;
        bm.gettree(jj).getbots(bnv);
        J += bnv.size();
        for(int nn=0;nn!=bnv.size();nn++) {
          double nodeval = bnv[nn]->gettheta();
          ssmu += nodeval*nodeval;
        }
      }
      //to print number of terminal nodes and sum of squared mu's
      //if(i == (nd+burn - 1)) printf("J: %u, sum-of-squares mu: %f\n", J, ssmu);
      //double gamma = 20.0 / m; // needs to be added as an argument!
      double astar = (n + J + nu) / 2;
      //double foo = (1 / (nu * lambda)) * rss; //debugging
      //double bar = 1 / (gamma * nu * lambda) * ssmu; //debugging
      //printf("foo: %f, bar: %f\n", foo, bar); //debugging
      double bstar = 0.5 * ((1.0 / (nu * lambda)) * rss + (1.0 / (gamma * nu * lambda)) * ssmu + 1.0);
      //if(i == (nd+burn - 1))
        //printf("gamma: %.3f, rss: %.3f, ssmu: %.3f, astar: %.1f, bstar: %.3f\n", gamma, rss, ssmu, astar, bstar);
      //end mine

      //theirs
      //sigma = sqrt((nu*lambda + rss)/gen.chi_square(n+nu));
      // end theirs
      //mine
      sigma = sqrt(nu * lambda * 1.0 / gen.gamma(astar, bstar));
      //end mine
      for(size_t k=0;k<n;k++) svec[k]=iw[k]*sigma;
      sdraw[i]=sigma;
      if(i>=burn) {
         for(size_t k=0;k<n;k++) trmean[k]+=bm.f(k);
         if(nkeeptrain && (((i-burn+1) % skiptr) ==0)) {
            //index = trcnt*n;;
            //for(size_t k=0;k<n;k++) trdraw[index+k]=bm.f(k);
            for(size_t k=0;k<n;k++) TRDRAW(trcnt,k)=bm.f(k);
            trcnt+=1;
         }
         keeptest = nkeeptest && (((i-burn+1) % skipte) ==0) && np;
         keeptestme = nkeeptestme && (((i-burn+1) % skipteme) ==0) && np;
         if(keeptest || keeptestme) bm.predict(p,np,ixp,fhattest);
         if(keeptest) {
            //index=tecnt*np;
            //for(size_t k=0;k<np;k++) tedraw[index+k]=fhattest[k];
            for(size_t k=0;k<np;k++) TEDRAW(tecnt,k)=fhattest[k];
            tecnt+=1;
         }
         if(keeptestme) {
            for(size_t k=0;k<np;k++) temean[k]+=fhattest[k];
            temecnt+=1;
         }
         keeptreedraw = nkeeptreedraws && (((i-burn+1) % skiptreedraws) ==0);
         if(keeptreedraw) {
//	   #ifndef NoRcpp
//	   Rcpp::List lists(m*treesaslists);
//	   #endif

            for(size_t j=0;j<m;j++) {
	      treess << bm.gettree(j);
/*
	      #ifndef NoRcpp
	      varcount.row(treedrawscnt)=varcount.row(treedrawscnt)+bm.gettree(j).tree2count(p);
	      if(treesaslists) lists(j)=bm.gettree(j).tree2list(xi, 0., 1.);
	      #endif
*/
	    }
            #ifndef NoRcpp
//	    if(treesaslists) list_of_lists(treedrawscnt)=lists;
	    ivarcnt=bm.getnv();
	    ivarprb=bm.getpv();
	    size_t k=(i-burn)/skiptreedraws;
	    for(size_t j=0;j<p;j++){
	      varcnt(k,j)=ivarcnt[j];
	      //varcnt(i-burn,j)=ivarcnt[j];
	      varprb(k,j)=ivarprb[j];
	      //varprb(i-burn,j)=ivarprb[j];
	    }
            #else
	    varcnt.push_back(bm.getnv());
	    varprb.push_back(bm.getpv());
	    #endif

            treedrawscnt +=1;
         }
      }

     //mine, trying to get bottom nodes
     //if(i == (nd+burn - 1)) {
      //for(int jj = 0; jj < m; jj++) {
        //tree::npv bnv;
        //bm.gettree(jj).getbots(bnv);
        //int ss;
        //ss = bnv.size();
        //for(int foo = 0; foo < bnv.size(); foo ++) {printf("%u, %f\n", foo, bnv[foo]);}
        //printf("number of bottom nodes in tree %u: %u\n", jj + 1, bnv.size());
        //for(tree::npv::size_type nn=0;nn!=bnv.size();nn++) {
          //printf("tree: %u, node: %u, value: %f\n", jj + 1, nn + 1, bnv[nn]->gettheta());
        //}
      //}
     //}
     //end mine

   }

   //mine, compute posterior sd of terminal nodes
   int nnodes = 0;
   double nodesumsq = 0.0;
   for(int jj = 0; jj < m; jj++) {
     tree::npv bnv;
     bm.gettree(jj).getbots(bnv);
     nnodes += bnv.size();
     for(int nn=0;nn!=bnv.size();nn++) {
       double mu = bnv[nn]->gettheta();
       nodesumsq += mu*mu;
     }
   }
   double node_sd = sqrt(nodesumsq / (nnodes - 1));
   //end mine

   int time2 = time(&tp);
   //printf("time: %ds\n",time2-time1);
   for(size_t k=0;k<n;k++) trmean[k]/=nd;
   for(size_t k=0;k<np;k++) temean[k]/=temecnt;
   //printf("check counts\n");
   //printf("trcnt,tecnt,temecnt,treedrawscnt: %zu,%zu,%zu,%zu\n",trcnt,tecnt,temecnt,treedrawscnt);
   //--------------------------------------------------
   //PutRNGstate();

   if(fhattest) delete[] fhattest;
   if(svec) delete [] svec;

   //--------------------------------------------------
   //return
#ifndef NoRcpp
   Rcpp::List ret;
   ret["sigma"]=sdraw;
   ret["yhat.train.mean"]=trmean;
   ret["yhat.train"]=trdraw;
   ret["yhat.test.mean"]=temean;
   ret["yhat.test"]=tedraw;
   //ret["varcount"]=varcount;
   ret["varcount"]=varcnt;
   ret["varprob"]=varprb;

   //mine
   ret["node_sd"] = node_sd;
   //end mine

   //for(size_t i=0;i<m;i++) {
    //  bm.gettree(i).pr();
   //}

   Rcpp::List xiret(xi.size());
   for(size_t i=0;i<xi.size();i++) {
      Rcpp::NumericVector vtemp(xi[i].size());
      std::copy(xi[i].begin(),xi[i].end(),vtemp.begin());
      xiret[i] = Rcpp::NumericVector(vtemp);
   }

   Rcpp::List treesL;
   //treesL["nkeeptreedraws"] = Rcpp::wrap<int>(nkeeptreedraws); //in trees
   //treesL["ntree"] = Rcpp::wrap<int>(m); //in trees
   //treesL["numx"] = Rcpp::wrap<int>(p); //in cutpoints
   treesL["cutpoints"] = xiret;
   treesL["trees"]=Rcpp::CharacterVector(treess.str());
//   if(treesaslists) treesL["lists"]=list_of_lists;
   ret["treedraws"] = treesL;

   return ret;
#else

#endif

}
