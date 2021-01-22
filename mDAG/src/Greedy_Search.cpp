#include <RcppArmadillo.h>
#include <vector>
#include <iostream>
#include <random>
// [[Rcpp::depends("RcppArmadillo")]]
using namespace Rcpp;
using std::log;
// [[Rcpp::export]]
List gaussian_cpp(){
  // Obtain environment containing function
  Rcpp::Environment stats("package:stats"); 
  // Make function callable from C++
  Rcpp::Function  gaussian("gaussian");
  List out=gaussian();
  return out;
}

// [[Rcpp::export]]
List binomial_cpp(){
  // Obtain environment containing function
  Rcpp::Environment stats("package:stats"); 
  // Make function callable from C++
  Rcpp::Function  binomial("binomial");
  List out=binomial();
  return out;
}


// [[Rcpp::export]]
List poisson_cpp(){
  // Obtain environment containing function
  Rcpp::Environment stats("package:stats"); 
  // Make function callable from C++
  Rcpp::Function  poisson("poisson");
  List out=poisson();
  return out;
}


// [[Rcpp::export]]
List glm_fit_cpp( arma::mat x, arma::mat y, std::vector<double> weights ,List family  ){
  
  // Obtain environment containing function
  Rcpp::Environment stats("package:stats"); 
  // Make function callable from C++
  Rcpp::Function  glm_fit("glm.fit");
  List l=glm_fit(_["x"]=x,_["y"]=y,_["weights"]=weights,_["family"]=family);
  return l;
}


// [[Rcpp::export]]
double logLik_cpp(List fit,int samplesize){
  //logLik_cpp is calculating -0.5*BIC
  // AIC=-2*log_likelihood+2*(number of parameters)
  // BIC=-2*log_likelihood+(number of parameters)*log(sample size)
  // -0.5*BIC=-0.5*AIC+(number of parameters)-0.5*(number of parameters)*log(sample size)
  List fam = fit["family"];
  int p = fit["rank"];
  std::string family=fam["family"];
  if (family=="gaussian") p=p+1;
  double AIC=fit["aic"];
  double val = p - AIC/2;
  return val-0.5*log((double)samplesize)*p;
}

// [[Rcpp::export]]
List logistf_control_cpp(){
  // Obtain environment containing function
  Rcpp::Environment stats("package:logistf");
  // Make function callable from C++
  Rcpp::Function  logistf_control("logistf.control");
  List out=logistf_control();
  return out;
}

// [[Rcpp::export]]
List logistf_fit_cpp(arma::mat x, arma::mat y, std::vector<double> weights ,List control){
  Rcpp::Environment utils("package:utils");
  Rcpp::Environment logistf("package:logistf");
  Rcpp::Function getFromNamespace("getFromNamespace");
  Rcpp::Function f=getFromNamespace("logistf.fit", "logistf");
  List l=f(_["x"]=x,_["y"]=y,_["weight"]=weights,_["control"]=control);
  return l;
}

// [[Rcpp::export]]
double multinom_BIC_cpp(arma::mat x, arma::mat y, std::vector<double> weights) {
  // Make function callable from C++
  Rcpp::Function  multinom_BIC("multinom_BIC");
  List l=multinom_BIC(_["x"]=x,_["y"]=y,_["weights"]=weights);
  return(l["out"]);
}


// [[Rcpp::export]]
double  ScoreNodeWithNoneParents(std::vector< std::string > type, std::vector<int> level, int v, arma::mat data,std::vector<double> weights){
  List fam,fit;
  double LL_null;
  if (type.at(v)=="c"){
    if (level.at(v)==2){
      fam=binomial_cpp();
    }
  }
  if (type.at(v)=="g") fam=gaussian_cpp();
  if (type.at(v)=="p") fam=poisson_cpp();
  arma::mat one;
  one.ones(data.n_rows,1);
  if (type.at(v)=="c"){
    if (level.at(v)==2){
      List control=logistf_control_cpp();
      fit=logistf_fit_cpp(one,data.col(v),weights,control);
      LL_null=fit["loglik"];
      LL_null=LL_null-0.5*log((double)data.n_rows)*1;
    }else{
      LL_null=multinom_BIC_cpp(one,data.col(v),weights);
    }
    
  }else{
    fit=glm_fit_cpp(one,data.col(v),weights,fam);
    LL_null=logLik_cpp(fit,data.n_rows);
  }
  
  return LL_null;
}

// 'InitScore' gives us the initial score for each node based on empty graph
// [[Rcpp::export]]
NumericVector InitScore(std::vector< std::string > type, std::vector<int> level, arma::mat data,std::vector<double> weights){
  NumericVector scores(data.n_cols);
  for (int v=0;v<data.n_cols;v++){
    scores[v]=ScoreNodeWithNoneParents(type,level,v,data,weights);
    
  }
  return scores;
}

// 'ReturnParents' returns the parents of the i_th node given graph 'AdjMat'
// [[Rcpp::export]]
NumericVector ReturnParents(int i,IntegerMatrix AdjMat){
  NumericVector parents;
  for (int j=0;j<AdjMat.ncol();j++){
    if (AdjMat(j,i)==1)
      parents.push_back(j);
  }
  return parents;
}


// [[Rcpp::export]]
arma::mat subcolMatrix(arma::mat matrix, NumericVector index){
  arma::mat out;
  int i,j;
  out.zeros(matrix.n_rows,index.length());
  for (i=0;i<index.length();i++){
    j=index(i);
    out.col(i)=matrix.col(j);
  }
  return out;
}

// 'ScoreGraph' calculates the score for each node given graph 'AdjMat' 
// [[Rcpp::export]]
NumericVector ScoreGraph(std::vector< std::string > type, std::vector<int> level, arma::mat data,std::vector<double> weights,IntegerMatrix AdjMat){
  NumericVector tmp(data.n_cols);
  List fam,fit;
  for (int v=0;v<data.n_cols;v++){
    NumericVector allparents=ReturnParents(v,AdjMat);
    arma::mat childvector=data.col(v);
    if (allparents.size()==0){
      tmp[v]=ScoreNodeWithNoneParents(type,level,v,data,weights);
    }else{
      if (type.at(v)=="c") {
        arma::mat one;
        one.ones(data.n_rows,1);
        arma::mat x=subcolMatrix(data,allparents);
        x.insert_cols(0,one);
        
        if (level.at(v)==2){
          List control=logistf_control_cpp();
          fit=logistf_fit_cpp(x,data.col(v),weights,control);
          int p=x.n_cols;
          double ll=fit["loglik"];
          tmp[v]=ll-0.5*log((double)data.n_rows)*p;
        }else{
          tmp[v]=multinom_BIC_cpp(x,data.col(v),weights);
        }

      }else{
        if (type.at(v)=="g") fam=gaussian_cpp();
        if (type.at(v)=="p") fam=poisson_cpp();
        fit=glm_fit_cpp(subcolMatrix(data,allparents),data.col(v),weights,fam);
        tmp[v]=logLik_cpp(fit,data.n_rows);
      }
    }
  }
  return tmp;
}

// 'SettingEdges' sets edges for initializing the graph at the first iteration step
// [[Rcpp::export]]
List SettingEdges(NumericVector scores,arma::mat data,List rst,std::vector< std::string > type, std::vector<int> level, std::vector<int> SNP,std::vector<double> weights){
  IntegerMatrix graph(data.n_cols);
  List nodes=rst["nodes"];
  
  for (int i=0;i<data.n_cols;i++){
    List currentnode=nodes.at(i);
    IntegerVector nbr_index=currentnode["nbr_index"];
    for (int j=0;j<nbr_index.size();j++){
      if (graph(i,nbr_index[j]-1)==0 && graph(nbr_index[j]-1,i)==0 && SNP.at(nbr_index[j]-1)==0 ){
        graph(i,nbr_index[j]-1)=1;
        double before=sum(scores);
        NumericVector tmp=ScoreGraph(type,level,data,weights,graph);
        double after=sum(tmp);
        if (after>before){
          scores=tmp;
        }else{
          graph(i,nbr_index[j]-1)=0;
        }
      }
      
    }
  }
 
  return List::create(Named("graph",graph),Named("scores",scores));
}


// [[Rcpp::export]]
void AddReverseDelete(IntegerMatrix AdjMat,NumericVector scores,arma::mat data ,List rst, std::vector< std::string > type, std::vector<int> level, std::vector<int> SNP,std::vector<double> weights){
  List nodes=rst["nodes"];
  double before, after;
  NumericVector tmp;
  for (int i=0;i<data.n_cols;i++){
    List currentnode=nodes.at(i);
    IntegerVector nbr_index=currentnode["nbr_index"];
    for (int j=0;j<nbr_index.size();j++){
        
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_int_distribution<int> distribution(1,10);
      int rnd= distribution(gen);
    
      //int rnd=(rand()%10)+1;
      
      if (AdjMat(i,nbr_index[j]-1)==0 && AdjMat(nbr_index[j]-1,i)==0 && SNP.at(nbr_index[j]-1)==0 ){
        // we only allow SNP --> other nodes
        AdjMat(i,nbr_index[j]-1)=1; //add an edge
        before=sum(scores);
        tmp=ScoreGraph(type,level,data,weights,AdjMat);
        after=sum(tmp);
        if (after>before){
          scores=tmp;
        }else{
          AdjMat(i,nbr_index[j]-1)=0;
        }
      }else if( AdjMat(i,nbr_index[j]-1)==1 && rnd>5 ){
        //delete an edge
        AdjMat(i,nbr_index[j]-1)=0; 
        before=sum(scores);
        tmp=ScoreGraph(type,level,data,weights,AdjMat);
        after=sum(tmp);
        if (after>before){
          scores=tmp;
        }else{
          AdjMat(i,nbr_index[j]-1)=1;
        }
      }else if(AdjMat(nbr_index[j]-1,i)==1 && rnd<=5 && SNP.at(nbr_index[j]-1)==0 ){
        //reverse an edge
        AdjMat(nbr_index[j]-1,i)=0; 
        AdjMat(i,nbr_index[j]-1)=1;
        before=sum(scores);
        tmp=ScoreGraph(type,level,data,weights,AdjMat);
        after=sum(tmp);
        if (after>before){
          scores=tmp;
        }else{
          AdjMat(nbr_index[j]-1,i)=1;
          AdjMat(i,nbr_index[j]-1)=0;
        }
      }
    }
  }
  
  
  
}

// [[Rcpp::export]]
IntegerMatrix GreedySearch(arma::mat data,std::vector< std::string > type, std::vector<int> level, std::vector<int> SNP , List rst,std::vector<double> weights){
  
  NumericVector scores=InitScore(type,level,data,weights);
  List initialgraph=SettingEdges(scores,data,rst,type,level,SNP,weights);
  scores=initialgraph["scores"];
  IntegerMatrix  graph=initialgraph["graph"];
  int count=0; 
  int n_iteration=100;
  for (int k=0;k<n_iteration;k++ ){
    if (count==5) {return graph;}
    double bef=sum(scores);
    AddReverseDelete(graph,scores,data,rst,type,level,SNP,weights);
    if (sum(scores)==bef) {
      count++;
    }
  }
  
  return graph;
}


/*** R
if(F){
  # output of multinom_BIC is -0.5*BIC
  multinom_BIC=function(x,y,weights){
    formula=y~x
    fit=multinom(formula,weights=weights,trace=F)
    p=fit$rank
    AIC=fit$AIC
    out=-0.5*AIC+p-0.5*p*log(length(weights)) # -0.5*BIC
    return( list("out"=out) )
  }
}




*/
