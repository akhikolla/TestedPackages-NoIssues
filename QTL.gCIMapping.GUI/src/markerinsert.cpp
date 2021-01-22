#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
 NumericMatrix markerinsert(NumericMatrix mp,NumericMatrix geno,NumericMatrix map,int cl,int gg1,int gg2,int gg0,int flagRIL){  
  
  NumericMatrix markerpart(1,mp.nrow());
 
  for (int ire=0;ire<mp.nrow();ire++){
     int wp1=0;
     int wp2=0;

    if (mp(ire,1)==1)
    {
      wp1=1;
      wp2=0;
      for(int i=0;i<map.nrow();i++){
        if(map(i,0)==1){wp2+=1;}
      }
    }
    if (mp(ire,1)>1)
    {
      wp1=0;
      wp2=0;
      for(int i=0;i<map.nrow();i++){
        if(map(i,0)<=(mp(ire,1)-1)){wp1+=1;}
      }
      wp1=wp1+1;
      for(int i=0;i<map.nrow();i++){
        if(map(i,0)<=mp(ire,1)){wp2+=1;}
      }
    }
    int ichr=mp(ire,1);
    int q5count=0;
    for(int i=0;i<map.nrow();i++){
      if(map(i,0)==ichr){q5count+=1;}
    }
    IntegerMatrix q5(q5count,1);
    int j=0; 
    for(int i=0;i<map.nrow();i++){
      if(map(i,0)==ichr){
        q5(j,0)=i;
        j+=1; 
      }
    }
    NumericMatrix r0t(1,(q5count-1));
    for (int ii=0;ii<(q5count-1);ii++)
    {
      r0t(0,ii)=(1-exp(-2*(map(q5(ii+1,0),1)-map(q5(ii,0),1))/100))/2;
    }
    int mp2=(int)mp(ire,2);
    int mp21=(int)mp(ire,2)+1;
    double r01=(1-exp(-2*(mp(ire,0)-map(q5(mp2-1,0),1))/100))/2;
    double r02=(1-exp(-2*(-mp(ire,0)+map(q5(mp21-1,0),1))/100))/2;
    int mp3=(int)mp(ire,3);
    int mp1=(int)mp(ire,1);
    int abc=mp3+mp1-1;
    NumericMatrix ac4(1,2);
    NumericMatrix cn00(1,2);
    NumericMatrix dd2(2,2);
    dd2(0,0)=1.0;
    dd2(0,1)=0.0;
    dd2(1,0)=0.0;
    dd2(1,1)=0.0;
    int j_1=1;
    
    int uu1=wp1;
    int uu2=wp2;
    NumericMatrix az1(1,2);
    NumericMatrix r_tr(2,2);
    az1(0,0)=1.0;
    az1(0,1)=1.0;
    int nni=0;
    for (int itt =(uu1-1);itt<(uu2-1);itt++)
    {
      nni=nni+1;
      NumericMatrix dd1(2,2);
      dd1(0,0)=1.0;
      dd1(0,1)=0.0;
      dd1(1,0)=0.0;
      dd1(1,1)=1.0;
      if (geno(j_1-1,itt)==gg1){dd1(1,1)=0;}
      if (geno(j_1-1,itt)==gg2){dd1(0,0)=0;}
      if (geno(j_1-1,itt)==gg0){
        for (int id = 0;id < 2;id++){
          for (int jd = 0;jd < 2;jd++){
            dd1(id,jd) = dd1(id,jd)*0.5;
          }
        }
      }
      if ((abs(itt+1-abc))<1e-6)
      {
        NumericMatrix azdd(1,2);
        for (int i = 0; i < 1; i++)            
        {
          for (int j = 0; j < 2; j++)  
          {  
            azdd(i,j) = 0;  
            for (int k = 0; k < 2; k++){  
              azdd(i,j) += az1(i,k) * dd1(k,j);  
            }  
          }  
        }
        az1=azdd;
        if (flagRIL==0){
          r_tr(0,0)=1-r01;
          r_tr(0,1)=r01;
          r_tr(1,0)=r01;
          r_tr(1,1)=1-r01;
        }
        if (flagRIL==1)
        {
          r_tr(0,0)=1/(1+2*r01);
          r_tr(0,1)=2*r01/(1+2*r01);
          r_tr(1,0)=2*r01/(1+2*r01);
          r_tr(1,1)=1/(1+2*r01);
        }    
        NumericMatrix azrtd(1,2);
        NumericMatrix azrt(1,2);
        for (int i = 0; i < 1; i++)            
        {
          for (int j = 0; j < 2; j++)  
          {  
            azrt(i,j) = 0;  
            for (int k = 0; k < 2; k++){  
              azrt(i,j) += az1(i,k) * r_tr(k,j);  
            }  
          }  
        }
        for (int i = 0; i < 1; i++)            
        {
          for (int j = 0; j < 2; j++)  
          {  
            azrtd(i,j) = 0;  
            for (int k = 0; k < 2; k++){  
              azrtd(i,j) += azrt(i,k) * dd2(k,j);  
            }  
          }  
        }
        az1=azrtd;
        if (flagRIL==0){
          r_tr(0,0)=1-r02;
          r_tr(0,1)=r02;
          r_tr(1,0)=r02;
          r_tr(1,1)=1-r02;
        }
        if (flagRIL==1)
        {
          r_tr(0,0)=1/(1+2*r02);
          r_tr(0,1)=2*r02/(1+2*r02);
          r_tr(1,0)=2*r02/(1+2*r02);
          r_tr(1,1)=1/(1+2*r02);
        }    
        NumericMatrix az1rt(1,2);
        for (int i = 0; i < 1; i++)            
        {
          for (int j = 0; j < 2; j++)  
          {  
            az1rt(i,j) = 0;  
            for (int k = 0; k < 2; k++){  
              az1rt(i,j) += az1(i,k) * r_tr(k,j);  
            }  
          }  
        }
        az1=az1rt;
      }
      if ((abs(itt+1-abc))>1e-6)
      {
        NumericMatrix azdd1(1,2);
        for (int i = 0; i < 1; i++)            
        {
          for (int j = 0; j < 2; j++)  
          {  
            azdd1(i,j) = 0;  
            for (int k = 0; k < 2; k++){  
              azdd1(i,j) += az1(i,k) * dd1(k,j);  
            }  
          }  
        }
        az1=azdd1;
        if (flagRIL==0){
          r_tr(0,0)=1-r0t(0,(itt-uu1+1));
          r_tr(0,1)=r0t(0,(itt-uu1+1));
          r_tr(1,0)=r0t(0,(itt-uu1+1));
          r_tr(1,1)=1-r0t(0,(itt-uu1+1));
        }
        if (flagRIL==1)
        {
          r_tr(0,0)=1/(1+2*r0t(0,(itt-uu1+1)));
          r_tr(0,1)=2*r0t(0,(itt-uu1+1))/(1+2*r0t(0,(itt-uu1+1)));
          r_tr(1,0)=2*r0t(0,(itt-uu1+1))/(1+2*r0t(0,(itt-uu1+1)));
          r_tr(1,1)=1/(1+2*r0t(0,(itt-uu1+1)));
        } 
        NumericMatrix azrt2(1,2);
        for (int i = 0; i < 1; i++)            
        {
          for (int j = 0; j < 2; j++)  
          {  
            azrt2(i,j) = 0;  
            for (int k = 0; k < 2; k++){  
              azrt2(i,j) += az1(i,k) * r_tr(k,j);  
            }  
          }  
        }
        az1=azrt2; 
      }
    }
    NumericMatrix dd1(2,2);
    dd1(0,0)=1.0;
    dd1(0,1)=0.0;
    dd1(1,0)=0.0;
    dd1(1,1)=1.0;
    if (geno(j_1-1,uu2-1)==gg1){dd1(1,1)=0;}
    if (geno(j_1-1,uu2-1)==gg2){dd1(0,0)=0;}
    if (geno(j_1-1,uu2-1)==gg0){
      for (int id = 0;id < 2;id++){
        for (int jd = 0;jd < 2;jd++){ 
          dd1(id,jd)=dd1(id,jd)*0.5;
        }
      }
    }
    NumericMatrix az1dd1(1,2);
    for (int i = 0; i < 1; i++)            
    {
      for (int j = 0; j < 2; j++)  
      {  
        az1dd1(i,j) = 0;  
        for (int k = 0; k < 2; k++){  
          az1dd1(i,j) += az1(i,k) * dd1(k,j);  
        }  
      }  
    }
    az1=az1dd1;
    NumericMatrix vaz1(az1.ncol(),1);
    vaz1(0,0)=1.0;
    vaz1(1,0)=1.0;
    NumericMatrix az1va(1,1);
    for (int i = 0; i < 1; i++)            
    {
      for (int j = 0; j < 1; j++)  
      {  
        az1va(i,j) = 0;  
        for (int k = 0; k < 2; k++){  
          az1va(i,j) += az1(i,k) * vaz1(k,j);  
        }  
      }  
    }
    ac4(0,0)=az1va(0,0);
    az1(0,0)=1.0;
    az1(0,1)=1.0;
    nni=0;
    dd2(0,0)=0.0;
    dd2(1,1)=1.0;
    for (int itt =(uu1-1);itt<(uu2-1);itt++)
    {
      nni=nni+1;
      NumericMatrix dd1(2,2);
      dd1(0,0)=1.0;
      dd1(0,1)=0.0;
      dd1(1,0)=0.0;
      dd1(1,1)=1.0;
      if (geno(j_1-1,itt)==gg1){dd1(1,1)=0;}
      if (geno(j_1-1,itt)==gg2){dd1(0,0)=0;}
      if (geno(j_1-1,itt)==gg0){
        for (int id = 0;id < 2;id++){
          for (int jd = 0;jd < 2;jd++){
            dd1(id,jd) = dd1(id,jd)*0.5;
          }
        }
      }
      if ((abs(itt+1-abc))<1e-6)
      {
        NumericMatrix azdd(1,2);
        
        for (int i = 0; i < 1; i++)            
        {
          for (int j = 0; j < 2; j++)  
          {  
            azdd(i,j) = 0;  
            for (int k = 0; k < 2; k++){  
              azdd(i,j) += az1(i,k) * dd1(k,j);  
            }  
          }  
        }
        az1=azdd;
        if (flagRIL==0){
          r_tr(0,0)=1-r01;
          r_tr(0,1)=r01;
          r_tr(1,0)=r01;
          r_tr(1,1)=1-r01;
        }
        if (flagRIL==1)
        {
          r_tr(0,0)=1/(1+2*r01);
          r_tr(0,1)=2*r01/(1+2*r01);
          r_tr(1,0)=2*r01/(1+2*r01);
          r_tr(1,1)=1/(1+2*r01);
        }    
        NumericMatrix azrtd(1,2);
        NumericMatrix azrt(1,2);
        for (int i = 0; i < 1; i++)            
        {
          for (int j = 0; j < 2; j++)  
          {  
            azrt(i,j) = 0;  
            for (int k = 0; k < 2; k++){  
              azrt(i,j) += az1(i,k) * r_tr(k,j);  
            }  
          }  
        }
        for (int i = 0; i < 1; i++)            
        {
          for (int j = 0; j < 2; j++)  
          {  
            azrtd(i,j) = 0;  
            for (int k = 0; k < 2; k++){  
              azrtd(i,j) += azrt(i,k) * dd2(k,j);  
            }  
          }  
        }
        az1=azrtd;
        if (flagRIL==0){
          r_tr(0,0)=1-r02;
          r_tr(0,1)=r02;
          r_tr(1,0)=r02;
          r_tr(1,1)=1-r02;
        }
        if (flagRIL==1)
        {
          r_tr(0,0)=1/(1+2*r02);
          r_tr(0,1)=2*r02/(1+2*r02);
          r_tr(1,0)=2*r02/(1+2*r02);
          r_tr(1,1)=1/(1+2*r02);
        }    
        NumericMatrix az1rt(1,2);
        for (int i = 0; i < 1; i++)            
        {
          for (int j = 0; j < 2; j++)  
          {  
            az1rt(i,j) = 0;  
            for (int k = 0; k < 2; k++){  
              az1rt(i,j) += az1(i,k) * r_tr(k,j);  
            }  
          }  
        }
        az1=az1rt;  
      }
      if ((abs(itt+1-abc))>1e-6)
      {
        NumericMatrix azdd1(1,2);
        for (int i = 0; i < 1; i++)            
        {
          for (int j = 0; j < 2; j++)  
          {  
            azdd1(i,j) = 0;  
            for (int k = 0; k < 2; k++){  
              azdd1(i,j) += az1(i,k) * dd1(k,j);  
            }  
          }  
        }
        az1=azdd1;
        if (flagRIL==0){
          r_tr(0,0)=1-r0t(0,(itt-uu1+1));
          r_tr(0,1)=r0t(0,(itt-uu1+1));
          r_tr(1,0)=r0t(0,(itt-uu1+1));
          r_tr(1,1)=1-r0t(0,(itt-uu1+1));
        }
        if (flagRIL==1)
        {
          r_tr(0,0)=1/(1+2*r0t(0,(itt-uu1+1)));
          r_tr(0,1)=2*r0t(0,(itt-uu1+1))/(1+2*r0t(0,(itt-uu1+1)));
          r_tr(1,0)=2*r0t(0,(itt-uu1+1))/(1+2*r0t(0,(itt-uu1+1)));
          r_tr(1,1)=1/(1+2*r0t(0,(itt-uu1+1)));
        }
        NumericMatrix azrt2(1,2);
        for (int i = 0; i < 1; i++)            
        {
          for (int j = 0; j < 2; j++)  
          {  
            azrt2(i,j) = 0;  
            for (int k = 0; k < 2; k++){  
              azrt2(i,j) += az1(i,k) * r_tr(k,j);  
            }  
          }  
        }
        az1=azrt2;
      }
    }
    dd1(0,0)=1.0;
    dd1(0,1)=0.0;
    dd1(1,0)=0.0;
    dd1(1,1)=1.0;
    if (geno(j_1-1,uu2-1)==gg1){dd1(1,1)=0;}
    if (geno(j_1-1,uu2-1)==gg2){dd1(0,0)=0;}
    if (geno(j_1-1,uu2-1)==gg0){
      for (int id = 0;id < 2;id++){
        for (int jd = 0;jd < 2;jd++){ 
          dd1(id,jd)=dd1(id,jd)*0.5;
        }
      }
    }
    for (int i = 0; i < 1; i++)            
    {
      for (int j = 0; j < 2; j++)  
      {  
        az1dd1(i,j) = 0;  
        for (int k = 0; k < 2; k++){  
          az1dd1(i,j) += az1(i,k) * dd1(k,j);  
        }  
      }  
    }
    az1=az1dd1;
    vaz1(0,0)=1.0;
    vaz1(1,0)=1.0;
    for (int i = 0; i < 1; i++)            
    {
      for (int j = 0; j < 1; j++)  
      {  
        az1va(i,j) = 0;  
        for (int k = 0; k < 2; k++){  
          az1va(i,j) += az1(i,k) * vaz1(k,j);  
        }  
      }  
    }
    ac4(0,1)=az1va(0,0);
    cn00(0,0)=ac4(0,0)/(ac4(0,0)+ac4(0,1));
    cn00(0,1)=ac4(0,1)/(ac4(0,0)+ac4(0,1));
    markerpart(0,ire)=cn00(0,0)-cn00(0,1);
  } 
  return(markerpart);
}
