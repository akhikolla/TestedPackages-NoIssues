/********** C++ function: printPercMsgs **********/

/* Prints percentage-type progress messages. */

/* Last changed: 16 SEP 2020 */

#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]

int printPercMsgs(int msgCode,int loopSize,int iLoop,int percCnt)
{
   /* Declare all non-input variables: */

   int incmnt;
   double currPerc;

   /* Organise percentage message accoding 
      to the specified message code. */
   
   incmnt = 1;
   currPerc = 100*iLoop/loopSize;
   if (currPerc>=percCnt)
   {
      if (msgCode>0)
      { 
         if (percCnt>0)
         {
            if (percCnt<10) 
               {Rcout << percCnt << ",";};              
              
            if ((percCnt>=10)&(percCnt<=99)) 
               {Rcout << percCnt << ",";};              
        
            if (percCnt>=100) 
               {Rcout << percCnt << ".\n";};     

            if (msgCode==2)
            {
               if ((percCnt==20)||(percCnt==37)||(percCnt==54)
                   ||(percCnt==71)||(percCnt==88))
                  {Rcout << "\n   ";};     
            };            
         }
         if (msgCode==1) 
         {
            if (percCnt<10) {incmnt = 1;};
            if (percCnt>=10) {incmnt = 10;};
         }
         if (msgCode==2) {incmnt = 1;};
         if (msgCode==3) {incmnt = 10;};
         percCnt = percCnt + incmnt;
      }
   }

   return percCnt;
}

/************ End of printPercMsgs ************/

