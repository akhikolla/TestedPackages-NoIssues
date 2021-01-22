//
//  huge.cpp
//
//
//  Created by 张力翔 on 2018/9/4.
//
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <Rcpp.h>
using namespace Rcpp;


extern int VERSION2;

typedef struct
{
    int id;
    float value;
} SORT_FLOAT;

typedef struct
{
    int n;
    int *id;
} CLink;

typedef struct
{
    int id;
    int value;
} SORT_INT;

#define EPS 1.0e-3

void simp1(double **a, int mm, int *ll, int nll, int iabf,
           int *kp, double *bmax)
{
    int k;
    double test;

    if (nll <=0)
        *bmax = 0.0;
    else {
        *kp = ll[1];
        *bmax=a[mm+1][*kp+1];
        for (k=2; k<=nll; k++) {
            if (iabf == 0)
                test = a[mm+1][ll[k]+1]-(*bmax);
            else
                test = fabs(a[mm+1][ll[k]+1])-fabs(*bmax);
            if (test > 0.0) {
                *bmax = a[mm+1][ll[k]+1];
                *kp = ll[k];
            }
        }
    }
}

/*--------------------------------------------------------*/
/* Locate a pivot element, taking degeneracy into account */
/*--------------------------------------------------------*/
void simp2(double **a, int m, int n, int *ip, int kp)
{
    int k,i;
    double qp, q0, q, q1;

    *ip = 0;
    for (i=1; i<=m; i++)
        if (a[i+1][kp+1]<-EPS) break; /* Any possible pivots? */

    if (i>m) return;
    q1 = -a[i+1][1]/a[i+1][kp+1];
    *ip = i;
    for (i=*ip+1; i<=m; i++) {
        if (a[i+1][kp+1] < -EPS) {
            q = -a[i+1][1]/a[i+1][kp+1];
            if (q<q1) {
                *ip = i;
                q1 = q;
            }
            else {
                if (q==q1) {
                    for (k=1; k<=n; k++) {
                        qp = -a[*ip+1][k+1]/a[*ip+1][kp+1];
                        q0 = -a[i+1][k+1]/a[i+1][kp+1];
                        if (q0 != qp) break;
                    }
                    if (q0<qp) *ip = i;
                }
            }
        }
    }
}

/*-------------------------------------------------------------------*/
/* Matrix operations to exchange a left-hand and right-hand variable */
/*-------------------------------------------------------------------*/
void simp3(double **a, int i1, int k1, int ip, int kp)
{
    int kk, ii;
    double piv;

    piv = 1.0/a[ip+1][kp+1];
    for (ii=1; ii<=i1+1; ii++)
        if (ii-1 != ip ) {
            a[ii][kp+1] *= piv;
            for (kk=1; kk<=k1+1; kk++)
                if (kk-1!=kp)
                    a[ii][kk] -= a[ip+1][kk]*a[ii][kp+1];
        }

    for (kk=1; kk<=k1+1; kk++)
        if (kk-1 != kp) a[ip+1][kk] *= -piv;
    a[ip+1][kp+1] = piv;
}

/*-----------------------------------------------------------------------*/
/* Simplex method for linear programming. Input paramenters a, m, n, mp, */
/* np, m1, m2, and m3, and output parameters a, icase, izrov, and iposv. */
/*-----------------------------------------------------------------------*/

void simplx(double **a, int m, int n, int m1, int m2, int m3, int *icase,
            int *izrov, int *iposv)
{
    int i, ip, is, k, kh, kp, nl1;
    int *l1, *l3;
    double q1, bmax;

    if (m != (m1+m2+m3)) {
        Rcpp::stop("Bad input constraint counts in simplx");    }

    l1=(int *)calloc(n+2,sizeof(int));
    l3=(int *)calloc(m+1,sizeof(int));

    if (l1==NULL || l3==NULL) {
        Rcpp::stop("Can't allocate space in simplx");
    }

    nl1 = n;
    for (k=1; k<=n; k++) l1[k] = izrov[k] = k;


    /* Initialize index list of columns admissible for exchange, and make
     all variables initially right-hand. */

    for (i=1; i<=m; i++) {
        if (a[i+1][1] < 0.0) {   /* constants b_i must be nonnegative */
            Rcpp::stop("Bad input tableau in simplx");
        }
        iposv[i] = n+i;
    }

    if (m2+m3) {
        for (i=1; i<=m2; i++) l3[i] = 1;
        for (k=1; k<=n+1; k++) {
            q1 = 0.0;
            for (i=m1+1; i<=m; i++) q1 += a[i+1][k];
            a[m+2][k] = -q1;
        }
        for (; ;) {
            simp1(a, m+1, l1, nl1, 0, &kp, &bmax);
            if (bmax <= EPS && a[m+2][1] < -EPS) {
                *icase = -1;
                free(l3);
                free(l1);
                return;
            }
            else {
                if (bmax <=EPS && a[m+2][1]<=EPS) {
                    for (ip=m1+m2+1; ip <=m; ip++) {
                        if (iposv[ip] == (ip+n)) {
                            simp1(a, ip, l1, nl1, 1, &kp, &bmax);
                            if (bmax > EPS)
                                goto one;
                        }
                    }
                    for (i=m1+1; i<=m1+m2; i++)
                        if (l3[i-m1] == 1)
                            for (k=1; k<=n+1; k++)
                                a[i+1][k] = -a[i+1][k];
                    break;
                }
            }

            simp2(a, m, n, &ip, kp);
            if (ip==0) {
                *icase = -1;
                free(l1); free(l3);
                return;
            }
        one: simp3(a, m+1, n, ip, kp);
            if (iposv[ip] >= (n+m1+m2+1)) {
                for (k=1; k<=nl1; k++)
                    if (l1[k] == kp) break;
                --nl1;
                for (is=k; is <=nl1; is++) l1[is]=l1[is+1];
            }
            else {
                kh = iposv[ip]-m1-n;
                if (kh >= 1 && l3[kh]) {
                    l3[kh] = 0;
                    ++a[m+2][kp+1];
                    for (i=1; i<=m+2; i++)
                        a[i][kp+1]= -a[i][kp+1];
                }
            }
            is = izrov[kp];
            izrov[kp]=iposv[ip];
            iposv[ip] = is;
        }
    }

    for (; ; ) {
        simp1(a, 0, l1, nl1, 0, &kp, &bmax);
        if (bmax <= EPS) {
            *icase = 0;
            free (l1); free(l3);
            return;
        }
        simp2(a, m, n, &ip, kp);
        if (ip == 0) {
            *icase = 1;
            free(l1); free(l3);
            return;
        }
        simp3(a,m,n,ip,kp);
        is = izrov[kp];
        izrov[kp] = iposv[ip];
        iposv[ip] = is;
    }
}


//simplex.c
extern void simplx(double **, int, int, int, int, int, int *, int *, int *);

int VERSION2=1;



/*--------------------------------------------------------*/
/* Compute distance between all the pairs of clusters     */
/* based on membership of all the data samples.           */
/*--------------------------------------------------------*/

float dist2cls(int *cls1, int *cls2, int len, int id1, int id2)
{
    int i;
    float v1;

    for (i=0,v1=0.0; i<len;i++){
        if ((cls1[i]==id1 && cls2[i]!=id2)||(cls1[i]!=id1 && cls2[i]==id2))
            v1+=1.0;
    }
    return(v1);
}

float dist2cls_normalized(int *cls1, int *cls2, int len, int id1, int id2)
{
    int i;
    float v1,v2,v3,r1,r2;

    v1=v2=v3=0.0;
    for (i=0; i<len;i++){
        if (cls1[i]==id1 && cls2[i]!=id2)
            v1+=1.0; //in cluster 1 but not in cluster 2
        if (cls1[i]!=id1 && cls2[i]==id2)
            v2+=1.0; //in cluster 2 but not in cluster 1
        if (cls1[i]==id1 && cls2[i]==id2)
            v3+=1.0; //in both clusters
    }

    if (v1+v3==0.0) r1=1.0; else r1=v1/(v1+v3);//ratio of subtracted set
    if (v2+v3==0.0) r2=1.0; else r2=v2/(v2+v3);

    //return((r1+r2)/2.0);
    //Revised on Nov 7, 2017 to Jacaad index
    if (v3+v1+v2==0.0) r1=1.0; else r1=(v1+v2)/(v3+v1+v2);
    return(r1);
}

//Revised from dist2cls_normalized() above
//Compute Jacaard for two clusters represented by set of point ids
float Jacaard_pts(int *ptid1, int len1, int *ptid2, int len2)
{
    int i;
    float v1,v2,v3,r1;
    int *cls1, *cls2;
    int id1, id2;
    int len;

    len=0;
    for (i=0;i<len1;i++) if (len<ptid1[i]) len=ptid1[i];
    for (i=0;i<len2;i++) if (len<ptid2[i]) len=ptid2[i];
    len++;

    cls1=(int *)calloc(len,sizeof(int));
    cls2=(int *)calloc(len,sizeof(int));
    for (i=0;i<len;i++) cls1[i]=cls2[i]=0;
    for (i=0;i<len1;i++) cls1[ptid1[i]]=1;
    for (i=0;i<len2;i++) cls2[ptid2[i]]=1;

    id1=1; id2=1;
    v1=v2=v3=0.0;
    for (i=0; i<len;i++){
        if (cls1[i]==id1 && cls2[i]!=id2)
            v1+=1.0; //in cluster 1 but not in cluster 2
        if (cls1[i]!=id1 && cls2[i]==id2)
            v2+=1.0; //in cluster 2 but not in cluster 1
        if (cls1[i]==id1 && cls2[i]==id2)
            v3+=1.0; //in both clusters
    }

    //Jacaard index
    if (v3+v1+v2==0.0) r1=1.0; else r1=(v1+v2)/(v3+v1+v2);

    free(cls1);
    free(cls2);

    return(r1);

}



void allpairs(int *cls1, int *cls2, int len, int n1, int n2, float *distmat)
//n1 is the #cluster in first clustering result, n2 is that for the second result
//dist[n1*n2], dist[0:n2-1] is the distance between cluster 0 in result 1 and cluster
//0, ..., n2-1 in result 2, so on so forth
{
    int i,j;

    if (!VERSION2) {
        for (i=0; i<n1; i++) {
            for (j=0; j<n2; j++)
                distmat[i*n2+j]=dist2cls(cls1,cls2,len,i,j);
        }
    } else {//Normalized distance with respect to each cluster's size, Jacaard index used
        for (i=0; i<n1; i++) {
            for (j=0; j<n2; j++)
                distmat[i*n2+j]=dist2cls_normalized(cls1,cls2,len,i,j);
        }
    }
}


/*--------------------------------------------------------*/
/* IRM matching scheme.                                   */
/*--------------------------------------------------------*/
float match_fast(float *dist, float *p1_in, float *p2_in,
                 int num1, int num2, float *wt)
/* wt[num1*num2] stores the computed weights by optimization */
{
    int i,j, ii,jj;
    float minval;
    float *p1, *p2;
    int sum1, sum2;
    float res;
    float TINY=1.0e-8;

    p1=(float *)calloc(num1,sizeof(float));
    p2=(float *)calloc(num2,sizeof(float));

    for (i=0; i<num1; i++) p1[i]=p1_in[i];
    for (i=0; i<num2; i++) p2[i]=p2_in[i];

    for (i=0;i<num1*num2;i++) wt[i]=0.0;

    sum1=sum2=0;

    while(sum1<num1 && sum2<num2) {
        minval=HUGE_VAL;
        ii=0, jj=0;
        for (i=0; i<num1; i++) {
            if (p1[i]<TINY)
                continue;
            for (j=0; j<num2; j++) {
                if (p2[j]<TINY)
                    continue;
                if (dist[i*num2+j] < minval)
                {
                    ii=i;
                    jj=j;
                    minval = dist[i*num2+j];
                }
            }
        }

        if (p1[ii]<=p2[jj]) {
            wt[ii*num2+jj] = p1[ii];
            p2[jj] -= p1[ii];
            p1[ii]=0.0;
            sum1++;
            if (p2[jj]<TINY) sum2++;
        }
        else {
            wt[ii*num2+jj] = p2[jj];
            p1[ii] -= p2[jj];
            p2[jj]= 0.0;
            sum2++;
            if (p1[ii]<TINY) sum1++;
        }
    }

    for (i=0,res=0.0;i<num1*num2;i++) res+=(wt[i]*dist[i]);

    free(p1);
    free(p2);

    return(res);

}

/*--------------------------------------------------------*/
/* Compute E(||X-Y||^p) by mallows distance matching.     */
/* The distance matrix ||X-Y||^p is assumed given by      */
/* *dist.                                                 */
/*--------------------------------------------------------*/
float match(float *dist, float *p1, float *p2, int num1, int num2, float *wt)
/* wt[num1*num2] stores the computed weights by optimization */
{
    int i,j,k,m,n;
    int nconstraint, nvar;
    int icase;
    float res,v1;
    double **a;
    int *iposv,*izrov;

    nvar = num1*num2;
    nconstraint = num1+num2;

    /* The unknowns are elements of wt[num1*num2] */

    /*-------------------------------------------*/
    /* Allocate space                            */
    /*-------------------------------------------*/

    m=(num1>num2)?num1:num2;
    a=(double **)calloc(m*2+3,sizeof(double *));
    for (i=0; i<m*2+3; i++) a[i]=(double *)calloc(m*m+2,sizeof(double));

    iposv=(int *)calloc(m*2+1,sizeof(int));
    izrov=(int *)calloc(m*m+1,sizeof(int));

    /*-------------------------------------------*/
    /* Use a to store all the constrains in the  */
    /* linear programming problem.               */
    /*-------------------------------------------*/

    for (i=0; i<nconstraint+3; i++)
        for (j=0; j<nvar+2; j++)
            a[i][j]=0.0;

    for (i=2, k=0; k<num1; i++, k++) {
        a[i][1] = p1[k];
        if (a[i][1]<0.0) a[i][1]=0.0;

        m = 2+(k+1)*num2;
        for (j=2+k*num2; j<m; j++)
            a[i][j] = -1.0;
    }

    for (i=num1+2,k=0; k<num2; i++, k++)
    {
        a[i][1] = p2[k];
        if (a[i][1]<0.0) a[i][1]=0.0;

        for (j=2+k, n=0; n<num1; j+=num2, n++)
            a[i][j] = -1.0;
    }

    a[1][1] = 0.0;
    for (i=0; i<nvar; i++)
        a[1][i+2] = -dist[i];  // negative sign due to we want to minimize
    // while simplx is for maximize

    simplx(a, nconstraint, nvar, 0, 0, nconstraint, &icase, izrov, iposv);

    if (icase!=0) {
        // icase=1, unbounded function, icase=-1, no feasible solution
        // use IRM distance
        Rcpp::warning("Warning: Mallows distance replaced by IRM");
        res=match_fast(dist, p1, p2, num1, num2,wt);
    }
    else {
        res=-a[1][1];

        for (i=0;i<num1*num2;i++) wt[i]=0.0;
        for (j=0; j<nconstraint; j++) {
            if (iposv[j+1]-1<nvar)
                wt[iposv[j+1]-1]=a[j+2][1];
        }
    }


    for (i=0,v1=0.0;i<num1*num2;i++) v1+=wt[i]*dist[i];
    //fprintf(stderr, "res=%e, v1=%e\n",res,v1);

    m=(num1>num2)?num1:num2;
    for (i=0; i<m*2+3; i++) free(a[i]);
    free(a);
    free(iposv);
    free(izrov);

    return(res);
}

/*--------------------------------------------------------*/
/* Compute Mallows distance for two bags of vectors.      */
/*--------------------------------------------------------*/
float alignclusters(int *cls1, int *cls2, int len, int n1, int n2, float *wt)
{
    int i,m;
    float *dist,mdist;
    float *p1, *p2, v1,v2;

    //Determine n1 or n2 automatically
    if (n1<=0) {
        n1=0;
        for (i=0;i<len;i++) {if (cls1[i]>n1) n1=cls1[i];}
        n1++;
    }
    if (n2<=0) {
        n2=0;
        for (i=0;i<len;i++) {if (cls2[i]>n2) n2=cls2[i];}
        n2++;
    }

    //compute all pairwise distances
    dist = (float *)calloc(n1*n2,sizeof(float));
    allpairs(cls1,cls2,len,n1,n2,dist);

    //compute weight for each cluster
    p1=(float *)calloc(n1,sizeof(float));
    p2=(float *)calloc(n2,sizeof(float));

    for (i=0,v1=v2=0.0;i<len;i++) {
        if (cls1[i]>=0) { p1[cls1[i]]+=1.0; v1+=1.0;}
        if (cls2[i]>=0) { p2[cls2[i]]+=1.0; v2+=1.0;}
    }

    //showinfo(dist,n1,n2,p1,p2);

    for (i=0;i<n1;i++) p1[i]/=v1;
    for (i=0;i<n2;i++) p2[i]/=v2;

    //Optimal transport matching (Wasserstein matching)
    mdist=match(dist, p1, p2, n1, n2,wt);

    //Forcefully set negative values in wt[] into 0.0
    for (i=0,m=n1*n2; i<m;i++){  if (wt[i]<0.0) wt[i]=0.0;}

    free(dist);
    free(p1); free(p2);

    return(mdist);
}

//The first len elements of cls[] are the reference clustering result
//ns is # bootstrap samples, including the reference clustering result
//len is the number of points in the data set
void align(int *cls, int ns, int len, float **wt, int **clsct_pt, float **dc_pt, int equalcls)
{
    int i,j,m;
    float *dc;
    int **clsarray, *clsct;

    if (ns<=1) {
        Rcpp::stop("Wrong input: number of clustering results < 2");
    }

    clsarray=(int **)calloc(ns,sizeof(int *));
    for (i=0;i<ns;i++) clsarray[i]=cls+i*len;

    clsct=(int *)calloc(ns,sizeof(int));
    for (i=0;i<ns;i++){
        clsct[i]=0;
        for (j=0;j<len;j++) {if (clsarray[i][j]>clsct[i]) clsct[i]=clsarray[i][j];}
        clsct[i]++;
    }

    if (equalcls) {//force to have same number of clusters for each sample
        for (i=0,m=0;i<ns;i++){ if (m<clsct[i]) m=clsct[i];}
        if (clsct[0]<m) {
            Rcpp::warning("The reference clustering has empty cluster");
        }
        for (i=0;i<ns;i++){ clsct[i]=m;}
    }

    dc=(float *)calloc(ns,sizeof(float));
    dc[0]=0.0; //distance to itself
    for (i=1,m=0;i<ns;i++) m+=clsct[i];
    m*=clsct[0];

    *wt=(float *)calloc(m,sizeof(float));
    for (i=1,m=0;i<ns;i++) {
        dc[i]=alignclusters(clsarray[i], clsarray[0], len, clsct[i], clsct[0], *wt+m);
        m+=clsct[i]*clsct[0];
    }

    *clsct_pt=clsct;
    *dc_pt=dc;
    free(clsarray);
    return;
}

void convertcode(int *cls, int len, int minsz)
{
    int i,m;
    int nc,*clsct,*code;

    for (i=0,nc=0;i<len;i++) {
        if (cls[i]>nc) nc=cls[i];
    }
    nc++;

    clsct=(int *)calloc(nc,sizeof(int));
    code=(int *)calloc(nc,sizeof(int));
    for (i=0;i<nc;i++) clsct[i]=0;
    for (i=0;i<len;i++) clsct[cls[i]]++;

    for (i=0,m=0;i<nc;i++) {
        if (clsct[i]<minsz) {code[i]=-1;}
        else {code[i]=m; m++;}
    }

    for (i=0;i<len;i++) cls[i]=code[cls[i]];

    free(clsct);
    free(code);
}

//wt[n1*n2], 1D array for matrix wt[n1][n2], each column corresponds
// to one cluster in the second result, and each row corresponds to
// one cluster in the first result
//reference cluster result is the seconde cluster result
//*code: 0 for match, 1 for split (a cluster in the second result is
//split into mulitple clusters in the first result), 2 for merge, 3 for others
//*nf: 1 for match, #clusters split into for split, #clusters needed to merge
//for merge
//code[n2], nf[n2] have been allocated with space
float maxentry(float *wt, int n1,int *id)
{
    int j;
    float v3;

    *id=0;
    v3=wt[0];
    for (j=1;j<n1;j++) {
        if (wt[j]>v3) {
            v3=wt[j];
            *id=j;
        }
    }
    return(v3);
}

//*nf is the number of clusters covered by the cluster in the reference
//return the value of the covered percentage of the cluster in the reference by
//all the clusters in the compared clustering that are covered by that cluster
//*maxv is the maxmum coverage by those clusters in the compared clustering that
//covered
float covercmp(float *wtcmp, float *wtref, int n1, int n2, int *nf, float *maxv,
               int *maxid, float thred, float *coverage)
{
    int j,m,n;
    float v1,v2;

    //*maxv=v2 is the maximum of wtref[] of those clusters covered in the compared clustering
    //returned value v1 is the sum of such wtref[]'s. Hence v1>=v2.
    n=0;
    m=0;
    v1=v2=0.0;
    for (j=0;j<n1;j++) {
        if (wtcmp[j]>=thred) {
            if (coverage!=NULL) coverage[j]=wtref[j];
            n++;
            v1+=wtref[j];
            if (wtref[j]>v2) {
                v2=wtref[j];
                m=j;
            }
        } else {
            if (coverage!=NULL) coverage[j]=-1.0;
        }
    }

    *nf=n;
    *maxv=v2;
    *maxid=m;
    return(v1);
}



//The matching,splitting, merging matrix res[n1*n2] cannot be of arbitrary pattern.
//It has to satisfy certain clear patterns. This is a parity check.
void paritycheck(float *res, int n1, int n2)
{
    int i,j;
    int m1,m2,m3;

    for (i=0;i<n2;i++) {
        m1=m2=m3=0;
        for (j=0;j<n1;j++) {
            if (res[j*n2+i]<0.0) m1++;
            if (res[j*n2+i]>=0.0 && res[j*n2+i]<=1.0) m2++;
            if (res[j*n2+i]>=2.0 && res[j*n2+i]<=3.0) m3++;
        }
        if (m1+m2+m3<n1) {
            Rcpp::warning("m1+m2+m3<n1");
        }

        if (m3>1) {
            Rcpp::warning("Merge to more than 1");
        }
        else {
            if (m3==1) {
                if (m3+m1<n1)
                    Rcpp::warning("m3+m1<n1");
            }
            else {//m3==0 case
                if (m1+m2<n1)
                    Rcpp::warning("m2+m1<n1");
            }
        }
    }
}

//Modified from assess2. Differ by having an extra output of the same dimension as
//wt[n1*n2], stores matching, splitting, and merging results res[n1*n2]
void assess2(float *wt, float *res, int n1, int n2, int *code, int *nf, float thred)
{
    int i,j,m,ii;
    float *wtcol, *wtrow; //column and row wise normalization of wt
    float *wtcmp, *wtref;
    float v1,v2,v3,v4,v5;
    int maxid;
    float *coverage;

    wtcol=(float *)calloc(n1*n2,sizeof(float));
    wtrow=(float *)calloc(n1*n2,sizeof(float));
    wtref=(float *)calloc(n1,sizeof(float));
    wtcmp=(float *)calloc(n1,sizeof(float));
    coverage=(float *)calloc(n1>n2 ? n1:n2, sizeof(float));

    //row-wise normalization
    for (i=0;i<n1;i++) {
        for (j=0,v1=0.0;j<n2;j++)
            v1+=wt[i*n2+j];
        if (v1>0.0) {
            for (j=0;j<n2;j++)
                wtrow[i*n2+j]=wt[i*n2+j]/v1;
        }
        else {//treating empty clusters
            for (j=0;j<n2;j++)
                wtrow[i*n2+j]=0.0;
        }
    }

    //column-wise normalization
    for (i=0;i<n2;i++) {
        for (j=0,v2=0.0;j<n1;j++)
            v2+=wt[j*n2+i];
        if (v2>0.0) {
            for (j=0;j<n1;j++)
                wtcol[j*n2+i]=wt[j*n2+i]/v2;
        }
        else {//treating empty clusters
            for (j=0;j<n1;j++)
                wtcol[j*n2+i]=0.0;
        }
    }

    for (ii=0;ii<n2;ii++) {
        for (j=0;j<n1;j++) wtref[j]=wtcol[j*n2+ii];
        for (j=0;j<n1;j++) wtcmp[j]=wtrow[j*n2+ii];
        v1=covercmp(wtcmp, wtref, n1, n2, nf+ii, &v2,&maxid,thred,coverage);
        for (j=0;j<n1;j++) { res[j*n2+ii]=coverage[j]; }

        if (v2>=thred) {//match
            code[ii]=0;
        }
        else {
            if (v1>=thred) {//split
                code[ii]=1;
            } else {//possible merge or none of the three cases
                //It's possible that several components in reference merge to a component
                //in the bootstrap sample and one of merged members also matches with the
                //component in the bootstrap sample
                v3=maxentry(wtref,n1,&m);
                for (j=0;j<n1;j++) { res[j*n2+ii]=-1.0; }
                if (v3>=thred) {//possible merge
                    //m is the cluster which the reference clusters possibly merge to
                    v4=covercmp(wtcol+m*n2,wtrow+m*n2,n2,n1,nf+ii,&v5,&maxid, thred,coverage);

                    if (v4>=thred) {//merge
                        code[ii]=2;
                        res[m*n2+ii]=2.0+coverage[ii];
                        if (coverage[ii]<0.0) { Rcpp::warning("Paradox in assess2()");}
                    }
                    else {
                        code[ii]=3;
                        nf[ii]=0;
                    }
                }
                else {
                    code[ii]=3;
                    nf[ii]=0;
                }
            }
        }
    }

    //Purely for testing purpose
    //Cases for a cluster in the reference: Pure match, impure match, split, merge, none-above
    paritycheck(res,n1,n2);

    free(wtcol);
    free(wtrow);
    free(wtcmp);
    free(wtref);
    free(coverage);
    return;
}

//Based on the matching weight matrix, generate the match,split,merge, etc.
//status for all the bootstrap samples
//Output: res, codect, nfave
//Basically, calls assess2() for all the bootstrap samples and summarize
//some information for each cluster in the reference, that is, codect, nfave
void MatchSplit(float *wt, float *res, int *numcls, int nbs,
                int **codect, float **nfave, float thred)
{
    int i,j,k,m, n0;
    int *code, *nf;

    n0=numcls[0];
    if (thred<=0.5){
        Rcpp::warning("Coverage threshold is too small");
    }

    code=(int *)calloc(n0,sizeof(int));
    nf=(int *)calloc(n0,sizeof(int));

    for (j=0;j<n0;j++) {
        for (k=0;k<4;k++) {
            codect[j][k]=0;
            nfave[j][k]=0.0;
        }
    }

    //Matching, split, merge assessment for each bootstrap sample
    for (i=1,m=0;i<nbs;i++){
        assess2(wt+m*numcls[0],res+m*numcls[0],numcls[i],n0,code,nf,thred);
        m+=numcls[i];

        for (j=0;j<n0;j++) {
            codect[j][code[j]]++;
            nfave[j][code[j]]+=(float)nf[j];
        }
    }

    //Summarize the result for each cluster in the reference result
    for (j=0;j<n0;j++) {
        for (k=0;k<4;k++) {
            if (codect[j][k]>0)
                nfave[j][k]/=(float)codect[j][k];
        }
    }

    free(code);
    free(nf);
    return;
}

//-------- For computing confidence set for all the matched/split clusters for --
//-------- each cluster in the reference clustering with #clusters=n2
//??? Insert confset2.c here
static int CompFcn(const void *a,const void *b)
{
    SORT_INT c=*(SORT_INT*)a;
    SORT_INT d=*(SORT_INT*)b;
    if (c.value > d.value)
        return (1);
    if (c.value < d.value)
        return (-1);
    return (0);
}

void SortInt(int *org, int *buf, int *invid, int sz)
{
    int j;
    SORT_INT *score;

    score=(SORT_INT *)calloc(sz,sizeof(SORT_INT));
    if (score==NULL) {
        Rcpp::stop("Unable to allocate space in SortInt");
    }
    for (j=0;j<sz;j++) {
        score[j].id=j;
        score[j].value=org[j];
    }
    qsort(score, sz, sizeof(SORT_INT), CompFcn);

    for (j=0;j<sz;j++) {
        buf[j]=org[score[j].id];
        invid[j]=score[j].id; //store the original id for the new order
    }

    free(score);
}

//Sort the ids in each clist[] so that they are in ascending order
void Sortcls(CLink *clist, int ns)
{
    int i,j,m;
    int *buf, *invid;

    m=0;
    for (i=0;i<ns;i++) {
        if (m<clist[i].n) m=clist[i].n;
    }

    buf=(int *)calloc(m,sizeof(int));
    invid=(int *)calloc(m,sizeof(int));
    for (i=0;i<ns;i++) {
        SortInt(clist[i].id, buf, invid, clist[i].n);
        for (j=0;j<clist[i].n;j++)
            clist[i].id[j]=buf[j];
    }
    free(buf);
}


//cls[0] has space allocated
void NewCLink(CLink *cls, int np)
{
    cls->n=np;
    cls->id=(int *)calloc(np,sizeof(int));
}

void FreeCLink(CLink *cls)
{
    if (cls->n >0) free(cls->id);
}



//Generate a list of numbers from 0 to *nids-1 for the original IDs
//in clist[].
void MapIds(CLink *clist, int ns, int *maxid, int *nids, int **id2num, int **num2id)
{
    int i,j,k,m;

    m=0;
    for (i=0;i<ns;i++)
        for (j=0;j<clist[i].n;j++)
            if (clist[i].id[j]>m)
                m=clist[i].id[j];
    m++;
    *maxid=m; //maximum value of id plus 1

    *id2num=(int *)calloc(m,sizeof(int));
    for (i=0;i<m;i++) (*id2num)[i]=0;

    for (i=0;i<ns;i++)
        for (j=0;j<clist[i].n;j++)
            (*id2num)[clist[i].id[j]]++;

    k=0;
    for (i=0;i<m;i++) {
        if ((*id2num)[i]==0) (*id2num)[i]=-1; //non-existing ids
        else {
            (*id2num)[i]=k;
            k++;
        }
    }

    *nids=k;
    *num2id=(int *)calloc(k,sizeof(int));
    for (i=0;i<m;i++)
        if ((*id2num)[i]>=0) {
            (*num2id)[(*id2num)[i]]=i;
        }
}

//For a given point with id, compute the number of clusters that
//contain this point. Whether a cluster should be considered is indicated
//by keepcls[ns].
int ClusterInclude(CLink *clist, int ns, unsigned char *keepcls, int id, unsigned char *touched)
{
    int i,j,n;

    n=0;
    for (i=0;i<ns;i++) {
        touched[i]=0;
        if (keepcls[i]==0) continue;
        for (j=0;j<clist[i].n;j++) {
            if (clist[i].id[j]==id) {
                n++;
                touched[i]=1;
                break;
            }
            else {
                if (clist[i].id[j]>id) break; //We can do this assuming the ids are sorted
            }
        }
    }

    return(n);
}

//Find the confidence set. The only output is pts[nids], keepcls[ns]
//The space for pts[] and keepcls[] has been allocated.
//alpha is the percentage of clusters that can be excluded at most.
void ConfidenceSet(CLink *clist, int ns, int nids, int *id2num, int *num2id,
                   unsigned char *pts, unsigned char *keepcls, float alpha)
{
    int i,j,m;
    int nexclude, ncv;
    unsigned char *touched, *buf;
    int *p, mv;

    p=(int *)calloc(nids,sizeof(int));
    touched =(unsigned char *)calloc(ns,sizeof(char));
    buf =(unsigned char *)calloc(ns,sizeof(char));

    nexclude=(int)(alpha*(float)ns);
    ncv=ns;

    for (i=0;i<nids;i++) pts[i]=1; //indicator for the inclusion of each point
    for (i=0;i<ns;i++) keepcls[i]=1; //indicator for each cluster being covered

    while (ncv>ns-nexclude) {
        //For each existing point and existing cluster, compute p[]
        for (i=0;i<nids;i++) p[i]=0;
        m=-1; mv=ns+1;
        for (i=0;i<nids;i++) {
            if (pts[i]==0) continue; //already excluded
            p[i]=ClusterInclude(clist, ns, keepcls, num2id[i], buf);
            if (p[i]<mv) {
                mv=p[i];
                m=i;
                for (j=0;j<ns;j++) touched[j]=buf[j];
            }
        }

        if (ncv-mv>=ns-nexclude) {
            //exclude points and clusters
            for (j=0;j<ns;j++)
                if (touched[j]) {
                    keepcls[j]=0;
                }
            ncv-=mv;

            for (i=0;i<nids;i++) pts[i]=0; //indicator for the inclusion of each point
            for (i=0;i<ns;i++) {
                if (keepcls[i]==0) continue;
                for (j=0;j<clist[i].n;j++) pts[id2num[clist[i].id[j]]]=1;
            }
        } else { //terminate the loop
            break;
        }
    }

    free(p);
    free(buf);
    free(touched);
}


//  int nbs,len
//  float alpha=0.1, v1,v2,v3;
//  CLink *clist;
void confset(CLink *clist, int nbs, float alpha, int **confpts, int *npts, unsigned char **keepcls_pt, float **cvp_pt, float *tightness, float *coverage)
{
    int i,j,k,m,i1,i2;
    float v1,v2;

    /*----------------------------------------------------------------*/
    /*----------------- Read in data ---------------------------------*/
    /*----------------------------------------------------------------*/

    Sortcls(clist,nbs);

    /*------------------------------------------------------------*/
    /*- Find Confidence Set                                       */
    /*------------------------------------------------------------*/

    int nids, maxid, *id2num, *num2id;
    unsigned char *pts, *keepcls;
    float *cvp;

    MapIds(clist, nbs, &maxid, &nids, &id2num, &num2id);

    pts=(unsigned char *)calloc(nids, sizeof(char));
    keepcls=(unsigned char *)calloc(nbs, sizeof(char));
    ConfidenceSet(clist, nbs, nids, id2num, num2id, pts, keepcls, alpha);

    /*------------------------------------------------------------*/
    /*- Output the result                                         */
    /*------------------------------------------------------------*/
    for (i=0,m=0;i<nids;i++) if (pts[i]) m++;
    *npts=m;
    *confpts=(int *)calloc(m,sizeof(int));
    for (i=0,k=0;i<nids;i++)
        if (pts[i]) {
            (*confpts)[k]=num2id[i];
            k++;
        }

    //For included clusters, compute the percentage w.r.t. the confidence set
    //For excluded clusters, compute the percentage of points in each cluster
    //that are included in the confidence set
    cvp=(float *)calloc(nbs,sizeof(float));
    i1=i2=0;
    v1=v2=0.0;
    for (i=0;i<nbs;i++) {
        if (keepcls[i]) {
            cvp[i]=((float)clist[i].n)/((float)m);
            v1+=cvp[i];
            i1++;
        }
        else {
            for (j=0,k=0;j<clist[i].n;j++)
                if (pts[id2num[clist[i].id[j]]]) k++;
            cvp[i]=((float)k)/((float)clist[i].n); //percentage of covered points
            v2+=cvp[i];
            i2++;
        }
    }

    if (i1>0) v1/=(float)i1;
    if (i2>0) v2/=(float)i2; else v2=1.0; //Every is covered, no missing set

    *tightness=v1;
    *coverage=v2;

    *keepcls_pt=keepcls;
    *cvp_pt=cvp;
    free(pts);
    free(id2num);
    free(num2id);
}


//============================================
void MatchCluster(float *res, int *numcls, int nbs, float thred, int *cls, int len, CLink ***clist2, int **nsamples, int usesplit, int *matchID)
{
    int i,j,k,m,n0,n1,n2,ii,k2,kk;
    int *matched, nmatch;
    int m1,m2,m3;
    float *res_cur;
    float v1;


    for (i=0;i<len*nbs;i++) matchID[i]=-1;

    *nsamples=(int *)calloc(numcls[0],sizeof(int));
    *clist2=(CLink **)calloc(numcls[0], sizeof(CLink *));
    for (ii=1,n0=0;ii<nbs;ii++){  if (n0<numcls[ii]) n0=numcls[ii];}
    matched=(int *)calloc(n0,sizeof(int));

    n2=numcls[0];
    for (i=0;i<n2;i++) {//process each cluster in the reference
        for (ii=1,m=0,nmatch=0;ii<nbs;ii++){
            n1=numcls[ii];
            res_cur=res+m*n2;

            m1=m2=m3=0;
            v1=-1.0;
            for (j=0;j<n1;j++) {
                if (res_cur[j*n2+i]<0.0) m1++;
                if (res_cur[j*n2+i]>=0.0 && res_cur[j*n2+i]<=1.0) {
                    m2++;
                    if (v1<res_cur[j*n2+i]) v1=res_cur[j*n2+i];
                }
            }

            if (m2>=1){//match or split
                if (m2==1 || v1>=thred) {//match
                    nmatch++;
                }
                else { //split
                    if (usesplit) nmatch++;
                }
            }
            m+=numcls[ii];
        }

        //Reference cluster itself is always recorded
        (*clist2)[i]=(CLink *)calloc(nmatch+1,sizeof(CLink));
        (*nsamples)[i]=nmatch+1; //nsamples[i] is the number of samples used for confset for cluster i

        //Record the point ids for the reference cluster
        for (j=0,k=0;j<len;j++) {
            if (cls[j]==i){ k++; }
        }

        NewCLink((*clist2)[i], k);

        for (j=0,k=0;j<len;j++) {
            if (cls[j]==i){
                (*clist2)[i][0].id[k]=j; //point id
                matchID[j]=i;
                k++;
            }
        }

        //Record point ids for the matched clusters in the other partitions
        for (ii=1,m=0,kk=0;ii<nbs;ii++){
            n1=numcls[ii];
            res_cur=res+m*n2;
            for (j=0;j<n1;j++) matched[j]=0;

            m1=m2=m3=0;
            v1=-1.0; k2=0;
            for (j=0;j<n1;j++) {
                if (res_cur[j*n2+i]<0.0) m1++;
                if (res_cur[j*n2+i]>=0.0 && res_cur[j*n2+i]<=1.0) {
                    m2++;
                    matched[j]=1;
                    if (v1<res_cur[j*n2+i]) { v1=res_cur[j*n2+i]; k2=j;}
                }
            }

            if (usesplit!=1 && m2>1) {
                for (j=0;j<n1;j++) matched[j]=0;
                if (v1>=thred) //reset matched for impure match case
                    matched[k2]=1;//The single 1 in matched[]
            }

            //List point ids for the matched cluster(s)
            if (m2>=1 && (m2==1 || v1>=thred || usesplit)) {
                for (j=ii*len,k=0;j<ii*len+len;j++) {
                    if (matched[cls[j]]){
                        k++;
                    }
                }

                NewCLink((*clist2)[i]+kk+1, k);

                for (j=ii*len,k=0;j<ii*len+len;j++) {
                    if (matched[cls[j]]){
                        (*clist2)[i][kk+1].id[k]=j-ii*len; //point id
                        matchID[j]=i;
                        k++;
                    }
                }

                kk++; //index for nmatch
            }

            m+=numcls[ii];
        }
    }

    free(matched);

}

//confpts[numcls], npts[numcls], avetight[numcls], avecov[numcls], avejacaard[numcls], rinclude[numcls] have space allocated
//confpts[numcls] is an empty link to be allocated
//confpts[i] is a list of point ids that are in the confidence set
//The size of confpts[i] is output npts[i]
//avetight[i] and avecov[i] are outputs: average tightness and average coverage
//for the confidence set for the ith reference cluster
//avejacaard[i] is the average jacaard index between the reference
//and the other mathced clusters.
//rinclude[i] is the percentage of clusters fully covered by the confidence set
//csetdist[numcls*numcls] is a square matrix recording the pairwise Jacaard distance
//between the confidence set of two reference clusters, space allocated assumed
//Only output: confpts, npts, avetight, avecov, avejacaard, rinclude, csetdist
void AveConfset(CLink **clist2, int numcls, int *nsamples, float alpha, int **confpts, int *npts, float *avetight, float *avecov, float *avejacaard, float *rinclude, float *csetdist)
{
    int i,j;
    unsigned char *keepcls;
    float *cvp, v1;

    for (i=0;i<numcls;i++) {
        if (nsamples[i]>1) {
            confset(clist2[i], nsamples[i], alpha, confpts+i, npts+i, &keepcls, &cvp, avetight+i, avecov+i);
            for (j=0, rinclude[i]=0.0;j<nsamples[i];j++)
                if (keepcls[j]) rinclude[i]+=1.0;
            rinclude[i]/=(float)nsamples[i];

            v1=0.0;
            for (j=1;j<nsamples[i];j++)
                v1+=Jacaard_pts(clist2[i][0].id, clist2[i][0].n, clist2[i][j].id, clist2[i][j].n);
            avejacaard[i]=v1/(float)(nsamples[i]-1);

            free(keepcls);
            free(cvp);
        }
        else {
            //artificial setup, confidence set is the reference set
            //which is not matched with any
            avetight[i]=avecov[i]=1.0;
            rinclude[i]=1.0;
        }
    }

    for (i=0;i<numcls;i++) {
        csetdist[i*numcls+i]=0.0;
        for (j=i+1;j<numcls;j++){
            csetdist[i*numcls+j]=Jacaard_pts(confpts[i], npts[i], confpts[j],npts[j]);
            csetdist[j*numcls+i]=csetdist[i*numcls+j];
        }
    }

}

// [[Rcpp::export]]
List ACPS(IntegerVector x, int nbs, int reference)
{
    int i,j,k,n,m=x.size();
    int len;
    int *cls;
    float *wt, *res, *dist;
    int *numcls;
    float thred=0.8;
    int equalcls=0;
    int usesplit=0;
    float alpha=0.1;

    if (m%nbs!=0) {
        Rcpp::stop("Wrong input value");
    } else {
        len=m/nbs;
    }

    cls=(int *)calloc(len*nbs,sizeof(int));

    for (m=0,n=len*nbs;m<n;m++) {
        cls[m]=x[m];
    }

    if (reference!=0 && reference!=1) {
        Rcpp::stop("Wrong reference value");
    }
    else if (reference == 0){
        /*------------------------------------------------------------*/
        /*- Iteration for nbs times                                  -*/
        /*------------------------------------------------------------*/
    
        float *avedist;
        int iter;

        avedist=(float *)calloc(nbs,sizeof(float));
        for (iter=0;iter<nbs;iter++) {
            // Insert the bootstrap sample as the current reference
            for(m=0;m<len;m++) {
                cls[m] = cls[len*iter+m];
            }

        /*------------------------------------------------------------*/
        /*- Alignment done here to get the key matching matrix wt[]  -*/
        /*------------------------------------------------------------*/
    
            align(cls,nbs,len,&wt,&numcls,&dist,equalcls);
  
        /*------------------------------------------------------------*/
        /*- Output the weight matrix and the distances between        */
        /*- clustering results.                                       */
        /*------------------------------------------------------------*/
            for (i=0,avedist[iter]=0.0;i<nbs;i++){
                avedist[iter]+=dist[i];
            }
            avedist[iter]/=(float)nbs; //v1 is the average Distance to the other bootstrap samples

            free(wt);
            free(numcls);
            free(dist);
        }

        NumericVector y3(nbs);
        for (i=0;i<nbs;i++){
            y3[i]=avedist[i];
        }

        List L=List::create(Named("avedist")=y3);
        return L;
    }
    /*------------------------------------------------------------*/
    /*- Alignment done here to get the key matching matrix wt[]  -*/
    /*------------------------------------------------------------*/
    
    align(cls,nbs,len,&wt,&numcls,&dist,equalcls);
    
    /*------------------------------------------------------------*/
    /*- Analysis done based on the matching weight matrix.        */
    /*------------------------------------------------------------*/
    int n0=numcls[0], nr;
    int **codect, *clsct0;
    float **nfave;
    
    codect=(int **)calloc(n0,sizeof(int *));
    nfave=(float **)calloc(n0,sizeof(float *));
    for (i=0;i<n0;i++) {
        codect[i]=(int *)calloc(4,sizeof(int));
        nfave[i]=(float *)calloc(4,sizeof(float));
    }
    
    for (i=1,nr=0;i<nbs;i++){ nr+=numcls[i];}
    res=(float *)calloc(nr*n0,sizeof(float));
    
    //Convert weight matrix wt[] into match-split status matrix res[]
    //Summarize the reliability of all the clusters
    MatchSplit(wt, res, numcls, nbs, codect, nfave, thred);
    
    //Compute confidence sets and average ratios for matched/split clusters
    CLink **clist2;
    int *nsamples;
    int **confpts, *npts, *matchID;
    float *avetight, *avecov, *avejacaard, *rinclude, *csetdist;
    
    matchID=(int *)calloc(len*nbs,sizeof(int));
    confpts=(int **)calloc(numcls[0],sizeof(int *));
    npts=(int *)calloc(numcls[0],sizeof(int));
    avetight=(float *)calloc(numcls[0],sizeof(float));
    avecov=(float *)calloc(numcls[0],sizeof(float));
    avejacaard=(float *)calloc(numcls[0],sizeof(float));
    rinclude=(float *)calloc(numcls[0],sizeof(float));
    csetdist=(float *)calloc(numcls[0]*numcls[0],sizeof(float));
    
    MatchCluster(res,numcls,nbs,thred,cls,len,&clist2, &nsamples, usesplit, matchID);
    AveConfset(clist2, numcls[0], nsamples, alpha, confpts, npts, avetight, avecov, avejacaard, rinclude, csetdist);
    
    /*------------------------------------------------------------*/
    /*- Output results                                            */
    /*------------------------------------------------------------*/
    //Output per cluster summary
    clsct0=(int *)calloc(n0,sizeof(int));
    for (i=0;i<len;i++) {
        if (cls[i]>=0)  {clsct0[cls[i]]++;}
    }
    /*------------------------------------------------------------*/
    /*-                      RETURN                              -*/
    /*------------------------------------------------------------*/
    IntegerMatrix id(len,nbs);
    for(i=0;i<len;i++){
        for(j=0;j<nbs;j++){
            id(i,j)=matchID[i+j*len]+1;
        }
    }
    
    NumericVector y1(nbs);
    for (i=0;i<nbs;i++){
        y1[i]=dist[i];
    }
    
    NumericVector y2(nbs);
    for (i=0;i<nbs;i++){
        y2[i]=numcls[i];
    }
    
    NumericMatrix st(n0,6);
    for(i=0;i<n0;i++){
        st(i,0)=i;
        st(i,1)=clsct0[i];
        st(i,2)=rinclude[i];
        st(i,3)=avetight[i];
        st(i,4)=avecov[i];
        st(i,5)=1-avejacaard[i];
    }
    colnames(st) = CharacterVector::create("cls_id","num_obs","rinclude","avg_tight","avg_cov","avg_jaccard");
    
    
    IntegerMatrix ma(n0,4);
    for(i=0;i<n0;i++){
        for(j=0;j<4;j++){
            ma(i,j)=codect[i][j];
        }
    }
    colnames(ma) = CharacterVector::create("match","split","merge","l.c.");
    
    
    NumericMatrix cap(n0,n0);
    for (i=0;i<n0;i++){
        for (j=0;j<n0;j++)
            cap(i,j)=csetdist[i*n0+j];
    }
    
    //Output confidence set for each reference cluster
    //Every cluster takes one row. Cluster ID, Cluster size, points IDs one by one
    
    IntegerMatrix cps(n0,len);
    for (i=0;i<n0;i++) {
        for (j=0;j<npts[i];j++){
            cps(i,confpts[i][j])=1;
        }
    }
    
    NumericMatrix wtm(nr,n0);
    for (i=1,m=0;i<nbs;i++){
        for (j=0;j<numcls[i];j++){
            for (k=0;k<numcls[0];k++){
                wtm(m,k) = wt[m*numcls[0]+k];
            }
            m++;
        }
    }
    
    NumericMatrix resm(nr,n0);
    for (i=1,m=0;i<nbs;i++){
        for (j=0;j<numcls[i];j++){
            for (k=0;k<numcls[0];k++){
                resm(m,k) = res[m*numcls[0]+k];
            }
            m++;
        }
    }

    
    List L=List::create(Named("distance")=y1,
                        Named("numcls")=y2,
                        Named("statistics")=st,
                        Named("cap")=cap,
                        Named("id")=id,
                        Named("cps")=cps,
                        Named("match")=ma,
                        Named("weight")=wtm,
                        Named("topo_result")=resm);
    
    free(cls);
    free(wt);
    free(res);
    free(dist);
    free(numcls);
    free(codect);
    free(clsct0);
    free(nfave);
    free(nsamples);
    free(confpts);
    free(npts);
    free(matchID);
    free(avetight);
    free(avecov);
    free(avejacaard);
    free(rinclude);
    free(csetdist);
    
    return L;
}






