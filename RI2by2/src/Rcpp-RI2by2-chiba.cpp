#include <Rcpp.h>
#include <algorithm>
#include <map>
using namespace Rcpp;

// // Modular function to construct ECDF (for Blaker method of CIs) ===========
double ECDF_pv(std::vector<double> teststat, std::vector<double> prob,
               double teststat_obs) {

  std::map<double, double> m;
  int n=teststat.size();
  for (int i=0; i<n; ++i) {
    m[teststat[i]] += prob[i];
  }
  std::vector<double> out_teststat;
  std::vector<double> out_prob;
  std::vector<double> out_cumprob;
  std::vector<double> out_cumprob_rev;
  double cumprob=0.0;
  double critical_value=0.0;
  std::map<double, double>::iterator it = m.begin();
  while(it != m.end()) {
    out_teststat.push_back(it->first);
    out_prob.push_back(it->second);
    cumprob += it->second;
    out_cumprob.push_back(cumprob);
    double cumprob_rev=1.0-cumprob+(it->second);
    out_cumprob_rev.push_back(cumprob_rev);
    if (std::abs((it->first)-teststat_obs)<
      (std::numeric_limits<double>::epsilon()*1e1)) {
      critical_value = std::max(critical_value,
                                std::min(cumprob, cumprob_rev));
    }
    it++;
  }
  double pvalue=0.0;
  for (int j=0; j<out_teststat.size(); ++j) {
    if( std::min(out_cumprob[j], out_cumprob_rev[j]) <=
        (critical_value+std::numeric_limits<double>::epsilon()*1e1) ) {
      pvalue += out_prob[j];
    }
  }
  return( pvalue );
}


// [[Rcpp::export(name = ".CI_2by2_chiba_tau_v7_approx")]]
List CI_2by2_chiba_tau_v7_approx(int a, int b, int c, int d, double level,
                                 int total_tests) {
  int n = a+b+c+d;
  int n11i = a+c;
  int n10i = a+d;
  int n01i = b+c;
  int n00i = b+d;
  std::vector<double> RDtrue;
  std::vector<double> p_l;
  std::vector<double> p_u;
  std::vector<double> p_2;
  std::vector<double> p_21;

  double RDo = a/((double) a+b)-c/((double) c+d);

  int nCk_dim = n + 1;
  NumericMatrix nCk( nCk_dim, nCk_dim );
  std::fill( nCk.begin(), nCk.end(), -1.0 );
  double normconst = ::Rf_lchoose(n,a+b);

  List nvector_list_out; // Checking order of n vectors; remove in production!
  List pvector_list_out; // Checking order of n vectors; remove in production!

  // Vector of possible values for r-s
  std::vector<int> r_minus_s_values;
  // Map with vector of ordered (q,r,s) values for each r-s
  std::map<int, std::vector< std::vector<int> > > map_qrs;

  for (int r_minus_s=-n01i; r_minus_s<=n10i; ++r_minus_s) {
    r_minus_s_values.push_back(r_minus_s);
    std::vector< std::vector<int> > qrsvalues_inc;
    // range of values for s
    int s_min_tau0 = std::max(0,-r_minus_s);
    int s_max_tau0 = std::min(n01i,n10i-r_minus_s);
    for (int s = s_min_tau0; s<=s_max_tau0; ++s) {
      // in v6, (r,s) sorted jointly by r-s, so searching from upper limit
      // requires reordering (r,s) s.t. small values of s are tested first.
      // Here, the (r,s) have the same order for each fixed value of r-s,
      // since r-s is ordered separately from (r,s).
      int r = r_minus_s+s; // fixed (r, s)
      // linear space for a given (r,s)
      int q_max_rs = std::min(n11i, std::min(n-b-r,n-d-s));
      int q_min_rs = std::max(0, std::max(c-s,a-r));
      q_max_rs = std::max(q_min_rs,q_max_rs);
      for (int q=q_max_rs; q>=q_min_rs; --q) {
        // // start with larger values of q
        int t = n-(q+r+s);
        if (t>=0 && t<=n00i) {
          std::vector<int> qrs(3, q);
          qrs[1] = r;
          qrs[2] = s;
          qrsvalues_inc.push_back(qrs);
        }
      }
    }
    // save the vector of vectors for fixed tau0
    map_qrs[r_minus_s] = qrsvalues_inc;
  }
  bool find_lower_bound=true;

  int rs_i=0;
  int rs_size=r_minus_s_values.size();
  // number of tests for each value of tau: at least 2 points (first,last)
  int apx = std::max(2, total_tests/rs_size);
  do {
    std::vector<double> pv_max_rs(4, 0.0);
    int r_minus_s=r_minus_s_values[rs_i];
    std::vector< std::vector<int> > qrsvalues = map_qrs[r_minus_s];
    double RDrs = ((double) r_minus_s)/((double) n);

    // // Test only selected values of (q,r,s) w.r.t r-s i.e. approx!! ===
    // // // (max) p-value for given (r,s) \leq exact p-value
    int q_num_rs = qrsvalues.size();
    std::vector<int> qrs_indices;
    if (apx >= q_num_rs) {
      for (int qrs_i=0; qrs_i<q_num_rs; ++qrs_i) {
        qrs_indices.push_back(qrs_i);
      }
    } else {
      // // evenly-spaced points
      for (int qrs_i=0; qrs_i<apx; ++qrs_i) {
        qrs_indices.push_back(qrs_i*(q_num_rs-1)/(apx-1));
      }
    }
    // // ================================================================

    for (int qrs_i=0; qrs_i<qrs_indices.size(); ++qrs_i) {
      std::vector<int>  qrs = qrsvalues[qrs_indices[qrs_i]];
      int q = qrs[0];
      int r = qrs[1];
      int s = qrs[2];
      int t = n-(q+r+s);

      // Pv_2by2_chiba function with caching for lchoose ==============
      // std::vector<double> pvs = Pv_2by2_chiba(q,r,s,t,a,b,c,d);
      // Additional Blaker method from Agresti (2003)
      int n11=q; int n10=r; int n01=s; int n00=t;
      std::vector<int> nvector;
      nvector.push_back(n11);
      nvector.push_back(n10);
      nvector.push_back(n01);
      nvector.push_back(n00);

      double tau0 = RDrs;
      std::vector<double> p_chibaL;
      std::vector<double> p_chibaU;
      std::vector<double> p_chiba2;
      std::vector<double> p_all;
      std::vector<double> RDn_all;
      int m=a+b;
      for (int i=0; i<=n11; ++i) {
        for (int j=0; j<=n10; ++j) {
          for (int k=0; k<=n01; ++k) {
            int l = m-(i+j+k);
            if (l>=0 && l<=n00) {
              double p_tmp = 0.0;
              double RDn = (i+j)/((double) m)-
                ((n11-i)+(n01-k))/((double) n-m);

              bool Iz_chibaL=(RDn-RDo)>=
                (-std::numeric_limits<double>::epsilon()*1e1);
              bool Iz_chibaU=(RDn-RDo)<=
                (std::numeric_limits<double>::epsilon()*1e1);
              bool Iz_chiba2=std::abs(RDn-tau0)-std::abs(RDo-tau0)>=
                (-std::numeric_limits<double>::epsilon()*1e1);

              std::vector<int> z1;
              z1.push_back(i);
              z1.push_back(j);
              z1.push_back(k);
              z1.push_back(l);

              double p_tmp_z1 = 0.0;
              for (int ntype=0; ntype < 4; ++ntype) {
                int nn = nvector[ntype], kk = z1[ntype];
                // check cache
                if( nCk(nn, kk) < 0.0 ) {
                  // calculate the lchoose value
                  nCk(nn, kk) = ::Rf_lchoose( nn, kk );
                  nCk(nn, nn - kk) = nCk(nn, kk);
                }
                p_tmp_z1 += nCk(nn, kk);
              }
              p_tmp += exp(p_tmp_z1-normconst);

              if (Iz_chibaL) {
                p_chibaL.push_back(p_tmp);
              }
              if (Iz_chibaU) {
                p_chibaU.push_back(p_tmp);
              }
              if (Iz_chiba2) {
                p_chiba2.push_back(p_tmp);
              }
              // Blaker method
              p_all.push_back(p_tmp);
              RDn_all.push_back(RDn);
            }
          }
        }
      }
      std::vector<double> pvs(4);
      pvs[0] = std::accumulate(p_chibaL.begin(), p_chibaL.end(), 0.0);
      pvs[1] = std::accumulate(p_chibaU.begin(), p_chibaU.end(), 0.0);
      pvs[2] = std::accumulate(p_chiba2.begin(), p_chiba2.end(), 0.0);
      pvs[3] = ECDF_pv(RDn_all,p_all,RDo);
      nvector_list_out.push_back(nvector);
      pvector_list_out.push_back(pvs);
      // ==============================================================
      for (int pf=0; pf<4; pf++) {
        pv_max_rs[pf] = std::max(pv_max_rs[pf], pvs[pf]);
      }
      bool p_1b = std::min(pv_max_rs[0], pv_max_rs[1]) >=
        (level-std::numeric_limits<double>::epsilon()*1e1)/2;
      bool p_2b = pv_max_rs[2] >=
        (level-std::numeric_limits<double>::epsilon()*1e1);
      bool p_2nestb = pv_max_rs[3] >=
        (level-std::numeric_limits<double>::epsilon()*1e1);
      if (p_1b && p_2b && p_2nestb) {
        break;
      }
    }
    p_l.push_back( pv_max_rs[0] );
    p_u.push_back( pv_max_rs[1] );
    p_2.push_back( pv_max_rs[2] );
    p_21.push_back( pv_max_rs[3] );
    RDtrue.push_back( RDrs );

    bool p_2_bound = pv_max_rs[2] >=
      (level-std::numeric_limits<double>::epsilon()*1e1);
    bool p_21_bound = pv_max_rs[3] >=
      (level-std::numeric_limits<double>::epsilon()*1e1);
    bool p_1_bound = false;
    if (!find_lower_bound) {
      p_1_bound = pv_max_rs[1] >=
        (level-std::numeric_limits<double>::epsilon()*1e1)/2;
      if (p_1_bound && p_2_bound && p_21_bound) {
        rs_i = rs_size;
      } else {
        rs_i--;
      }
    } else {
      p_1_bound = pv_max_rs[0] >=
        (level-std::numeric_limits<double>::epsilon()*1e1)/2;
      if (p_1_bound && p_2_bound && p_21_bound) {
        find_lower_bound = false;
        rs_i = rs_size-1;
      } else {
        rs_i++;
      }
    }
  } while(rs_i<rs_size);

  std::vector<double> chiba1(2);
  chiba1[0]=1.0;
  chiba1[1]=-1.0;
  std::vector<double> perm2(2);
  perm2[0]=1.0;
  perm2[1]=-1.0;
  std::vector<double> perm21(2);
  perm21[0]=1.0;
  perm21[1]=-1.0;

  for (int i=0; i<RDtrue.size(); ++i) {
    if (p_l[i] >= (level-std::numeric_limits<double>::epsilon()*1e1)/2) {
      chiba1[0] = std::min(RDtrue[i], chiba1[0]);
    }
    if (p_u[i] >= (level-std::numeric_limits<double>::epsilon()*1e1)/2) {
      chiba1[1] = std::max(RDtrue[i], chiba1[1]);
    }
    if (p_2[i] >= (level-std::numeric_limits<double>::epsilon()*1e1)) {
      perm2[0] = std::min(RDtrue[i], perm2[0]);
      perm2[1] = std::max(RDtrue[i], perm2[1]);
    }
    if (p_21[i] >= (level-std::numeric_limits<double>::epsilon()*1e1)) {
      perm21[0] = std::min(RDtrue[i], perm21[0]);
      perm21[1] = std::max(RDtrue[i], perm21[1]);
    }
  }
  return(List::create(chiba1, perm2, perm21,
                      RDtrue, p_l, p_u, p_2, p_21,
                      nvector_list_out, pvector_list_out));
}
