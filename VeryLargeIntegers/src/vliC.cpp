
#include <Rcpp.h>
using namespace Rcpp;

int maxC(IntegerVector x);
int compareC(List x, List y);
int onebitsC(List x);
int negbinC(int r, double p);
int randC(int n);
String zerossC(String x);
String replaceC(std::string& s, const std::string& from, const std::string& to);
String collapseC(IntegerVector x);
String addzerosC(String x, int size);
String addrzerosC(String x, int n);
String printvliC(List x);
IntegerVector repC(int x, int n);
IntegerVector splitC(String x, int n);
IntegerVector appendC(IntegerVector a, IntegerVector b);
IntegerVector zerosvC(IntegerVector x);
IntegerVector selC(IntegerVector a, int n, int m);
IntegerVector divp2C(IntegerVector x, int k);
List equalC(IntegerVector x, IntegerVector y);
List vliC(String x);
List vlivC(int sign, IntegerVector x);
List sumC(List x, List y);
List subC(List x, List y);
List mulbaseC(List x, List y);
List mulC(List x, List y);
List karC(List x, List y);
List fkzerosC(List x, int k);
List divbaseC(List x, List y);
List expC(List x, List n);
List rootC(List x);
List rootkC(List x);
List factbaseC(List x);
List binomC(List n, List k);
List randomvliC(int n);
bool gtC(List x, List y);
bool ltC(List x, List y);
bool geqC(List x, List y);
bool leqC(List x, List y);
bool eqC(List x, List y);
bool neqC(List x, List y);

List zero = vlivC(1, IntegerVector::create(0));
List one = vlivC(1, IntegerVector::create(1));
List two = vlivC(1, IntegerVector::create(2));


// [[Rcpp::export]]
String zerossC(String x){
  std::string a = x;
  while (a.substr(0,1) == "0"){
    a = a.substr(1);
  }
  if (a == "") a = "0";
  return a;
}

// [[Rcpp::export]]
IntegerVector zerosvC(IntegerVector x){
  IntegerVector a = x;
  int n = a.size();
  while ((n>1) & (a[0] == 0)){
    IntegerVector b(n-1);
    for (int i = 0; i < (n-1); i = i + 1){
      b[i] = a[i+1];
    }
    a = b;
    n = a.size();
  }
  return a;
}

// [[Rcpp::export]]
IntegerVector repC(int x, int n) {
  IntegerVector out(n);
  for (int i = 0; i < n; i = i + 1) {
    out[i] = x;
  }
  return out;
}

// [[Rcpp::export]]
String replaceC(std::string& s, const std::string& from, const std::string& to){
  if(!from.empty())
    for(size_t pos = 0; (pos = s.find(from, pos)) != std::string::npos; pos += to.size())
      s.replace(pos, from.size(), to);
  return s;
}

// [[Rcpp::export]]
int maxC(IntegerVector x){
  int t = x[0];
  for (int i = 1; i < x.size(); i = i+1){if (x[i] > t){t=x[i];}}
  return t;
}

// [[Rcpp::export]]
IntegerVector appendC(IntegerVector a, IntegerVector b){
  int sa = a.size();
  int sb = b.size();
  IntegerVector out(sa+sb);
  for (int i = 0; i < sa; i = i+1){
    out[i] = a[i];
  }
  for (int i = 0; i < sb; i = i+1){
    out[i+sa] = b[i];
  }
  return(out);
}

// [[Rcpp::export]]
IntegerVector selC(IntegerVector a, int n, int m){
  int la = a.size();
  if (n<0){n=0;}
  else if(n>la){n=la;}
  if (m<0){m=0;}
  else if(m>la){m = la;}
  if ((n == m) & (n == 0)){
    return(IntegerVector::create(0));
  }
  else if (n>=m){m=n+1;}
  int lo = m-n;
  IntegerVector out(lo);
  for (int i = 0; i < lo; i = i+1){
    out[i] = a[i+n];
  }
  return(out);
}

// [[Rcpp::export]]
String addzerosC(String x, int size){
  std::string a;
  std::string b = x;
  unsigned s = size;
  int l = b.length();
  while ((a.length()+l) < s){
    a = a + "0";
  }
  return(a+b);
}

// [[Rcpp::export]]
String collapseC(IntegerVector x){
  std::string a;
  for (int i=0; i < x.size(); i = i+1){
    String b = x[i];
    b = addzerosC(b,4);
    std::string c = b;
    a = a + c;
    a = zerossC(a);
  }
  return(a);
}

// [[Rcpp::export]]
String addrzerosC(String x, int n){
  std::string a = x;
  for (int i = 0; i < n; i = i + 1){
    a = a + "0";
  }
  return(a);
}

// [[Rcpp::export]]
List equalC(IntegerVector x, IntegerVector y){
  int lx, ly, m;
  IntegerVector z;
  List out(2);
  lx = x.size();
  ly = y.size();
  if (lx<ly){
    m = ly-lx;
    z = appendC(repC(0,m), x);
    x = z;
  }
  else if (lx>ly){
    m = lx-ly;
    z = appendC(repC(0,m), y);
    y = z;
  }
  out[0] = x;
  out[1] = y;
  return out;
}

// [[Rcpp::export]]
List fkzerosC(List x, int k){
  int n = x[1];
  IntegerVector t = x[2];
  int s = t.size();
  IntegerVector a(s+k);
  for (int i = 0; i < s; i = i+1){a[i] = t[i];}
  for (int i = s; i < (s+k); i = i+1){a[i] = 0;}
  List out(3);
  out[0] = x[0];
  out[1] = n + (4*k);
  out[2] = a;
  out.names() = CharacterVector::create("sign", "length", "value");
  out.attr("class") = "vli";
  return(out);
}

// [[Rcpp::export]]
IntegerVector splitC(String x, int n){
  std::string a = x;
  if (n<4){
    IntegerVector out(1);
    out = atoi(a.c_str());
    return out;
  }
  else{
    int min_ = 1;
    if (n%4 == 0){
      min_ = 0;
    }
    int parts = trunc(n/4) + min_;
    IntegerVector out(parts);
    IntegerVector begin(parts);
    begin[0] = 0;
    int from = n-((parts-1)*4); int to = n-4;
    int k = 1;
    for (int i = from; i <= to; i = i + 4) {
      begin[k] = i;
      k += 1;
    }
    out[0] = atoi(a.substr(begin[0], n-((parts-1)*4)).c_str());
    for (int l = 1; l < parts; l = l + 1){
      out[l] =  atoi(a.substr(begin[l], 4).c_str());
    }
    return out;
  }
}

// [[Rcpp::export]]
List vliC(String x){
  List out(3);
  std::string a = x;
  a = replaceC(a,",","");
  a = replaceC(a,".","");
  String sign;
  out[0] = 1;
  if (a.substr(0,1) == "-"){
    out[0] = -1;
    a = a.substr(1);
  }
  else if (a.substr(0,1) == "+"){
    a = a.substr(1);
  }
  a = zerossC(a);
  if(a=="0"){out[0]=1;}
  int n = a.length();
  IntegerVector t = splitC(a,n);
  out[1] = n;
  out[2] = t;
  out.names() = CharacterVector::create("sign", "length", "value");
  out.attr("class") = "vli";
  return out;
}

// [[Rcpp::export]]
List vlivC(int sign, IntegerVector x){
  List out(3);
  IntegerVector a = zerosvC(x);
  String f = a[0];
  std::string m = f;
  int n = (a.size()*4) - (4-m.size());
  out[0] = sign;
  out[1] = n;
  out[2] = a;
  out.names() = CharacterVector::create("sign", "length", "value");
  out.attr("class") = "vli";
  return out;
}

// [[Rcpp::export]]
String printvliC(List x){
  std::string sign = "";
  IntegerVector v = x[2];
  int signx = x[0];
  std::string value = collapseC(x[2]);
  IntegerVector z = 0;
  if (signx == -1){sign = "-";}
  return (sign + value);
}

// [[Rcpp::export]]
List sumC(List x, List y){
  int xsign = x[0];
  int ysign = y[0];
  bool fgt = TRUE;
  if (gtC(y,x)){ fgt = FALSE;}
  if (xsign != ysign){
    if (xsign == 1){
      if (fgt){return(subC(x, vlivC(1 ,y[2])));}
      else {return(subC(y, vlivC(-1 ,x[2])));}
    }
    else{
      if (fgt){return(subC(x, vlivC(-1 ,y[2])));}
      else {return(subC(y, vlivC(1 ,x[2])));}
    }
  }
  else{
    IntegerVector a = x[2];
    IntegerVector b = y[2];
    List ab = equalC(a,b);
    a = ab[0];
    b = ab[1];
    int m = a.size();
    IntegerVector sum(m);
    int c = 0;
    for (int i = m-1; i>=0; i = i-1){
      sum[i] = a[i] + b[i] + c;
      if (sum[i] >= 10000){
        sum[i] = sum[i] - 10000;
        c = 1;
      }
      else c = 0;
    }
    if (c == 1){sum = appendC(IntegerVector::create(1),sum);}
    return(vlivC(xsign, sum));
  }
}

// [[Rcpp::export]]
List subC(List x, List y){
  int xsign = x[0];
  int ysign = y[0];
  if (xsign != ysign){
    if (xsign == 1){return(sumC(x, vlivC(1 ,y[2])));}
    else {return(vlivC(-1, sumC(vlivC(1,x[2]), y)[2]));}
  }
  else{
    IntegerVector a = x[2];
    IntegerVector b = y[2];
    List out(3);
    int sign = xsign;
    if (ltC(vlivC(1,a),vlivC(1,b))){
      sign = -1 * xsign;
      IntegerVector c=a;
      a=b;
      b=c;
    }
    int la = a.size();
    int lb = b.size();
    int dif = la-lb;
    if (dif>0){
      IntegerVector b_(la);
      for (int i=0; i<dif; i++){
        b_[i] = 0;
      }
      for (int i=dif; i<la; i++){
        b_[i] = b[i-dif];
      }
      b = b_;
    }
    IntegerVector sub(la);
    int carry = 0;
    for (int i = la-1; i >= 0; i = i-1){
      if (a[i]<b[i]){sub[i] = a[i]+(10000) - (b[i]+carry); carry=1;}
      else {sub[i]= a[i] - (b[i]+carry); carry=0;}
    }
    return(vlivC(sign, sub));
  }
}

// [[Rcpp::export]]
int compareC(List x, List y){
  int sx = x[0];
  int sy = y[0];
  IntegerVector x2 = x[2];
  IntegerVector y2 = y[2];
  List ab = equalC(x2, y2);
  IntegerVector a = ab[0];
  IntegerVector b = ab[1];
  int l = a.size();
  int t = 0;
  int i=0;
  while ((t == 0) & (i<l)){
    if (sx*a[i] > sy*b[i]) {
      t = 1;
    }
    else if (sx*a[i] < sy*b[i]){
      t = 2;
    }
    i = i+1;
  }
  a[0] = abs(a[0]);
  b[0] = abs(b[0]);
  return(t);
}

// [[Rcpp::export]]
bool gtC(List x, List y){
  int t = compareC(x,y);
  return(t == 1);
}

// [[Rcpp::export]]
bool ltC(List x, List y){
  int t= compareC(x,y);
  return(t == 2);
}

// [[Rcpp::export]]
bool geqC(List x, List y){
  int t = compareC(x,y);
  return((t == 0) | (t == 1));
}

// [[Rcpp::export]]
bool leqC(List x, List y){
  int t = compareC(x,y);
  return((t == 0) | (t == 2));
}

// [[Rcpp::export]]
bool eqC(List x, List y){
  int t = compareC(x,y);
  return(t == 0);
}

// [[Rcpp::export]]
bool neqC(List x, List y){
  int t = compareC(x,y);
  return((t == 1) | (t == 2));
}

// [[Rcpp::export]]
List mulbaseC(List x, List y){
  int a = x[0]; int b = y[0];
  IntegerVector l = IntegerVector::create(a,b);
  if ( maxC(l) <= 40 ){return(mulC(x,y));}
  else return(karC(x,y));
}

// [[Rcpp::export]]
List mulC(List x, List y){
  int signx = x[0];
  int signy = y[0];
  IntegerVector a = x[2];
  IntegerVector b = y[2];
  IntegerVector len = IntegerVector::create(a.size(), b.size());
  IntegerVector parc(len[0]+1);
  List prod(len[1]);
  for (int i = 0; i < len[1]; i = i + 1){
    int carry = 0;
    for (int j = len[0]; j > 0; j = j - 1){
      int p = a[j-1]*b[len[1]-(i+1)];
      parc[j] = ((p + carry)%(10000));
      carry = trunc((p + carry)/(10000));
    }
    parc[0] = carry;
    prod[i] = addrzerosC(collapseC(parc),i*4);
  }
  List out = zero;
  for(int i = 0; i < len[1]; i = i + 1){
    out = sumC(out, vliC(prod[i]));
  }
  out[0] = signx * signy;
  return(out);
}

 // [[Rcpp::export]]
List karC(List x, List y){
  IntegerVector a = x[2];
  IntegerVector b = y[2];
  int signx = x[0];
  int signy = y[0];
  IntegerVector l = IntegerVector::create(a.size(),b.size());
  int sign = signx * signy;
  if (maxC(l) < (10)){
      List out = mulC(x,y);
    out[0] = sign;
    return(out);
  }
  else{
    int size = maxC(IntegerVector::create(
    trunc(l[0]/2) + l[0]%2,
    trunc(l[1]/2) + l[1]%2));
    List a1, a2, b1, b2;
    if (l[0] > 1){
      a1 = vlivC(1,selC(a, 0, l[0]-size));
      a2 = vlivC(1,selC(a, l[0]-size, l[0]));
    }
    else{
      a1 = zero;
      a2 = vlivC(1,a);
    }
    if (l[1] > 1){
      b1 = vlivC(1,selC(b, 0, l[1]-size));
      b2 = vlivC(1,selC(b, l[1]-size, l[1]));
    }
    else{
      b1 = zero;
      b2 = vlivC(1,b);
    }
    List z2 = karC(a1,b1);
    List z0 = karC(a2,b2);
    List z1 = subC(subC(karC(sumC(a1,a2),sumC(b1,b2)),z2) , z0);
    List r =  sumC(sumC(fkzerosC(z2,2*size), fkzerosC(z1,size)) , z0);
    r[0] = sign;
    return(r);
  }
}

// [[Rcpp::export]]
IntegerVector divp2C(IntegerVector x, int k){
  IntegerVector v = x;
  for (int n = 0; n < k; n++){
    IntegerVector out(v.size());
    int c = 0;
    for (int i = 0; i < v.size(); i++){
      out[i] = trunc(v[i]/2) + (c*5000);
      c = v[i]%2;
    }
    v = out;
  }
  return(zerosvC(v));
}

// [[Rcpp::export]]
List divbaseC(List x, List y){
  List out(2); List qu(3); List r(3);
  r = x;
  List b = y;
  List bcopy = y;
  int sx = r[0];
  int sy = b[0];
  /* choosing output sign */
  int s = 1;
  if (sx != sy){s = -1;}
  r[0] = 1;
  b[0] = 1;
  IntegerVector r2 = r[2];
  IntegerVector b2 = b[2];
  /* trivial case */
  if ((r2.size() == 1) & (b2.size() == 1)){
    int r20 = r2[0];
    int b20 = b2[0];
    int q = trunc(r20/ b20);
    if(q == 0){qu = zero;}
    else {qu = vlivC(s, IntegerVector::create(q));}
    r = vlivC(1, IntegerVector::create(r20 % b20));
    /* negative dividend subcase  */
    if (sx == -1){
      r = subC(vlivC(1, bcopy[2]), r);
      qu = sumC(qu, one);
    }
  }
  /* greater divisor case */
  else if(ltC(x,y)){
    qu = zero;
    /* negative dividend subcase  */
    if (sx == -1){
      r = subC(vlivC(1, bcopy[2]), r);
      qu = sumC(qu, one);
    }
    if (neqC(qu, zero)){qu[0] = s;}
  }
  /* main case  */
  else{
    /* normalizing divisor */
    int b20 = b2[0];
    int k=0;
    int mul = 1;
    while(b20 < 5000){
      k = k+1;
      mul = mul*2;
      b = mulbaseC(two ,b);
      b2 = b[2];
      b20 = b2[0];
    }
    List vlimul = vlivC(1, IntegerVector::create(mul));
    r =  mulbaseC(vlimul, r);
    r2 = r[2];
    /* division algorithm */
    int n = b2.size(); int m = r2.size() - n;
    IntegerVector q(m+1);
    List t = fkzerosC(b,m);
    if (geqC(r,t)){q[0]=1;r=(subC(r,t));}
    for (int j=0; j<m; j++){
      r2 = r[2];
      q[j+1] = trunc(((r2[0]*10000)+(r2[1]))/b2[0]);
      if (9999 < q[j+1]){q[j+1]=9999;}
      r = subC(r, fkzerosC((mulbaseC(vliC(collapseC(IntegerVector::create(q[j+1]))),b)), m-(j+1)));
      while (ltC(r,zero)){
        q[j+1] = q[j+1] - 1;
        r = sumC(r, fkzerosC(b,m-(j+1)));
      }
    }
    qu = vliC(collapseC(q));
    r = vlivC(1, divp2C(r[2],k));
    /* negative dividend subcase  */
    if (sx == -1){
      r = subC(vlivC(1, bcopy[2]), r);
      qu = sumC(qu, one);
    }
    if (neqC(qu, zero)){qu[0] = s;}
  }
  /* restoring original signs */
  x[0] = sx;
  y[0] = sy;
  /* exporting */
  out[0] = qu;
  out[1] = r;
  return(out);
}

// [[Rcpp::export]]
List expC(List x, List n){
  List a = x;
  List k = n;
  int sa = a[0];
  int sk = k[0];
  if (sk == -1){ /* negative exponent case */
    return(one);
  }
  else{ /* positive exponent case */
    List r = (divbaseC(k, two))[1];
    IntegerVector r2 = r[2];
    int rest = r2[0];
    List out = one;
    if (rest == 1){
      List exp = divbaseC((subC(k, one)), two)[0];
      List base = mulbaseC(a,a);
      for (List i=one; leqC(i, exp); i = sumC(i,one)){
        out = mulbaseC(out, base);
      }
      out = mulbaseC(a,out);
    }
    else{
      List exp = divbaseC(k, two)[0];
      List base = mulbaseC(a,a);
      for (List i=one; leqC(i, exp); i = sumC(i,one)){
        out = mulbaseC(out,base);
      }
    }
    /* sign */
    IntegerVector k2 = k[2];
    if ((sa == -1) & (k2[k2.size()-1] % 2 == 1)){
      out[0] = -1;
    }
    return(out);
  }
}

// [[Rcpp::export]]
List rootC(List x){
  List a = x;
  List out(2);
  if (ltC(a,two)){
    out[0] = a;
    out[1] = zero;
    return(out);
  }
  else{
    List u = divbaseC(a, two)[0];
    bool rep = TRUE;
    List s(3);
    List t(3);
    while (rep){
      s = u;
      t = sumC(s, divbaseC(a,s)[0]);
      u = divbaseC(t, two)[0];
      if (geqC(u, s)){rep = FALSE;}
    }
    out[0] = s;
    out[1] = subC(a, (mulbaseC(s,s)));
    return(out);
  }
}

// [[Rcpp::export]]
List rootkC(List x, List k){
  List a = x;
  List n = k;
  List n_1 = subC(n,one);
  List out(2);
  if (ltC(a,two)){
    out[0] = a;
    out[1] = zero;
    return(out);
  }
  else{
    List u = divbaseC(a, two)[0];
    bool rep = TRUE;
    List s(3);
    List t(3);
    while (rep){
      s = u;
      t = sumC( (mulbaseC(n_1, s)), divbaseC(a, expC(s, n_1))[0]);
      u = divbaseC(t, n)[0];
      if (geqC(u, s)){rep = FALSE;}
    }
    out[0] = s;
    out[1] = subC(a, (expC(s, n)));
    return(out);
  }
}

// [[Rcpp::export]]
int onebitsC(List x){
  List a = x;
  int ones = 0;
  List ar;
  List r;
  while (gtC(a, zero)){
    ar = divbaseC(a, two);
    a = ar[0];
    r = ar[1];
    if(eqC(r,one)){ones++;}
  }
  return(ones);
}

// [[Rcpp::export]]
List factbaseC(List x){
  List a = x;
  int la = a[1];
  List out = one;
  if (la > 40){
    for (List i = one; leqC(i,a); i = sumC(i,one)){
      out = karC(i, out);
    }
  }
  else{
    for (List i = one; leqC(i,a); i = sumC(i,one)){
      out = mulbaseC(i, out);
    }
  }
  return(out);
}

int randC(int n){
  NumericVector r = runif(1);
  int p = 1;
  for (int i = 0; i < n; i++){p *= 10;}
  int out = r[0] * p;
  return(out);
}

// [[Rcpp::export]]
List randomvliC(int n){
  int l = trunc(n/4);
  int n4 = n%4;
  if (n4 > 0) l = l + 1;
  IntegerVector t(l);
  int b = 0;
  if (n4 > 0){
    t[0] = randC(n4);
    b=1;
  }
  for (int i = b; i < l; i++){
    t[i] = randC(4);
  }
  return(vlivC(1, t));
}

// [[Rcpp::export]]
int negbinC(int r, double p){
  int succ = 0;
  int k = 0;
  while (succ < r){
    double ra = runif(1)[0];
    if ( ra < p ){
      succ = succ + 1;
    }
    k = k + 1;
  }
  return(k);
}

// [[Rcpp::export]]
List binomC(List n, List k){
  List out = one;
  for( List i = zero; ltC(i, k); i = sumC(i, one) ){
    out = mulbaseC(out, subC(n, i));
    out = divbaseC(out, sumC(i, one))[0];
  }
  return(out);
}


