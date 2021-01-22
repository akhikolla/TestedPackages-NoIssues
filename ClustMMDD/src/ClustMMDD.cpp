/** Begin include */

#include <Rcpp.h>
/** // [[Rcpp::depends(RcppProgress)]]*/
/** #include <progress.hpp> */

#include <Rmath.h>
//#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>

#include "macros.h"

#define PRECISION 9



/** End include */

using namespace Rcpp;

/** ------------------ */
/**
 * @brief	Global variables
 * @note	For R package ClustMMDD
 */
/**
  * @brief  EM parameters.
  * @param  _EPSI Current value for stoping EM;
  * @param  _NBRE_PTS_INITIAUX Current number of small EM;
  * @param  _NBRE_ITER_EM Current number of iterations of small EM;
  * @param  _STOCH  Current indicator if stochastique EM should be performed or not. Default value id false.
  */
const double MIN_EPSI = 1e-20;
const double EPSI = 1e-8;
const double MAX_EPSI = 1e-5;
double _EPSI = EPSI;

const int MIN_TYPE_SMALL_EM = 0;
const int TYPE_SMALL_EM = 0;
const int MAX_TYPE_SMALL_EM = 2;
int _TYPE_SMALL_EM = TYPE_SMALL_EM; // 0 = classic, 1 = SEM, 2 = CEM

const int MIN_TYPE_EM = 0;
const int TYPE_EM = 0;
const int MAX_TYPE_EM = 2;
int _TYPE_EM = TYPE_EM; // 0 = classic, 1 = SEM, 2 = CEM

const int MIN_NBER_SMALL_EM = 10;
const int NBER_SMALL_EM = 20;
const int MAX_NBER_SMALL_EM = 50;
int _NBER_SMALL_EM = NBER_SMALL_EM;

const int MIN_NBER_ITER_EM = 10;
const int NBER_ITER_EM = 15;
const int MAX_NBER_ITER_EM = 50;
int _NBER_ITER_EM = NBER_ITER_EM;

const int MIN_NBER_ITER_LONG_EM = 2000;
const int NBER_ITER_LONG_EM = 5000;
const int MAX_NBER_ITER_LONG_EM = 7000;
int _NBER_ITER_LONG_EM = NBER_ITER_LONG_EM;

Rcpp::CharacterVector _TYPE_EM_NAMES = Rcpp::CharacterVector::create("EM", "CEM", "SEM");

const int NBER_OCCURRENCES_MAX = 6;

/**
  * @brief  Differents available criteria.
  * @param  NBRE_CRITERES Constant given the number of available criteria;
  * @param  CRITERES  List of available criteria;
  * @param  CHOIX_CRIT Current choice of a criterion;
  * @param  CRITERES_WITH_CALIBRATION Number of criteria with calibration.
  */
Rcpp::CharacterVector CRITERIA_NAMES = Rcpp::CharacterVector::create("BIC", "AIC", "ICL", "CteDim");
const int NBER_CRITERIA = 4;

// [[Rcpp::export()]]
int getNberCriteria_Rcpp()
{
    return NBER_CRITERIA;
}

/**
 * @brief	Get the implemented criteria names
 */
// [[Rcpp::export()]]
Rcpp::CharacterVector getCriteriaNames_Rcpp()
{
    return CRITERIA_NAMES;
}

/**
 * @brief	NBER_OCCURRENCES_MAX = ploidy 
 */
//	[[Rcpp::export()]]
int getNberOccurrencesMax()
{
  return NBER_OCCURRENCES_MAX;
}

/**
 * @brief	Choice of the penalized criteria
 */
int _CRITERION_CHOICE;
bool _EXPLORE_WITH_MANY_CRITERIA;
int _EXPLORATION_CRITERION;

/**
    *	@brief	Put threshold on frequancies
    */
bool _PutTHRESHOLD = false;
const double _THRESHOLD = 4.0;

/**
 * @brief _OUTPUT initiate as a pointer on textEdit attribute of the MainWindow class.
 */

bool _FORCED_EXCLUSION = true;

bool _WITH_FORWARD = false;

/** Number of extra column in Explored models file */
const int NBER_MIN_COL_EXPLORED_MODELS_FILE = 6;

/** End definition of global variables */

/** Basic functions --------------------*/

/**
 * @brief Sum
 */

template <class T> T sum(const T *x, int n)
{
  T s = (T) 0;
  for(int i = 0; i<n; i++)
    s += x[i];
  return s;
}

template<class T> T prod(const T *x, int n)
{
  T s = (T) 1;
  for(int i = 0; i<n; i++)
    s *= x[i];
  return s;
}

/**
 * @brief	Return the index of the maximum element
 */

template <class T> int whichMax(T const *vect, int n)
{
  int i;
  int indice=0;
  for(i=0; i<n; i++)
  {
    if(vect[i]>vect[indice])
      indice=i;
  }
  return indice;
}

// template <class T> int whichMax(vector<T> const &vect)
// {
//   int i;
//   int indice=0;
//   for(i=0; i<(int) vect.size(); i++)
//   {
//     if(vect[i]>vect[indice])
//       indice=i;
//   }
//   return indice;
// }


template <class T> int whichMax(T *first, T *last)
{
  int i=0, j=1;
  T *largest = first;
  if (first==last) return i;
  while (++first!=last)
  {
    if (*largest<*first)   // or: if (comp(*largest,*lowest)) for the comp version
    {
      largest=first;
      i=j;
    }
    j++;
  }
  return i;
}

template<class T> T myMax(T *x, int n)
{
  T xmax = x[0];
  for(int i = 1; i < n; i++)
	if(x[i]>xmax) xmax = x[i];
	
  return xmax;
}


/**
 * @brief	Return the index of the minimum element
 */

template <class T> int whichMin(const T *vect, int n)
{
  int indice = 0 ;
  for(int j=0; j<n; j++)
  {
    if(vect[j]<vect[indice])
      indice = j ;
  }
  return indice ;
}

// template <class T> int whichMin(vector<T> const &vect)
// {
//   int indice = 0 ;
//   for(int j=0; j<vect.size(); j++)
//   {
//     if(vect[j]<vect[indice])
//       indice = j ;
//   }
//   return indice ;
// }


template <class T> int whichMin(T *first, T *last)
{
  int i = 0, j=1;
  T *smallest = first;
  if (first==last) return i;
  while (++first!=last)
  {
    if (*smallest>*first)   // or: if (comp(*largest,*lowest)) for the comp version
    {
      smallest=first;
      i=j;
    }
    j++;
  }
  return i;
}


/**
 * @brief	Convert integer to string
 */
std::string itos(int i)
{
  std::stringstream s;
  s << i;
  return s.str();
}

/**
 * @brief	Simulate a simplex using uniform distribution
 */
//	[[Rcpp::export]]
Rcpp::DoubleVector simulProb(int n)
{
  double s;
  Rcpp::DoubleVector unif;
  unif = runif(n, 0, 1);// From Rmath.h
  s = std::accumulate(unif.begin(), unif.end(), 0.0);
  return unif/s;
}

/**
 * @brief factorial_recursion
 * @param n
 * @return
 */
int myFactorialRecursion(int n)
{
    if(n < 0)
        throw Rcpp::exception("Argument is a positive integer");
    return (n < 2)? 1 : n*myFactorialRecursion(n-1);
}

int myFactorial(int n)
{
    if(n < 0)
        throw Rcpp::exception("Argument is a positive integer");
    int fact = 1;
    while(n >= 1)
    {
        fact *= n;
        n--;
    }

    return fact;
}

//[[Rcpp::export]]
void testFactorial()
{
    Rcpp::Rcout << "\n" << myFactorial(0) << "\n";
    Rcpp::Rcout << "\n" << myFactorial(1) << "\n";
    Rcpp::Rcout << "\n" << myFactorial(2) << "\n";
    Rcpp::Rcout << "\n" << myFactorial(3) << "\n";
    Rcpp::Rcout << "\n" << myFactorial(4) << "\n";
    Rcpp::Rcout << "\n" << myFactorial(5) << "\n";
    Rcpp::Rcout << "\n" << myFactorial(6) << "\n";
    ErrorMemoryAlloc();
}

/**
 * @brief sample  Sample a vector 0:(K-1) using probability in proba
 * @param n
 * @param K
 * @param proba
 * @param Z
 */
void sample (int n, int K, const double *proba, int *Z)
{
  //  K = nombre de labels
  //  n = taille de l'echantillon a simuler
  //  proba=vecteur de proba
  std::vector<double> som_prob_vect(K+1);
  double *som_prob = som_prob_vect.data();
  int i, j;
  double s = 0;

  som_prob[0] = 0;
  for(i = 1; i < K+1 ; i++)
  {
    s += proba[i-1];
    som_prob[i] = s;
  }

  Rcpp::DoubleVector vect_proba = runif(n, 0.0, 1.0);

  for(i=0; i<n; i++)
  {
    for(j = 0 ; j < K ; j++)
    {
      if((vect_proba[i]>som_prob[j]) && (vect_proba[i]<=som_prob[j+1]))
        Z[i] = j;
    }
  }
}

/**
 * @brief stochastique  Sample data using probabilities in Tik
 * @param n_lig
 * @param K
 * @param Tik
 * @param vect
 */
void stochastique(int n_lig, int K, const double *Tik, int *vect)
{
  int i;
  for(i=0; i<n_lig; i++)
    sample(1, K, Tik + (i*K), (vect+i));
}

/**
 * @brief cutInMiddle Cut a std::string in 2 strings of the same length.
 * @param gen
 * @param a1
 * @param a2
 * @return
 */
bool cutInMiddle(std::string gen, std::string &a1, std::string &a2 )
{
  int long_genotype = gen.size();
  if(long_genotype %2 !=0)
  {
      MyError("the length of the string to cut is not a multiple of 2");
      return false;
  }
  int long_all = long_genotype / 2;
  a1.assign(gen, 0, long_all);
  a2.assign(gen, long_all, long_all);
  return true;
}

/**
 * @brief cutInN_Cpp Cut a std::string in N strings of the same length.
 * @param x A std::string
 * @param N Number of the output strings.
 * @param vect Vector that will containts the output strings.
 * @return
 */
bool cutInN_Cpp(std::string x, int N, std::string *vect)
{
    if(N < 1)
    {
        MyError("Not positive desired number of strings");
        return false;
    }

    if(x.size() % N != 0)
    {
        MyError("the length of the string to cut is not a multiple of N");
        return false;
    }

    int sizeOne = x.size()/N, i, i0;
    for(i = 0; i < N; i++)
    {
        i0 = i*sizeOne;
        vect[i].assign(x, i0, sizeOne);
    }

    return true;
}

/**
 * @brief cutInN
 * @param x A vector of strings.
 * @param N
 * @return
 */
// [[Rcpp::export]]
Rcpp::CharacterMatrix cutInN(Rcpp::CharacterVector x, int N)
{
    if(N < 1)
    {
        MyError("Not positive desired number of strings");
        return false;
    }

    std::string xx(Rcpp::as<std::string>(x[0]));
    int n_xx = xx.size();

    if(xx.size() % N != 0)
    {
        MyError("incompatible length of strings and N");
        return false;
    }

    int sizeOne = xx.size()/N;
    int nx = x.length(), i0, i, j;
    std::string xxx;
    Rcpp::CharacterMatrix M(nx, N);

    for(i = 0; i < nx; i++)
    {
        xx = Rcpp::as<std::string>(x[i]);
        if(((int) xx.size()) != n_xx)
        {
            MyError("Incompatible length");
            throw Rcpp::exception("Verifie that all string has the same length");
        }
        for(j = 0; j < N; j++)
        {
            i0 = j*sizeOne;
            xxx.assign(xx, i0, sizeOne);
            M(i, j) = xxx;
        }
    }

    return M;
}

/**
 * @brief cutEachColInN
 * @param tab
 * @param N
 * @return
 */
// [[Rcpp::export]]
Rcpp::CharacterMatrix cutEachColInN(Rcpp::CharacterMatrix tab, int N)
{
    int P = N*tab.ncol(), j, L = tab.nrow(), k;
    Rcpp::CharacterMatrix M(L, P);
    Rcpp::CharacterMatrix Mj(L, N);
    Rcpp::CharacterVector vect(L);

    for(j = 0; j < tab.ncol(); j++)
    {
        for(k = 0; k < L; k++)
            vect[k] = tab(k, j);
        Mj = cutInN(vect, N);
        for(k = 0; k < N; k++)
            M(_, j*N + k) = Mj(_, k);
    }

    return M;
}

/**
 * @brief howmanyWords
 * @param line
 * @return
 */
// [[Rcpp::export]]
int howmanyWords(std::string line)
{
  int nb_mots = 0;
  std::istringstream tampon (line);
  std::string mot;
  while(!tampon.eof())
  {
    tampon >> mot;
    if(tampon.bad() || tampon.fail())
    {
      // l'erreur peut survenir si la ligne finit par un espace : alors il n'y a plus de mot suivant
        if(!tampon.eof())
        {
            MyError("Error while counting words !");
            return 0;
        }
      break;
    }
    else
      nb_mots++;
  }
  return nb_mots;
}

/**
 * @brief isComment Renvoie true si la ligne commence par un caractère spécial.
 * @param ligne
 * @return
 */
// [[Rcpp::export]]
bool isComment(std::string ligne)
{
  if(howmanyWords(ligne)<=0) return false;
  std::string s;
  std::istringstream tampon (ligne);
  tampon >> s; //recuperer le premier mot en se debarrassant des espaces qui le precedent
  // tester si c'est un commentaire, qui commence par //, #, >, *
  if((s.size()>=2 && s[0] == '/' && s[1] == '/') ||
          (s.size()>=2 && s[0] == '/' && s[1] == '*') ||
          (s.size()>=1 && (s[0] == '#' || s[0] == '>' ||
                           s[0] == '*'))) return true;
  else return false;
}

/**
 * @brief nextLine
 * @param str
 * @param ligne
 * @param accept_empty
 * @return
 */
bool nextLine(std::istream &str, std::string &ligne, bool accept_empty = false)
{
  while(getline(str, ligne))
  {
    // si accept_empty, on garde les lignes vides
    if(accept_empty && howmanyWords(ligne)<=0) return true;
    if(howmanyWords(ligne)>0 && !isComment(ligne)) return true;
  }
  // s'il y a eu une erreur ou si on est arrive au bout
  return false;
}

/**
 * @brief lastLineComment
 * @param fichier
 * @return
 */
bool lastLineComment(std::string fichier)
{
  std::string ligne="";
  int n=0;
  std::ifstream fp(fichier.c_str(), std::ios::in);
  while(getline(fp, ligne)) n++; // aller a la fin du fichier
  if(n==0)
  {
    // fichier vide ou inexistant
    return false;
  }
  if(howmanyWords(ligne)<=0 || isComment(ligne))
  {
    // empty line or commented line
    return true;
  }
  
  fp.close();
  return false;
}

/**
 * @brief extractWords
 * @param line
 * @param words
 * @param line_length
 * @return
 */
// int extractWords(std::string &line, std::string *words, int line_length)
// {
//   // retourne le nombre de mot de la ligne, et enregistre les ncol premiers mots dans le tableau mots
//   std::istringstream tampon (line);
//   int nb_mots=0;
//   std::string mot;
//   while(!tampon.eof())
//   {
//     tampon >> mot;
//     if(tampon.bad() || tampon.fail())
//     {
//         MyError("while cutting a word");
//         return 0;
//         break;
//     }
// 
//     if(nb_mots<line_length)
//     {
//       words[nb_mots]=mot;
//       //clog << "word: "<< mots[nb_mots] << ", length: " << mots[nb_mots].size() << endl;
//     }
//     nb_mots++;
//   }
//   return nb_mots ;
// }
// 
// /**
//  * @brief extractWords_R
//  * @param line
//  * @return
//  */
// // [[Rcpp::export]]
// Rcpp::CharacterVector extractWords_R(std::string line)
// {
//     int n = howmanyWords(line);
//     //std::string words_0[n];
// 	std::vector<std::string> words_0(n);
//     int m = extractWords(line, words_0, n);
//     Rcpp::CharacterVector words(m);
//     if(m > 0)
//     {
//         for(int i = 0; i < m; i++)
//             words[i] = words_0[i];
//     }
// 
//     return words;
// }

/**
 * @brief nberOfLines
 * @param fichier
 * @return
 */
// [[Rcpp::export]]
int nberOfLines (std::string fichier)
{
  int n=0;
  std::string ligne;

  std::ifstream fp(fichier.c_str(), std::ios::in);
  while(nextLine(fp, ligne))
  {
    n++;
  }
  fp.close();
  return n;
}

/**
 * @brief nberOfColumns
 * @param fichier
 * @return
 */
// [[Rcpp::export]]
int nberOfColumns (std::string fichier)
{
  int i=0 ;
  int n = 0 ;
  std::string ligne_lue;

  std::ifstream fp(fichier.c_str(), std::ios::in);
  while(nextLine(fp, ligne_lue))
  {
    if(i == 0)
    {
      n = howmanyWords(ligne_lue) ; // Nombre de colonnes
    }
    if(n != howmanyWords(ligne_lue))
    {
      MyError("Incomplete line");
      fp.close();
      return 0 ;
    }
    i++ ;
  }
  fp.close();
  return n ;
}

/**
 * @brief readUntil
 * @param is
 * @param word
 * @param line
 * @param with_error
 * @param error_message
 * @return
 */
bool readUntil(std::istream &is, std::string word, std::string & line, bool with_error, std::string error_message)
{
  is.clear();
  is.seekg (0, std::ios::beg);
  if(!is.good())
  {
      MyError("reading stream");
      return false;
  }
  std::string lu, lu_1;
  while(nextLine(is, lu))
  {
    std::istringstream tampon (lu);
    lu_1.erase();
    tampon >> lu_1;
    if(lu_1==word)
    {
      line.erase();
      line=lu;
      return true;
    }
    //lu.erase();
  }
  if(with_error) Rcpp::Rcout << error_message << "not found word " << word << " in stream.\n";
  return false;
}

/**
 * @brief readLineN
 * @param ligne
 * @param fichier
 * @param n
 * @return
 */
bool readLineN(std::string fichier, int n, std::string & ligne)
{
  if(n<0 || n>=nberOfLines(fichier))
  {
      ErrorIndexOutOfRange();
      return false;
  }

  std::ifstream fp(fichier.c_str(), std::ios_base::in); 
  
  if(fp.bad())
  {
      ErrorOpeningFile();
      return false;
  }
  
//   try 
//   {
// 		fp.open(fichier.c_str(), std::ios_base::in);
// 		fp.exceptions(fp.failbit);
//   } catch (const std::ios_base::failure& e)
//   {
// 		Rcpp::Rcout << "Caught an ios_base::failure.\n"
//                   << "Explanatory string: " << e.what() << '\n';
// 		return false;
//   }
  
  std::string ligne_lue;
  for(int i=0; i<=n; i++)
      if(!nextLine(fp, ligne_lue))
      {
          MyError("cannot read line some line ");
          return false;
      }
  ligne=ligne_lue;
  fp.close();
  return true;
}

/**
 * @brief readLineN_R
 * @param fichier
 * @param n
 * @return
 */
// [[Rcpp::export]]
std::string readLineN_R(std::string fichier, int n)
{
    std::string ligne;
    readLineN(fichier, n, ligne);
    return ligne;
}



/** End definition of basic functions */

/** Set EM and selection paramters */

/**
 * @brief	Set EM settings
 * @code	Use EmSettings() for default settings
 */
// [[Rcpp::export()]]
void initialiseEmSettings()
{
  _EPSI = EPSI;
  _NBER_SMALL_EM = NBER_SMALL_EM;
  _NBER_ITER_EM = NBER_ITER_EM;
  _NBER_ITER_LONG_EM = NBER_ITER_LONG_EM;
  _TYPE_SMALL_EM = TYPE_SMALL_EM;
  _TYPE_EM = TYPE_EM;
  _PutTHRESHOLD = false;
}


// [[Rcpp::export()]]
void EmOptionsDefault()
{
  _EPSI = EPSI;
  _NBER_SMALL_EM = NBER_SMALL_EM;
  _NBER_ITER_EM = NBER_ITER_EM;
  _NBER_ITER_LONG_EM = NBER_ITER_LONG_EM;
  _TYPE_SMALL_EM = TYPE_SMALL_EM;
  _TYPE_EM = TYPE_EM;
  _PutTHRESHOLD = false;
}



/**
 * @brief EmSettings
 * @param epsi
 * @param nberSmallEM
 * @param nberIerEM
 * @param typeEM
 * @param putThreshold
 * @param launchAnalisys
 */
template<class T> bool EmSettings_conformity(T param, T minParam, T maxParam)
{
    if(param < minParam || param > maxParam)
    {
        Rcpp::Rcout << "\n > Give EM parameter in [ " << minParam << ", " << maxParam << " ]\n";
        return false;
    }

    return true;
}

/** FIXME Error message */
//	[[Rcpp::export]]
void EmSettings(double xepsi = -1.0,
                int xnberSmallEM = -1,
                int xnberIterations = -1,
                int xtypeEM = -1,
                int xtypeSmallEM = -1,
                int xnberIterLongEM = -1,
                bool xputThreshold = false)
{
    if(xepsi == -1.0) xepsi = EPSI;
    if(xnberSmallEM == -1) xnberSmallEM = NBER_SMALL_EM;
    if(xnberIterations == -1) xnberIterations = NBER_ITER_EM;
    if(xnberIterLongEM == -1) xnberIterLongEM = NBER_ITER_LONG_EM;
    if(xtypeEM == -1) xtypeEM = TYPE_EM;
    if(xtypeSmallEM == -1) xtypeSmallEM = TYPE_SMALL_EM;

    if(!EmSettings_conformity(xepsi, MIN_EPSI, MAX_EPSI) ||
            !EmSettings_conformity(xnberSmallEM, MIN_NBER_SMALL_EM, MAX_NBER_SMALL_EM) ||
            !EmSettings_conformity(xnberIterations, MIN_NBER_ITER_EM, MAX_NBER_ITER_EM) ||
            !EmSettings_conformity(xtypeEM, MIN_TYPE_EM, MAX_TYPE_EM) ||
            !EmSettings_conformity(xtypeSmallEM, MIN_TYPE_SMALL_EM, MAX_TYPE_SMALL_EM) ||
            !EmSettings_conformity(xnberIterLongEM, MIN_NBER_ITER_LONG_EM, MAX_NBER_ITER_LONG_EM))
    {
        initialiseEmSettings();
        MyWarning("Some of the EM options are out of their range; default options were considered");
        return;
    }

    _EPSI = xepsi;
    _NBER_SMALL_EM = xnberSmallEM;
    _NBER_ITER_EM = xnberIterations;
    _TYPE_EM = xtypeEM;
    _TYPE_SMALL_EM = xtypeSmallEM;
    _NBER_ITER_LONG_EM = xnberIterLongEM;
    _PutTHRESHOLD = xputThreshold;
}

// [[Rcpp::export]]
void EmOptionsDisplay()
{
    Rcpp::Rcout << "\n > EPSI = " << _EPSI;
    Rcpp::Rcout << "\n > NBER_SMALL_EM = " << _NBER_SMALL_EM;
    Rcpp::Rcout << "\n > NBER_ITERATIONS_EM = " << _NBER_ITER_EM;
    Rcpp::Rcout << "\n > NBER_ITER_LONG_EM = " << _NBER_ITER_LONG_EM;
    Rcpp::Rcout << "\n > TYPE_SMALL_EM = " << _TYPE_EM_NAMES[_TYPE_SMALL_EM];
    Rcpp::Rcout << "\n > TYPE_EM = " << _TYPE_EM_NAMES[_TYPE_EM];

	if(_PutTHRESHOLD)
	  Rcpp::Rcout << "\n > Put THRESHOLD = TRUE";
	else Rcpp::Rcout << "\n > Put THRESHOLD = FALSE" ;
    Rcpp::Rcout <<"\n";
}

// [[Rcpp::export()]]
Rcpp::List getEmOptions_Rcpp()
{
    Rcpp::List outList = Rcpp::List::create(_["epsi"] = _EPSI,
            _["nberSmallEM"] = _NBER_SMALL_EM,
            _["nberIterations"] = _NBER_ITER_EM,
            _["typeSmallEM"] = _TYPE_SMALL_EM,
            _["typeEM"] = _TYPE_EM,
            _["nberMaxIterations"] = _NBER_ITER_LONG_EM,
            _["putThreshold"] = _PutTHRESHOLD);
    return outList;
}



/** End setting EM and selection parameters */

/**
 * @brief	Class PAR_KS
 * @Warning	TODO : Change the names of variables
 */

class PAR_KS
{
private:
  int _N_OfPAR_KS;
  int _K_OfPAR_KS;
  Rcpp::LogicalVector _S_OfPAR_KS;
  int _Dim;
  Rcpp::DoubleVector _PI_K;
  Rcpp::NumericMatrix _PROB;
  double _LOG_LIK;
  Rcpp::NumericMatrix _Tik;
  Rcpp::IntegerVector _POST_CLASSIF;
  double _ENT;
  Rcpp::DoubleVector _CRITERIA;

  Rcpp::CharacterVector _LEVELS; // Set only with setLEVELS
  Rcpp::IntegerVector _N_LEVELS;

public:
  PAR_KS()
  {
    _N_OfPAR_KS = 0;
    _K_OfPAR_KS = 0;
    _Dim = 0;
    _LOG_LIK = 0.0;
    _ENT = 0.0;
  }

  PAR_KS(int N, int K, Rcpp::LogicalVector S)
  {
    _N_OfPAR_KS = N;
    _K_OfPAR_KS = K;
    _S_OfPAR_KS = S;
    _Dim = 0;
    _LOG_LIK = 0.0;
    _ENT = 0.0;
  }

  PAR_KS(int N, int K,
         Rcpp::LogicalVector S,
         Rcpp::IntegerVector n_levels,
         Rcpp::DoubleVector levels_freq)
  {
    randomInitialise(N, K, S, n_levels, levels_freq);
  }

  /**
   * @brief PAR_KS
   * @param par_ks
   */
  PAR_KS(const PAR_KS &par_ks)
  {
    _N_OfPAR_KS = par_ks._N_OfPAR_KS;
    _K_OfPAR_KS = par_ks._K_OfPAR_KS;
    _S_OfPAR_KS = par_ks._S_OfPAR_KS;
    _Dim = par_ks._Dim;
    _PI_K = par_ks._PI_K;
    _PROB = par_ks._PROB;
    _LOG_LIK = par_ks._LOG_LIK;
    _Tik = par_ks._Tik;
    _POST_CLASSIF = par_ks._POST_CLASSIF;
    _ENT = par_ks._ENT;
    _CRITERIA = par_ks._CRITERIA;

    _LEVELS = par_ks._LEVELS;
    _N_LEVELS = par_ks._N_LEVELS;
  }

  /**
   * @brief PAR_KS
   * @param K
   * @param S
   * @param pi_k
   * @param prob
   * @param n_levels
   * @param levels_freq
   */
  PAR_KS(int N, int K,
         Rcpp::LogicalVector S,
         Rcpp::DoubleVector pi_k,
         Rcpp::NumericMatrix prob,
         Rcpp::IntegerVector n_levels,
         Rcpp::DoubleVector levels_freq)
  {
    // Set _LEVELS in set()
    set(N, K, S, pi_k, prob, n_levels, levels_freq);
  }

  /**
   * @brief PAR_KS
   * @param listParKS
   * @warning _Tik, _POST_CLASSIF are not defined
   */
  PAR_KS(Rcpp::List listParKS)
  {
    _N_OfPAR_KS = Rcpp::as<int>(listParKS["N"]);

    _K_OfPAR_KS = Rcpp::as<int>(listParKS["K"]);

    _S_OfPAR_KS = Rcpp::as<Rcpp::LogicalVector>(listParKS["S"]);

    _N_LEVELS = Rcpp::as<Rcpp::IntegerVector>(listParKS["NbersLevels"]);
    _Dim = Rcpp::as<int>(listParKS["dim"]);
    _PI_K = Rcpp::as<Rcpp::DoubleVector>(listParKS["pi_K"]);
    int j, a0, P = _S_OfPAR_KS.length(), a, k, a1 = std::accumulate(_N_LEVELS.begin(), _N_LEVELS.end(), 0);

    Rcpp::NumericMatrix PROB(a1, _K_OfPAR_KS);
    Rcpp::List probList = Rcpp::as<Rcpp::List>(listParKS["prob"]);
    for(j = 0; j < P; j++)
      {
          a0 = std::accumulate(_N_LEVELS.begin(), _N_LEVELS.begin() + j, 0);
          Rcpp::NumericMatrix probj = Rcpp::as<Rcpp::NumericMatrix>(probList[j]);
          for(k = 0; k < _K_OfPAR_KS; k++)
              for(a = 0; a < _N_LEVELS[j]; a++)
                  PROB(a0 + a, k) = probj(a, k);
      }
    _PROB = PROB;

    setLEVELS_default();
    _LOG_LIK = Rcpp::as<double>(listParKS["logLik"]);
    _ENT = Rcpp::as<double>(listParKS["entropy"]);

    _CRITERIA = Rcpp::as<Rcpp::DoubleVector>(listParKS["criteria"]);
    _Tik = Rcpp::as<Rcpp::NumericMatrix>(listParKS["Tik"]);

      //MyDebug("3");
    _POST_CLASSIF = Rcpp::as<Rcpp::IntegerVector>(listParKS["mapClassif"]);
  }

  ~PAR_KS(){}

  /**
   * @brief operator =
   * @param par_ks
   * @return
   */
  PAR_KS & operator= (PAR_KS par_ks)
  {
    _N_OfPAR_KS = par_ks._N_OfPAR_KS;
    _K_OfPAR_KS = par_ks._K_OfPAR_KS;
    _S_OfPAR_KS = par_ks._S_OfPAR_KS;
    _Dim = par_ks._Dim;
    _PI_K = par_ks._PI_K;
    _PROB = par_ks._PROB;
    _LOG_LIK = par_ks._LOG_LIK;
    _Tik = par_ks._Tik;
    _POST_CLASSIF = par_ks._POST_CLASSIF;
    _ENT = par_ks._ENT;
    _CRITERIA = par_ks._CRITERIA;

    _LEVELS = par_ks._LEVELS;
    _N_LEVELS = par_ks._N_LEVELS;

    return *this;
  }

  /**
   * @brief setDim
   */
  void setDim()
  {
    int j, s, sc;
    s = 0.0;
    sc = 0.0;

    for(j = 0; j < _S_OfPAR_KS.length(); j++)
    {
      if(_S_OfPAR_KS[j]) s += _N_LEVELS[j] - 1;
      else sc += _N_LEVELS[j] - 1;
    }

    _Dim = _K_OfPAR_KS - 1 + _K_OfPAR_KS*s + sc;
  }

  /**
   * @brief set
   * @param K
   * @param S
   * @param pi_k
   * @param prob
   * @param n_levels
   * @param levels_freq
   */
  void set(int N, int K,
           Rcpp::LogicalVector S,
           Rcpp::DoubleVector pi_k,
           Rcpp::NumericMatrix prob,
           Rcpp::IntegerVector n_levels,
           Rcpp::DoubleVector levels_freq)
  {
    _N_OfPAR_KS = N;
    _K_OfPAR_KS = K;
    _S_OfPAR_KS = S;
    _PI_K = pi_k;
    _PROB = prob;
    _N_LEVELS = n_levels;
    _LOG_LIK = 0.0;
    _ENT = 0.0;

    int j, k = 0, a0, a;

    for(j = 0; j < _S_OfPAR_KS.length(); j++)
      {
          if(!_S_OfPAR_KS[j])
          {
              a0 = std::accumulate(_N_LEVELS.begin(), _N_LEVELS.begin() + j, 0.0);
              for(a = 0; a < _N_LEVELS[j]; a++)
                  for(k = 0; k < _K_OfPAR_KS; k++)
                      _PROB(a0 + a, k) = levels_freq[a0 + a];
          }
      }

    setLEVELS_default();
    Rcpp::CharacterVector colNames(K);
    for(k = 0; k<K; k++)
        colNames[k] = "pop" + itos(k+1);
    _PROB.attr("dimnames") = Rcpp::List::create(_LEVELS, colNames);
    setDim();
  }

  /**
   * @brief setFromList
   * @param listParKS
   */
  void setFromList(Rcpp::List listParKS)
  {
      _N_OfPAR_KS = Rcpp::as<int>(listParKS["N"]);
      _K_OfPAR_KS = Rcpp::as<int>(listParKS["K"]);
      _S_OfPAR_KS = Rcpp::as<Rcpp::LogicalVector>(listParKS["S"]);

      _N_LEVELS = Rcpp::as<Rcpp::IntegerVector>(listParKS["NbersLevels"]);
      _Dim = Rcpp::as<int>(listParKS["dim"]);
      _PI_K = Rcpp::as<Rcpp::DoubleVector>(listParKS["pi_K"]);
      int j, a0, P = _S_OfPAR_KS.length(), a, k, a1 = std::accumulate(_N_LEVELS.begin(), _N_LEVELS.end(), 0);

      Rcpp::NumericMatrix PROB(a1, _K_OfPAR_KS);
      for(j = 0; j < P; j++)
      {
          a0 = std::accumulate(_N_LEVELS.begin(), _N_LEVELS.begin() + j, 0);
          Rcpp::List probList = Rcpp::as<Rcpp::List>(listParKS["prob"]);
          Rcpp::NumericMatrix probj = Rcpp::as<Rcpp::NumericMatrix>(probList[j]);
          for(k = 0; k < _K_OfPAR_KS; k++)
              for(a = 0; a < _N_LEVELS[j]; a++)
                  PROB(a0 + a, k) = probj(a, k);
      }
      _PROB = PROB;

      setLEVELS_default();
      _LOG_LIK = Rcpp::as<double>(listParKS["logLik"]);
      _ENT = Rcpp::as<double>(listParKS["entropy"]);

      _CRITERIA = Rcpp::as<Rcpp::DoubleVector>(listParKS["criteria"]);
      _Tik = Rcpp::as<Rcpp::NumericMatrix>(listParKS["Tik"]);
      _POST_CLASSIF = Rcpp::as<Rcpp::IntegerVector>(listParKS["mapClassif"]);
  }

  void setN(int N){_N_OfPAR_KS = N;}

  /**
   * @brief setK
   * @param K
   */
  void setK(int K){_K_OfPAR_KS = K;}

  /**
   * @brief setS
   * @param S
   */
  void setS(Rcpp::LogicalVector S){_S_OfPAR_KS = S;}

  /**
   * @brief setKS
   * @param K
   * @param S
   */
  void setKS(int K, Rcpp::LogicalVector S)
  {
      _K_OfPAR_KS = K;
      _S_OfPAR_KS = S;
  }

  void setNPKS(int N, int K, Rcpp::LogicalVector S)
  {
    _N_OfPAR_KS = N;
    _K_OfPAR_KS = K;
    _S_OfPAR_KS = S;
  }

  /**
   * @brief setPI_K
   * @param pi_k
   */
  void setPI_K(Rcpp::DoubleVector pi_k){_PI_K = pi_k;}

  /**
   * @brief setPI_k
   * @param p
   * @param k
   */
  void setPI_k(double p, int k) {_PI_K[k] = p;}

  /**
   * @brief setPROB
   * @param prob
   */
  void setPROB(Rcpp::NumericMatrix prob){_PROB = prob;}

  /**
   * @brief setPROBkja
   * @param p
   * @param k
   * @param j
   * @param a
   */
  void setPROBkja(double p, int k, int j, int a)
  {
    int a0 = std::accumulate(_N_LEVELS.begin(), _N_LEVELS.begin() + j, 0);
    _PROB(a0 + a, k) = p;
  }

  void setLOG_LIK(double log_lik){_LOG_LIK = log_lik;}

  void setTik(Rcpp::NumericMatrix Tik){_Tik = Tik;}

  void setPOST_CLASSIF(Rcpp::IntegerVector classif){_POST_CLASSIF = classif;}

  void setCRITERIA(double lv, double Cte)
  {
    if(_N_OfPAR_KS == 0)
      {
        throw Rcpp::exception("N is equal to 0");
        return;
      }
    double dim = (double) _Dim;
    _CRITERIA = Rcpp::DoubleVector::create(_["BIC"] = -lv + ((double) dim*log((double) _N_OfPAR_KS))/2,
                                           _["AIC"] = -lv + dim,
                                        _["ICL"] = -lv + (dim*log((double) _N_OfPAR_KS))/2 + _ENT,
                                           _["CteDim"] = -lv + Cte*dim);
  }

  void setLEVELS(Rcpp::CharacterVector levels) {_LEVELS = levels;}

  void setLEVELS_default()
  {
      int j, a, a1 = std::accumulate(_N_LEVELS.begin(), _N_LEVELS.end(), 0);
      Rcpp::CharacterVector vect(a1);
      int k = 0;
      for(j = 0; j < _N_LEVELS.length(); j++)
          for(a = 0; a < _N_LEVELS[j]; a++)
          {
              vect[k] = a+1;
              k++;
          }
      _LEVELS = vect;
  }

  void setN_LEVELS(Rcpp::IntegerVector n_levels) {_N_LEVELS = n_levels;}

  void setTik2(double *Tik)
  {
    _Tik = Rcpp::NumericMatrix(_N_OfPAR_KS, _K_OfPAR_KS);
    int i, k;
    for(i = 0; i < _N_OfPAR_KS; i++)
      for(k = 0; k < _K_OfPAR_KS; k++)
        _Tik(i, k) = Tik[i*_K_OfPAR_KS + k];
  }

  //void setENT(double ent){_ENT = ent;}

  void setENT()
  {
    int i, k;
    _ENT = 0.0;
    for(i = 0; i < _N_OfPAR_KS; i++)
        for(k = 0; k < _K_OfPAR_KS; k++)
            _ENT -= _Tik[i*_K_OfPAR_KS+k] <= 0? 0 : _Tik[i*_K_OfPAR_KS+k]*log(_Tik[i*_K_OfPAR_KS+k]);
  }


  /**
   * @brief	Getters
   */
  int getN(){return _N_OfPAR_KS;}

  int getK() {return _K_OfPAR_KS;}

  bool getS_j(int j) {return _S_OfPAR_KS[j];}

  Rcpp::LogicalVector getS() {return _S_OfPAR_KS;}

  int getDim() {return _Dim;}

  Rcpp::DoubleVector getPI_K() {return _PI_K;}

  Rcpp::NumericMatrix getPROB() {return _PROB;}

  double getLOG_LIK() {return _LOG_LIK;}

  Rcpp::NumericMatrix getTik() {return _Tik;}

  Rcpp::IntegerVector getPOST_CLASSIF() {return _POST_CLASSIF;}

  double getENT() {return _ENT;}

  Rcpp::CharacterVector getLEVELS() {return _LEVELS;}

  Rcpp::IntegerVector getN_LEVELS() {return _N_LEVELS;}

  /**
   * @brief	Get a product of Probabilities corresponding to a vector
   * @warning	FIXME in case of _PLOIDY > 2
   */
  double getProdProb1(int *x, int ploidy, int k, int j)
  {
    double p = 1.0;
    int a, a0 = std::accumulate(_N_LEVELS.begin(), _N_LEVELS.begin() + j, 0);
    for(a = 0; a < ploidy; a++)
      p *= _PROB(a0 + x[a], k);

    double b;
    if(ploidy == 2)
      b = (x[0] == x[1])? 1.0 : 2.0;
    else b = 1;

    return b*p;
  }

  /**
   * @brief getProdProb For all value of _PLOIDY
   * @param x
   * @param ploidy
   * @param k
   * @param j
   * @return
   */
  double getProdProb(int *x, int ploidy, int k, int j)
  {
    double p = 1.0;
    int i, a, a0 = std::accumulate(_N_LEVELS.begin(), _N_LEVELS.begin() + j, 0);
    for(a = 0; a < ploidy; a++)
      p *= _PROB(a0 + x[a], k);

    double num = myFactorial(ploidy), deno = 1;
    int s;
    for(a = 0; a < _N_LEVELS[j]; a++)
    {
        s = 0;
        for(i = 0; i < ploidy; i++)
            if(a == x[i]) s += 1;
        deno *= myFactorial(s);
    }

    double b = ((double) num)/deno;

    return b*p;
  }


  double getPI_k(int k) {return _PI_K[k];}

  double getPROBkja(int k, int j, int a)
  {
    int a0 = std::accumulate(_N_LEVELS.begin(), _N_LEVELS.begin() + j, 0);
    return _PROB(a0 + a, k);
  }

  /**
   * @brief	Return PAR_KS like a list for R
   */
  Rcpp::List getList()
  {
    int i0, i, j, p = _S_OfPAR_KS.length();
    Rcpp::List R_out;


    Rcpp::List ListProb(p);
    Rcpp::List ListLevels(p);
    Rcpp::CharacterVector vect1(_K_OfPAR_KS);
    for(i = 0; i<_K_OfPAR_KS; i++)
    {
      vect1[i] = "pop" + itos(i+1);
    }

    for(i = 0; i < p; i++)
    {
        i0 = std::accumulate(_N_LEVELS.begin(), _N_LEVELS.begin() + i, 0);

        Rcpp::NumericMatrix Probx = _PROB(Rcpp::Range(i0, i0 + _N_LEVELS[i] - 1), Rcpp::Range(0, _K_OfPAR_KS-1));

        Rcpp::CharacterVector vect(_N_LEVELS[i]);
        for(j = 0; j<_N_LEVELS[i]; j++)
            vect[j] = _LEVELS[i0 + j];

        ListLevels[i] = vect;
        Probx.attr("dimnames") = Rcpp::List::create(vect, vect1);
        ListProb[i] = Probx;
    }

    R_out["N"] = _N_OfPAR_KS;
    R_out["K"] = _K_OfPAR_KS;
    R_out["S"] = _S_OfPAR_KS;
    R_out["dim"] = _Dim;

    R_out["pi_K"] = _PI_K;
    R_out["prob"] = ListProb;
    R_out["logLik"] = _LOG_LIK;
    R_out["entropy"] = _ENT;
    R_out["criteria"] = _CRITERIA;
    R_out["Tik"] = _Tik;
    R_out["mapClassif"] = _POST_CLASSIF;
    R_out["NbersLevels"] = _N_LEVELS;
    R_out["levels"] = ListLevels;

    return R_out;
  }

  /**
   * @brief randomInitialise
   * @param K
   * @param S
   * @param n_levels
   * @param levels_freq
   */
  void randomInitialise(int N, int K,
                        Rcpp::LogicalVector S,
                        Rcpp::IntegerVector n_levels,
                        Rcpp::DoubleVector levels_freq)
  {
	  _N_OfPAR_KS = N;
      _K_OfPAR_KS = K;
      _S_OfPAR_KS = S;
      _LOG_LIK = 0.0;
      _ENT = 0.0;

      _PI_K = Rcpp::DoubleVector(K);
      _PI_K.fill(((double) 1/K));
      _N_LEVELS = n_levels;

      int j, a0, k, a;

      _PROB = Rcpp::NumericMatrix(std::accumulate(n_levels.begin(), n_levels.end(), 0), K);
      for(j = 0; j < n_levels.length(); j++)
      {
        a0 = std::accumulate(n_levels.begin(), n_levels.begin() + j, 0);
        for(k = 0; k < K; k++)
        {
          Rcpp::DoubleVector prob = simulProb(n_levels[j]);
          for(a = 0; a < n_levels[j]; a++)
            _PROB(a0 + a, k) = prob[a];
        }
      }

      if(levels_freq.length() != 0)
      {
        for(j = 0; j < n_levels.length(); j++)
        {
          if(!S[j])
          {
            a0 = std::accumulate(n_levels.begin(), n_levels.begin() + j, 0);
            for(k = 0; k < K; k++)
            {
              for(a = 0; a < n_levels[j]; a++)
              _PROB(a0 + a, k) = levels_freq[a0 + a];
            }
          }
        }
      }

      setLEVELS_default();
  }

  /**
   * @brief	Print PAR_KS
   */
  void print()
  {
    Rcpp::Rcout.precision(PRECISION);

    int i, j, a, a0, k;
    Rcpp::Rcout << "\n> PAR_KS print method\n";
    Rcpp::Rcout << "\n> Size of data N = " << _N_OfPAR_KS << "\n";
    Rcpp::Rcout << "\tNumber of populations K = " << _K_OfPAR_KS << "\n";
    Rcpp::Rcout << "\tSelected variables S = " ;
    for(j = 0; j < _S_OfPAR_KS.length(); j++)
      Rcpp::Rcout << (_S_OfPAR_KS[j]? 1 : 0) << " ";
    Rcpp::Rcout << "\n";
    Rcpp::Rcout << "\tMixing proportions : ";
    for(j = 0; j < _K_OfPAR_KS; j++)
      Rcpp::Rcout << _PI_K[j] << " ";
    Rcpp::Rcout << "\n";
    Rcpp::Rcout << "\tNumbers levels : ";
    for(j = 0; j < _S_OfPAR_KS.length(); j++)
        Rcpp::Rcout << " " << _N_LEVELS[j] ;

    Rcpp::Rcout << "\n\tEstimates of probabilities in different populations\n";
    for(j = 0; j<_S_OfPAR_KS.length(); j++)
    {
      Rcpp::Rcout << "\t X" << j+1 << "\n";

      a0 = 0;
      for(i = 0; i<j ; i++)
        a0 += _N_LEVELS[i];

      for(a = 0; a<_N_LEVELS[j]; a++)
      {
        Rcpp::Rcout << "\t\t  " << _LEVELS[a0+a] << "\t";
        for(k = 0; k<_K_OfPAR_KS; k++)
          Rcpp::Rcout << std::fixed << _PROB(a0 + a, k) << "\t";
        Rcpp::Rcout << "\n";
      }
    }
    Rcpp::Rcout << "\tDimension = " << _Dim << "\n";
    Rcpp::Rcout << "\tLog-likelihood = " << std::fixed << _LOG_LIK << "\n";
    Rcpp::Rcout << "\tEntropy = " << std::fixed << _ENT << "\n";
  }
  
  /**
   * @brief	Print PAR_KS
   */
  void printInFile(std::string file)
  {
		std::ofstream out ;
		out.open(file.c_str(), std::ios_base::out);
	
    out.precision(PRECISION);

    int i, j, a, a0, k;
    out << "#The set of parameters of a model (K,S) by the R package ClustMMDD\n";
    out << "#Size of data N = " << _N_OfPAR_KS << "\n";
    out << "K " << _K_OfPAR_KS << "\n";
    out << "S " ;
    for(j = 0; j < _S_OfPAR_KS.length(); j++)
      out << (_S_OfPAR_KS[j]? 1 : 0) << " ";
    out << "\n";
    out << "\tMixing proportions : ";
    for(j = 0; j < _K_OfPAR_KS; j++)
      out << _PI_K[j] << " ";
    out << "\n";
    out << "\tNumbers levels : ";
    for(j = 0; j < _S_OfPAR_KS.length(); j++)
        out << " " << _N_LEVELS[j] ;

    out << "\n\tEstimates of probabilities in different populations\n";
    for(j = 0; j<_S_OfPAR_KS.length(); j++)
    {
      out << "\t X" << j+1 << "\n";

      a0 = 0;
      for(i = 0; i<j ; i++)
        a0 += _N_LEVELS[i];

      for(a = 0; a<_N_LEVELS[j]; a++)
      {
        out << "\t\t  " << _LEVELS[a0+a] << "\t";
        for(k = 0; k<_K_OfPAR_KS; k++)
          out << std::fixed << _PROB(a0 + a, k) << "\t";
        out << "\n";
      }
    }
    out << "\tDimension = " << _Dim << "\n";
    out << "\tLog-likelihood = " << std::fixed << _LOG_LIK << "\n";
    out << "\tEntropy = " << std::fixed << _ENT << "\n";
	out << "END\n";
	
	out.close();
  }
  
  
  /**
   * @brief writeModelInFile
   * @param xfile
   */
  void writeModelInFile(std::string xfile)
  {
	std::ofstream out ;
	out.open(xfile.c_str(), std::ios_base::app);
	int P = _N_LEVELS.length();
	out.precision(PRECISION);
    out << _N_OfPAR_KS << " " << P << " " ;
	for(int i = 0; i < P; i++)
	  if(_S_OfPAR_KS[i]) out << 1 << " ";
	  else out << 0 << " ";
	out  << _LOG_LIK << " " << _Dim << " " << _ENT << "\n";
  
	out.close();
  }

};

/**
 * @brief	Module PAR_KS
 */
RCPP_MODULE(MODULE_PAR_KS)
{
  using namespace Rcpp;
  
  Rcpp::class_<PAR_KS>("PAR_KS")
  .constructor()
  .constructor<int, int, Rcpp::LogicalVector>()
  .constructor<int, int, Rcpp::LogicalVector,
          Rcpp::IntegerVector,
          Rcpp::DoubleVector>()
  .constructor<int, int,
          Rcpp::LogicalVector,
          Rcpp::DoubleVector,
          Rcpp::NumericMatrix,
          Rcpp::IntegerVector,
          Rcpp::DoubleVector>()
  .constructor<Rcpp::List>()

  .method("set", &PAR_KS::set)
  .method("setFromList", &PAR_KS::setFromList)
  .method("setN", &PAR_KS::setN)
  .method("setK", &PAR_KS::setK)
  .method("setS", &PAR_KS::setS)
  .method("setDim", &PAR_KS::setDim)
  .method("setPI_K", &PAR_KS::setPI_K)
  .method("setPROB", &PAR_KS::setPROB)
  .method("setLOG_LIK", &PAR_KS::setLOG_LIK)
  .method("setTik", &PAR_KS::setTik)
  .method("setPOST_CLASSIF", &PAR_KS::setPOST_CLASSIF)
  .method("setENT1", &PAR_KS::setENT)

  .method("getN", &PAR_KS::getN)
  .method("getK", &PAR_KS::getK)
  .method("getS", &PAR_KS::getS)
  .method("getDim", &PAR_KS::getDim)
  .method("getPI_K", &PAR_KS::getPI_K)
  .method("getPROB", &PAR_KS::getPROB)
  .method("getLOG_LIK", &PAR_KS::getLOG_LIK)
  .method("getTik", &PAR_KS::getTik)
  .method("getPOST_CLASSIF", &PAR_KS::getPOST_CLASSIF)
  .method("getENT1", &PAR_KS::getENT)
  .method("getList", &PAR_KS::getList)
  .method("randomInitialise", &PAR_KS::randomInitialise)

  .method("print", &PAR_KS::print)
  .method("writeModelInFile", &PAR_KS::writeModelInFile)
  ;
}


/**
 * @author	Wilson Toussile
 * @brief	Class DATA
 */
class DATA
{
protected:
    int _N_OfDATA;
    int _P_OfDATA;
    int _PLOIDY;
    int *_DATA;
    std::string *_LEVELS;
    int *_N_LEVELS;
    int  *_LEVELS_COUNT;
    double *_LEVELS_FREQ;
    int *_PRIOR_CLASSIF;

  void initialize()
  {
      _N_OfDATA = 0;
      _P_OfDATA = 0;
      _PLOIDY = 0;
      _DATA = NULL;
      _LEVELS = NULL;
      _N_LEVELS = NULL;
      _LEVELS = NULL;
      _LEVELS_COUNT = NULL;
      _LEVELS_FREQ = NULL;
      _PRIOR_CLASSIF = NULL;
  }

public:
  std::string name;

  DATA()
  {
      _N_OfDATA = 0;
      _P_OfDATA = 0;
      _PLOIDY = 0;
      _DATA = NULL;
      _LEVELS = NULL;
      _N_LEVELS = NULL;
      _LEVELS = NULL;
      _LEVELS_COUNT = NULL;
      _LEVELS_FREQ = NULL;
      _PRIOR_CLASSIF = NULL;
  }

  /**
   * @brief DATA
   * @param data
   * @param ploidy
   */
  DATA(Rcpp::IntegerMatrix data,
       int ploidy)
  {
      if(ploidy < 1 || ploidy > NBER_OCCURRENCES_MAX || data.nrow() % ploidy != 0)
          throw Rcpp::exception("Incompatible dimension or number of occurrences incorrect");
      _DATA = data.begin();
      _N_OfDATA = data.ncol(); // Warning : data is transposed and stored by row
      _P_OfDATA = data.nrow()/ploidy;
      _LEVELS = NULL;
      _N_LEVELS = NULL;
      _LEVELS_COUNT = NULL;
      _LEVELS_FREQ = NULL;
      _PRIOR_CLASSIF = NULL;
  }

  DATA(Rcpp::IntegerMatrix data,
       int ploidy,
       Rcpp::CharacterVector levels,
       Rcpp::IntegerVector n_levels,
       Rcpp::IntegerVector levels_count,
       Rcpp::DoubleVector levels_freq)
  {
      if(ploidy < 1 || ploidy > NBER_OCCURRENCES_MAX || data.nrow() % ploidy != 0)
          throw Rcpp::exception("Incompatible dimension or number of occurrences incorrect");
      _DATA = data.begin();
      _N_OfDATA = data.ncol();
      _PLOIDY = ploidy;
      _P_OfDATA = data.nrow()/_PLOIDY;
      //_LEVELS = levels.begin();
      _N_LEVELS = n_levels.begin();
      _LEVELS_COUNT = levels_count.begin();
      _LEVELS_FREQ = levels_freq.begin();
      _PRIOR_CLASSIF = NULL;
  }

  DATA(Rcpp::IntegerMatrix data,
       int ploidy,
       Rcpp::IntegerVector n_levels,
      Rcpp::DoubleVector levels_freq)
  {
      if(ploidy < 1 || ploidy > NBER_OCCURRENCES_MAX || data.nrow() % ploidy != 0)
          throw Rcpp::exception("Incompatible dimension or number of occurrences incorrect");
      _DATA = data.begin();
      _N_OfDATA = data.ncol();
      _P_OfDATA = data.nrow()/ploidy;
      _PLOIDY = ploidy;
      //_LEVELS = levels.begin();
      _N_LEVELS = n_levels.begin();
      //_LEVELS_COUNT = levels_count.begin();
      _LEVELS_FREQ = levels_freq.begin();
      _PRIOR_CLASSIF = NULL;
  }

  DATA(DATA &data){}

  ~DATA(){}

  /**
   * @brief	Get pointer on the first data of individual i
   */
  int * operator() (int i, Rcpp::internal::NamedPlaceHolder)
  {
    if(i >= _N_OfDATA)
    {
        MyError("Index out of bounds");
        throw Rcpp::exception("Index out of bounds");
    }
    return _DATA + i*_PLOIDY*_P_OfDATA;
  }

  /**
   * @brief	Get pointer on the firts data of individual i at variable j.
   */
  int * operator() (int i, int j)
  {
    if(i >= _N_OfDATA || j >= _P_OfDATA)
    {
        MyError("Index out of bounds");
        throw Rcpp::exception("Index out of bounds");
    }
    return _DATA + (_P_OfDATA*i + j)*_PLOIDY;
  }

  /**
   * @brief	Get the number of a in data of individual i at variable j?
   */
  int howMany(int a, int i, int j)
  {
    int s = 0;
    for(int k = 0; k < _PLOIDY; k++)
      s += (*(_DATA + (_P_OfDATA*i + j)*_PLOIDY+ k) == a)? 1 : 0;

    return s;
  }

  int * getDATA() {return _DATA;}

  int getN() {return _N_OfDATA;}

  int getP() {return _P_OfDATA;}

  std::string * getLEVELS() {return _LEVELS;}

  int *getN_LEVELS() {return _N_LEVELS;}

  double *getLEVELS_FREQ() {return _LEVELS_FREQ;}

  int *getLEVELS_COUNT() {return _LEVELS_COUNT;}

  int *getPRIOR_CLASSIF() {return _PRIOR_CLASSIF;}

  int getN_OCCURRENCES() {return _PLOIDY;}

  void setNAME(std::string new_name) {name = new_name;}

  void setDATA(Rcpp::IntegerMatrix mat)
  {
      // Summit a transposed matrix
      _N_OfDATA = mat.ncol();
      _P_OfDATA = mat.nrow();
      _DATA = mat.begin();
  }

  void set(Rcpp::IntegerMatrix data,
         int ploidy,
         Rcpp::IntegerVector n_levels,
         Rcpp::IntegerVector levels_count,
         Rcpp::DoubleVector levels_freq)
    {
        if(ploidy < 1 || ploidy > NBER_OCCURRENCES_MAX || data.nrow() % ploidy != 0)
            throw Rcpp::exception("Incompatible dimension or number of occurrences incorrect");
        _DATA = data.begin();
        _N_OfDATA = data.ncol();
        _PLOIDY = ploidy;
        _P_OfDATA = data.nrow()/_PLOIDY;
        //_LEVELS = levels.begin();
        _N_LEVELS = n_levels.begin();
        _LEVELS_COUNT = levels_count.begin();
        _LEVELS_FREQ = levels_freq.begin();
        _PRIOR_CLASSIF = NULL;
    }
};


/**
 * @brief	Create Rcpp module
 */
RCPP_EXPOSED_CLASS(DATA)
RCPP_MODULE(MODULE_DATA)
{
  using namespace Rcpp;
  
  Rcpp::class_<DATA>("DATA")
  .constructor()
  .constructor<Rcpp::IntegerMatrix, int>()
  .constructor<Rcpp::IntegerMatrix,
          int,
          Rcpp::CharacterVector,
          Rcpp::IntegerVector,
       Rcpp::IntegerVector,
          Rcpp::DoubleVector>()
  .constructor<Rcpp::IntegerMatrix,
          int,
          Rcpp::IntegerVector,
          Rcpp::DoubleVector>()


  .field("name", &DATA::name)
  .method("howMany", &DATA::howMany)
  .method("getDATA", &DATA::getDATA)
  .method("getN", &DATA::getN)
  .method("getP", &DATA::getP)
  .method("getN_OCCURENCES", &DATA::getN_OCCURRENCES)
  .method("getLEVELS", &DATA::getLEVELS)
  .method("getN_LEVELS", &DATA::getN_LEVELS)
  .method("getLEVELS_COUNT", &DATA::getLEVELS_COUNT)
  .method("getLEVELS_FREQ", &DATA::getLEVELS_FREQ)
  .method("getPRIOR_CLASSIF", &DATA::getPRIOR_CLASSIF)
  .method("set", &DATA::set)
  .method("setDATA", &DATA::setDATA)
  //.method("setName", &DATA::setName)
  ;
}




/**
  * @brief	Set observed frequencies for a vector of prior classification
  * @warning	Values in prior_classif must be in 0:K.
  * @warning   Set _K_OfPAR_KS and _S_OfPAR_KS before.
  */
void setParObs(DATA &data, int *prior_classif, PAR_KS &parObs)
{
  int K = parObs.getK();
  Rcpp::LogicalVector S = parObs.getS();
  if(K < 2)
  {
    MyError("The number of populations is < 2");
	return;
  }

  int i, j, s, a, k, a0, N = data.getN(), P = data.getP();
  int *N_LEVELS = data.getN_LEVELS(), N_OCCURRENCES = data.getN_OCCURRENCES();
  int a1 = std::accumulate(N_LEVELS, N_LEVELS + P, 0);
  double *LEVELS_FREQ = data.getLEVELS_FREQ();
  Rcpp::NumericMatrix prob(a1, K);

  //	Compute pi_k
  Rcpp::DoubleVector pi_k(K);
  Rcpp::IntegerVector Nk(K);
  for(k = 0; k < K; k++)
  {
    s = 0;
    for(i = 0; i < N; i++)
    {
      if(prior_classif[i] == k)
        s += 1;
    }
    Nk[k] = s;
    pi_k[k] = ((double) s/N);
  }

  Rcpp::IntegerVector xN_LEVELS(P);
  Rcpp::DoubleVector xLEVELS_FREQ(a1);
  
  // Compute prob
  for(j = 0; j < P; j++)
  {
    a0 = std::accumulate(N_LEVELS, N_LEVELS + j, 0);
    xN_LEVELS[j] = N_LEVELS[j];

    for(a = 0; a < (*(N_LEVELS +j)); a++)
    {
      for(k = 0; k < K; k++)
      {
        if(S[j])
        {
          s = 0;
          for(i = 0; i < N; i++)
          {
            if(prior_classif[i] == k) s += data.howMany(a, i, j);
          }
          prob(a0 + a, k) = ((double) s/(Nk[k]*N_OCCURRENCES));
        }
        else
          prob(a0 + a, k) = LEVELS_FREQ[a0 + a];
      }
    }
  }


  for(a = 0; a < a1; a++)
      xLEVELS_FREQ[a] = LEVELS_FREQ[a];

  //	Set PAR_KS
  parObs.set(N, K, S, pi_k, prob, xN_LEVELS, xLEVELS_FREQ);
}



/**
 * @brief obsFreq
 * @param data
 * @param ploidy
 * @param levels
 * @param n_levels
 * @param levels_freq
 * @param classif
 * @param S
 * @warning K = max(classif) + 1 : populations are labeled 0:(K - 1).
 * @return
 */
// [[Rcpp::export]]
Rcpp::List obsFreq(Rcpp::IntegerMatrix data,
                   int ploidy,
                   Rcpp::CharacterVector levels,
                   Rcpp::IntegerVector n_levels,
                   Rcpp::DoubleVector levels_freq,
                   Rcpp::IntegerVector classif,
                   Rcpp::LogicalVector S)
{
    DATA xdata(data, ploidy, n_levels, levels_freq);

    int K = max(classif) + 1, N = xdata.getN();
    PAR_KS parObs(N, K, S);

    ::setParObs(xdata, classif.begin(), parObs);
    parObs.setLEVELS(levels);

    return parObs.getList();
}

/**
 * @brief logLik
 * @param data
 * @param par
 * @return
 */
double logLik(DATA &data, PAR_KS &par)
{
    int i, j, k, ploidy = data.getN_OCCURRENCES();
    int N = data.getN();
    double prob, p, lv = 0.0;

    for(i = 0; i < N; i++)
    {
      prob = 0;
      for(k = 0; k < par.getK(); k++)
      {
        p = 1;
        for(j = 0; j < data.getP(); j++)
          p *= par.getProdProb(data(i, j), ploidy, k, j);

        prob += par.getPI_k(k)*p;
      }
      if(prob == 0)
      {
        throw Rcpp::exception("Null likelihood");
      }
      lv += log(prob);
    }
    par.setN(N);

    return lv;
}


/**
 * @brief EM1 only inside Rcpp
 * @param data
 * @param par
 * @return
 */
bool EM1_Cpp(DATA &data, PAR_KS &par, double Cte)
{
  int L = data.getP(), N = data.getN(), K = 1;
  Rcpp::LogicalVector S(L);
  S.fill(false);

  int *N_LEVELS = data.getN_LEVELS(), a1 = std::accumulate(N_LEVELS, N_LEVELS + L, 0), a, j;
  double *levels_freq = data.getLEVELS_FREQ();

  Rcpp::DoubleVector pi_k(1);
  pi_k[0] = 1.0;

  Rcpp::NumericMatrix prob(a1, 1);

  for(a = 0; a < a1; a++)
      prob(a, 0) = levels_freq[a];

  Rcpp::IntegerVector vectN_LEVELS(L);

  for(j = 0; j < L; j++)
  {
      vectN_LEVELS[j] = N_LEVELS[j];
  }

  Rcpp::DoubleVector vectLevels_freq(a1);

  for(a = 0; a < a1; a++)
  {
      vectLevels_freq[a] = levels_freq[a];
  }

  par.set(N, K, S, pi_k, prob, vectN_LEVELS, vectLevels_freq);

  double lv = logLik(data, par);// FIXME faire en sorte d'utiliser la methode de par_ks

  Rcpp::NumericMatrix Tik(N, 1);
  Tik.fill(1.0);
  par.setTik(Tik);

  par.setLOG_LIK(lv);
  par.setDim();
  par.setENT();
  par.setCRITERIA(lv, Cte);

  Rcpp::IntegerVector classif(N);
  classif.fill(1);
  par.setPOST_CLASSIF(classif);//FIXME

  return true;
}

/**
 * @brief	EM1 for R
 */
// [[Rcpp::export]]
Rcpp::List EM1_Rcpp(Rcpp::IntegerMatrix tab,
               int ploidy,
               Rcpp::CharacterVector levels,
               Rcpp::IntegerVector n_levels,
               Rcpp::IntegerVector levels_count,
               Rcpp::DoubleVector levels_freq, double Cte)
{
  DATA data(tab, ploidy, levels, n_levels, levels_count, levels_freq);
  PAR_KS par;
  EM1_Cpp(data, par, Cte);
  return par.getList();
}

/**
 * @brief Expectation_Cpp
 * @param data
 * @param par
 * @param Tik
 * @return
 */
bool Expectation_Cpp (DATA & data,
                      PAR_KS & par,
                      double *Tik)
{
  int N=data.getN(), L=data.getP(), K=par.getK();

  std::vector<double> som_Tik_k(N) ;
  int i, k, j;
  double s, p ;

  for(i=0; i<N ; i++)
  {
    for(k=0; k<K; k++)
    {
      p = 1.0 ;
      for(j = 0; j < L; j++)
      {
        p *= par.getProdProb(data(i, j), data.getN_OCCURRENCES(), k, j) ;
      }
      Tik[i*K+k] = p*par.getPI_k(k);
    }

    // renormalisation
    s = 0.0 ;
    for(k=0 ; k<K ;k++)
      s += Tik[i*K+k] ;
    som_Tik_k[i] = s;
    if(som_Tik_k[i] <= 0.0)
    {
      throw Rcpp::exception("Invalide value");
      return false;
    }

    for(k=0 ;k<K ; k++)
      Tik[i*K+k] = ((double)Tik[i*K+k])/som_Tik_k[i] ;
  }

  return true;
}

/**
 * @brief Expectation
 * @param tab
 * @param ploidy
 * @param levels
 * @param n_levels
 * @param levels_count
 * @param levels_freq
 * @param par_ks
 * @return
 */
Rcpp::NumericMatrix Expectation(Rcpp::IntegerMatrix tab,
                                int ploidy,
                                Rcpp::CharacterVector levels,
                                Rcpp::IntegerVector n_levels,
                                Rcpp::IntegerVector levels_count,
                                Rcpp::DoubleVector levels_freq,
                                Rcpp::List listParKS)
{
  DATA data(tab, ploidy, levels, n_levels, levels_count, levels_freq);

  PAR_KS par_ks(listParKS);

  int N = data.getN(), K = par_ks.getK(), i, k;
  std::vector<double> Tik(N*K);

  Expectation_Cpp(data, par_ks, Tik.data());

  Rcpp::NumericMatrix matTik(N, K);
  for(i= 0; i < N; i++)
  {
      for(k = 0; k < K; k++)
          matTik(i, k) = Tik[i*K+k];
  }

  return matTik;
}

/**
 * @brief Maximisation_Cpp
 * @param data
 * @param par
 * @param Tik
 * @return
 */
bool Maximisation_Cpp(DATA &data, PAR_KS &par, double * Tik)
{
  int N=data.getN(), L=data.getP(), K=par.getK();
  int *N_LEVELS = data.getN_LEVELS();

  if(K==1) 
  {
      MyError("The number of population is not > 1");
      throw Rcpp::exception("The number of population is not > 1");
      return false;
  }

  int i, k, j, a, ploidy = data.getN_OCCURRENCES();
  double s;
  std::vector<double> s1(K) ;

  // calcul des proportions
  for(k=0; k<K-1; k++)
  {
    s=0.0;
    for(i = 0; i < N; i++)
      s+= Tik[i*K+k] ; // Effectif attendu de la classe k

    s1[k] = s;
    par.setPI_k(((double) s)/N, k);
  }
  s = 1.0;
  for(k = 0; k < K-1; k++)
    s -= par.getPI_k(k);
  if(s < 0) s = 0.0;
  par.setPI_k(s, K-1);

  // Effectif attendu de la derniére classe s1[K-1] = N*par.getPI_k(K-1) = old
  s1[K-1] = N*s;

  double xThreshold = _PutTHRESHOLD? _THRESHOLD : 0.0;

  for(j=0; j<L; j++)
    if(par.getS_j(j))
    {
      for(k=0; k<K; k++)
      {
        for(a=0 ; a<N_LEVELS[j]-1; a++)
        {
          s=0.0;
          for(i = 0; i < N; i++)
            s += Tik[i*K + k] * data.howMany(a, i, j);
          s += ((double) xThreshold)/N_LEVELS[j] ; // _THRESHOLD for alleleic frequencies
          par.setPROBkja(((double) s)/(ploidy*s1[k] + xThreshold), k, j, a);
        }
        //  Assure que _PROB_ALL soit constitue de simplexes
        s = 1.0;
        for(a = 0; a < N_LEVELS[j]-1; a++)
          s -= par.getPROBkja(k, j, a);
        par.setPROBkja(s>0.0? s: 0.0, k, j, N_LEVELS[j]-1);
      }
    }

  return true;
}

/**
 * @brief mapClassification_Cpp
 * @param N
 * @param K
 * @param Tik
 * @param vect
 */
void mapClassification_Cpp(int N, int K, double *Tik, int *vect)
{
    for(int i = 0; i < N; i++)
      vect[i]=whichMax(Tik+(i*K), Tik+(i*K+K));
}

/**
 * @brief mapClassification_Cpp2
 * @param N
 * @param K
 * @param Tik
 * @param vect
 */
Rcpp::IntegerVector mapClassification_Cpp2(int N, int K, double *Tik)
{
    Rcpp::IntegerVector vect(N);
    for(int i = 0; i < N; i++)
      vect[i]=whichMax(Tik+(i*K), Tik+(i*K+K));
    return vect;
}

/**
 * @brief mapClassification_Rcpp
 * @param Tik
 * @return
 */
// [[Rcpp::export]]
Rcpp::IntegerVector mapClassification_Rcpp(Rcpp::NumericMatrix Tik)
{
    int k, K = Tik.ncol(), N = Tik.nrow();
    std::vector<double> prob(Tik.ncol());
    Rcpp::IntegerVector classif(N, K);

    for(int i = 0; i < N; i++)
    {
        for(k = 0; k < K; k++)
            prob[k] = Tik(i, k);
        classif[i] = whichMax(prob.data(), prob.data() +K);
    }

    return classif;
}



/**
 * @brief smallEM_Cpp
 * @param data
 * @param par
 * @param typeEM = 0 for classical EM, 1 for SEM and 2 for CEM.
 * @return
 */
bool smallEM_Cpp(DATA &data,
                  PAR_KS &par)
{
  int N=data.getN(), L=data.getP(), K = par.getK(), j, a;
  int *N_LEVELS = data.getN_LEVELS(), a1 = std::accumulate(N_LEVELS, N_LEVELS + L, 0);
  Rcpp::LogicalVector S = par.getS();
  double *LEVELS_FREQ = data.getLEVELS_FREQ();

  Rcpp::IntegerVector vectN_LEVELS(L);

  for(j = 0; j < L; j++)
  {
      vectN_LEVELS[j] = N_LEVELS[j];
  }

  Rcpp::DoubleVector vectLEVELS_FREQ(a1);
  for(a = 0; a < a1; a++)
      vectLEVELS_FREQ[a] = LEVELS_FREQ[a];

  if(K==1)
  {
    Rcpp::Rcout << "> K = 1\n> Empirical statistics are used\n";
    EM1_Cpp(data, par, 1.0);
    return true;
  }
  // Initialiser _Tik
  std::vector<double> Tik(K*N) ;
  double lv;
  int i;

  par.randomInitialise(N, K, S, vectN_LEVELS, vectLEVELS_FREQ);

  try
  {
    lv = ::logLik(data, par);
    par.setLOG_LIK(lv);
  }
  catch(std::string s)
  {
    MyError("Computing log-likelihood");
    par.print();
    return false;
  }

  PAR_KS par_0;

  double lv0;
  std::vector<int> vect(N);

  switch(_TYPE_SMALL_EM)
  {
    case 0:
      for(i = 0; i < _NBER_SMALL_EM; i++)
      {
        par_0.randomInitialise(N, K, S, vectN_LEVELS, vectLEVELS_FREQ);//FIXME
        for(j = 0 ; j<_NBER_ITER_EM; j++)
        {
          if(!Expectation_Cpp(data, par_0, Tik.data()))
            return false;

          if(!Maximisation_Cpp(data, par_0, Tik.data()))
            return false;
        }

        lv0 = ::logLik(data, par_0);
        par_0.setLOG_LIK(lv0);

        if(par_0.getLOG_LIK() > par.getLOG_LIK())
          par = par_0 ;
      }
      break;

    case 1: //FIXME
      for(i = 0; i < _NBER_SMALL_EM; i++)
      {
        par_0.randomInitialise(N, K, S, vectN_LEVELS, vectLEVELS_FREQ);//FIXME
        for(j = 0 ; j < _NBER_ITER_EM; j++)
        {
          if(!Expectation_Cpp(data, par_0, Tik.data()))
            return false;

          stochastique(N, K, Tik.data(), vect.data());

          setParObs(data, vect.data(), par_0);
        }

        lv0 = ::logLik(data, par_0);
        par_0.setLOG_LIK(lv0);

        if(par_0.getLOG_LIK() > par.getLOG_LIK())
          par = par_0 ;
      }
      break;

    case 2:
      for(i = 0; i < _NBER_SMALL_EM; i++)
      {
        par_0.randomInitialise(N, K, S, vectN_LEVELS, vectLEVELS_FREQ);
        for(j = 0 ; j<_NBER_ITER_EM; j++)
        {
          if(!Expectation_Cpp(data, par_0, Tik.data()))
            return false;

          mapClassification_Cpp(N, K, Tik.data(), vect.data());

          setParObs(data, vect.data(), par_0);
        }

        lv0 = ::logLik(data, par_0);
        par_0.setLOG_LIK(lv0);

        if(par_0.getLOG_LIK() > par.getLOG_LIK())
          par = par_0 ;
      }
      break;

    default:
      throw Rcpp::exception("Incorrect choice of the type of EM");
      return false;
  }

  Expectation_Cpp(data, par, Tik.data());
  par.setTik2(Tik.data());
  par.setENT();
  par.setDim();
  par.setCRITERIA(par.getLOG_LIK(), 1.0);

  return true;
}

/**
 * @brief smallEM_Rcpp
 * @param tab
 * @param ploidy
 * @param levels
 * @param n_levels
 * @param levels_count
 * @param levels_freq
 * @return
 */
// [[Rcpp::export]]
Rcpp::List smallEM_Rcpp(Rcpp::IntegerMatrix tab,
                    int ploidy,
                    Rcpp::CharacterVector levels,
                    Rcpp::IntegerVector n_levels,
                    Rcpp::IntegerVector levels_count,
                    Rcpp::DoubleVector levels_freq,
                   int K,
                   Rcpp::LogicalVector S)
{
    DATA data(tab, ploidy, levels, n_levels, levels_count, levels_freq);
    PAR_KS par(data.getN(), K, S, n_levels, levels_freq);

    if(!smallEM_Cpp(data, par))
    {
        MyError("Running small EM");
        throw Rcpp::exception("Running small EM");
    }

    return par.getList();
}


/**
 * @brief EM_Cpp
 * @param data
 * @param par
 * @param Cte
 * @param typeSmallEM
 * @param typeEM
 * @return
 */
bool EM_Cpp(DATA &data,
        PAR_KS &par,
        double Cte)
{
    int K = par.getK();

    if(K == 1)
    {
        return EM1_Cpp(data, par, Cte);
    }

    int N=data.getN(), L = data.getP(), j;
    int *N_LEVELS = data.getN_LEVELS();
    int a1 = std::accumulate(N_LEVELS, N_LEVELS + L, 0), a;
    double *levels_freq = data.getLEVELS_FREQ();
    Rcpp::IntegerVector vectN_LEVELS(L);
    Rcpp::DoubleVector vectLEVELS_FREQ(a1);
    Rcpp::LogicalVector S = par.getS();

    for(j = 0; j < L; j++)
    {
        vectN_LEVELS[j] = N_LEVELS[j];
    }

    for(a = 0; a < a1; a++)
        vectLEVELS_FREQ[a] = levels_freq[a];

    par.randomInitialise(N, K, S, vectN_LEVELS, vectLEVELS_FREQ);


    // Perform small EM
    Rprintf("\n ... Running %d small EM with %d iterations each...", _NBER_SMALL_EM, _NBER_ITER_EM);

    if(!smallEM_Cpp(data, par))
    {
        MyError("runing small EM");
        throw Rcpp::exception("");
        return false;
    }

    Rprintf("\n ... Runing a maximum of %d long run of EM...\n", _NBER_ITER_LONG_EM);

    // Initialiser _Tik
    std::vector<double> Tik(K*N) ;
    double D_lv, lv_1, lv_2;
    D_lv=10;
    lv_2 = par.getLOG_LIK();

    std::vector<int> vectClassif(N);

    switch(_TYPE_EM)
    {
    case 0:
        for(j = 0; D_lv>_EPSI && j < _NBER_ITER_LONG_EM; j++)
        {
            try
            {
              Expectation_Cpp(data, par, Tik.data());
              Maximisation_Cpp(data, par, Tik.data());
              lv_1=::logLik(data, par);
              par.setLOG_LIK(lv_1);
            }
            catch(std::string s)
            {
                MyError("in EM algorithm");
                par.print();
                return false;
            }
            D_lv = ((double) (lv_2 - lv_1)/lv_1) ; //erreur relative et non absolue
            D_lv = (D_lv < 0)? - D_lv : D_lv;
            lv_2 = lv_1;
        }
        break;

    case 1:
        for(j = 0; D_lv > _EPSI && j < _NBER_ITER_LONG_EM; j++)
        {
            try
            {
              Expectation_Cpp(data, par, Tik.data());

              stochastique(N, K, Tik.data(), vectClassif.data());

              ::setParObs(data, vectClassif.data(), par);

              lv_1 = ::logLik(data, par);
              par.setLOG_LIK(lv_1);
            }
            catch(std::string s)
            {
                MyError(" > runing EM ");
                par.print();
                return false;
            }
            D_lv = ((double) (lv_2 - lv_1)/lv_1) ; //erreur relative et non absolue
            D_lv = (D_lv < 0)? - D_lv : D_lv;
            lv_2 = lv_1;
        }
        break;

    case 2:
        for(j = 0; D_lv > _EPSI && j < _NBER_ITER_LONG_EM; j++)
        {
            try
            {
              Expectation_Cpp(data, par, Tik.data());

              mapClassification_Cpp(N, K, Tik.data(), vectClassif.data());

              ::setParObs(data, vectClassif.data(), par);

              lv_1=::logLik(data, par);

              par.setLOG_LIK(lv_1);
            }
            catch(std::string s)
            {
                MyError("runing EM ");
                par.print();
                return false;
            }
            D_lv = ((double) (lv_2 - lv_1)/lv_1) ; //erreur relative et non absolue
            D_lv = (D_lv < 0)? - D_lv : D_lv;
            lv_2 = lv_1;
        }
        break;

    default:
        MyError("Incorrect Argument typeEM");
        return true;
    }

    Rprintf("> Number of iterations = %d\n", j);

    Expectation_Cpp(data, par, Tik.data());
    par.setTik2(Tik.data());
    par.setENT();
    par.setDim();
    par.setCRITERIA(par.getLOG_LIK(), 1.0);

    Rcpp::IntegerVector classif(N);
    classif = mapClassification_Cpp2(N, K, Tik.data());
    par.setPOST_CLASSIF(classif);

    return true;
}


/**
 * @brief EM_Rcpp
 * @param tab
 * @param ploidy
 * @param levels
 * @param n_levels
 * @param levels_count
 * @param levels_freq
 * @param K
 * @param S
 * @param typeSmallEM
 * @param typeEM
 * @return
 */
// [[Rcpp::export]]
Rcpp::List EM_Rcpp(Rcpp::IntegerMatrix tab,
                    int ploidy,
                    Rcpp::CharacterVector levels,
                    Rcpp::IntegerVector n_levels,
                    Rcpp::IntegerVector levels_count,
                    Rcpp::DoubleVector levels_freq,
                   int K,
                   Rcpp::LogicalVector S,
				   double Cte = 1.0)
{
    DATA data(tab, ploidy, levels, n_levels, levels_count, levels_freq);
    PAR_KS par(data.getN(), K, S, n_levels, levels_freq);
    if(K == 1)
    {
        EM1_Cpp(data, par, Cte);
        return par.getList();
    }

    EM_Cpp(data, par, Cte);
	par.setLEVELS(levels);

    return par.getList();
}

/**
 * @brief readModelFromString_Rcpp Read a model from a string
 * @param mod string
 * @return Rcpp::List(P, K, S, logLik, dim, entropy)
 */
// [[Rcpp::export]]
Rcpp::List readModelFromString_Rcpp(std::string mod)
{
  int n = howmanyWords(mod), i, l;

  std::istringstream tampon (mod) ;

  int N;
  if(!(tampon >> N))
  {
      MyError("Incorrect value of N!");
      return false;
  }

  int L;
  if(!(tampon >> L))
  {
      MyError("Incorrect value of L!");
      return false;
  }

  if(n < L + NBER_MIN_COL_EXPLORED_MODELS_FILE)
  {
      MyError("Incorrect number of columns!");
      return false;
  }

  int K;
  if(!(tampon >> K))
  {
      MyError("Incorrect value of K !");
      return false;
  }

  Rcpp::LogicalVector S(L);

  for(l=0; l<L; l++)
  {
    if((tampon >> i) && (i == 0 || i == 1))
    {
      if(i == 1) S[l]=true;
      else S[l]=false;
    }
    else
    {
      MyError("Incorrect value in S");
      return false;
    }
  }

  double logLik;
  if(!(tampon >> logLik))
  {
    MyError("Incorrect value of log-likelihood!");
    return false;
  }
  
  // FIXME Problem avec cette lecture
  int d_KS;
  if(!(tampon >> d_KS))
  {
    MyError("Incorrect dimension of model");
    return false;
  }

  
  double entropy;
  if(!(tampon >> entropy))
  {
    MyError("Incorrect value of entropy!");
    return false;
  }

  Rcpp::List outList; // Create names TODO
  
  outList["N"] = N;
  outList["P"] = L;
  outList["K"] = K;
  outList["S"] = S;
  outList["logLik"] = logLik;
  outList["dim"] = d_KS;
  outList["entropy"] = entropy;

  return outList;
}

/**
 * @brief isInFile_Rcpp
 * @param parKS
 * @param fichier
 * @return
 */
// [[Rcpp::export]]
Rcpp::List isInFile_Rcpp(int K, Rcpp::LogicalVector S, std::string fichier, bool header)
{
  // verifie si un modele a ete enregistre dans un fichier
  // De meme que dans l'operateur == que cette fonction utilise,
  // les modeles avec K==1 sont identifies a ceux avec S vide
  int i = 0, num_ligne = 0, P = S.length(), P2, K2;
  std::string modele_ligne;
  std::vector<int> vect(P);
  int v;
  Rcpp::List outList;

	std::ifstream  file_in(fichier.c_str(), std::ios_base::in);
  
  if(file_in.bad())
  {
		MyError("opening file");
		outList["TrueFalse"] = false;
		return outList;
  }
  
//   std::ifstream file_in;
//   try 
//   {
// 		file_in.open(fichier.c_str(), std::ios_base::in);
// 		file_in.exceptions(file_in.eofbit |file_in.failbit);
//   } catch (const std::ios_base::failure& e)
//   {
// 		Rcpp::Rcout << "Caught an ios_base::failure.\n"
//                   << "Explanatory string: " << e.what() << '\n';
// 		outList["TrueFalse"] = false;
// 		return outList;
//   }

  if(header) nextLine(file_in, modele_ligne);

  while(nextLine(file_in, modele_ligne))
  {
    Rcpp::List parList = readModelFromString_Rcpp(modele_ligne);
    if(parList.length() == 0)
    {
			         file_in.close();
			MyError("Incorrect model at some line");
			outList["TrueFalse"] = false;
			return outList;
    }


    int N = Rcpp::as<int>(parList["N"]);
    P2 = Rcpp::as<int>(parList["P"]);
    K2 = Rcpp::as<int>(parList["K"]);
    Rcpp::LogicalVector S2 = Rcpp::as<Rcpp::LogicalVector>(parList["S"]);
    
	
    for(i = 0; i < P; i++)
        vect[i] = (S[i] == S2[i])? 1: 0;
	
    v = 1;
		for(i = 0; i < P; i++)
			v *= vect[i];


    if((P == P2) && (K == K2) && (v == 1))
    {
			         file_in.close() ;
			outList["TrueFalse"] = true;
			outList["line"] = num_ligne;
			outList["N"] = N;
      outList["logLik"] = Rcpp::as<double>(parList["logLik"]);
      outList["dim"] = Rcpp::as<int>(parList["dim"]);
      outList["entropy"] = Rcpp::as<double>(parList["entropy"]);
      return outList;
    }
    num_ligne++;
  }

  file_in.close() ;
	outList["TrueFalse"] = false;
  
  return outList;
}


/**
 */
//	[[Rcpp::export()]]
Rcpp::DoubleVector computeCriteria_Rcpp(double lv, int dim, int N, double entropy = 0, double Cte = 1.0)
{
	double xdim = (double) dim;
	Rcpp::DoubleVector criteria = Rcpp::DoubleVector::create(_["BIC"] = -lv + (xdim*log((double) N))/2,
                                           _["AIC"] = -lv + xdim,
                                        _["ICL"] = -lv + (xdim*log((double) N))/2 + entropy,
                                           _["CteDim"] = -lv + Cte*xdim);
	return criteria;
}
/**
 * @brief	Find indexes of the best models by criterion
 */
Rcpp::IntegerVector findBestModels(Rcpp::NumericMatrix criteria)
{
  Rcpp::IntegerVector bestModels(NBER_CRITERIA);
  int n = criteria.nrow(), j;
  
  for(j = 0; j < NBER_CRITERIA; j++)
	bestModels[j] = whichMin(criteria.begin() + j*n, n);
  
  return bestModels;
}

/**
 * @brief	Compute a matrix of criteria from a file containing explored models.
 * @return	A list of computed criteria matrix and a vector of indexes of best models fro different criteria
 */
//	[[Rcpp::export()]]
Rcpp::List computeCriteriaFromFile_Rcpp(std::string xfile,
                   double Cte,
                   bool header,
                   Rcpp::IntegerVector indexes = Rcpp::IntegerVector::create())
{
  int nLines;
  int i;
  if(indexes.length() != 0)
  {
		nLines = indexes.length();
  }
  else
  {
		nLines = header? nberOfLines(xfile) - 1 : nberOfLines(xfile);
		indexes = Rcpp::IntegerVector(nLines);
		for(i = 0; i < nLines; i++)
      indexes[i] = i;
  }
 
  std::ifstream fIn(xfile.c_str(), std::ios_base::in);
  
  Rcpp::NumericMatrix criteria(nLines, NBER_CRITERIA);
  Rcpp::List outPutList;
	
	if(fIn.bad())
	{
		ErrorOpeningFile();
		return outPutList;
	}
  
/*  try 
  {
		fIn.exceptions(fIn.failbit);
  } catch (const std::ios_base::failure& e)
  {
		Rcpp::Rcout << "Caught an ios_base::failure.\n"
                  << "Explanatory string: " << e.what() << '\n';
		return outPutList;
  }
 */ 
  if(nLines <= 0)
  {
		MyError("opening xfile");
		return outPutList;
  }
  
  std::string modelString;
  int dim, j, N;
  double logLik, entropy;
  
  if(header)
	nextLine(fIn, modelString);
  
  i = 0;
  int orderLine = 0;
  while(nextLine(fIn, modelString))
  {
    if(indexes[i] < 0)
          throw Rcpp::exception("Negative integer not allowed");

    if(orderLine == indexes[i])
	{
      Rcpp::List parList = readModelFromString_Rcpp(modelString);
      N = Rcpp::as<int>(parList["N"]);
      logLik = Rcpp::as<double>(parList["logLik"]);
      dim = Rcpp::as<int>(parList["dim"]);
      entropy = Rcpp::as<double>(parList["entropy"]);
	
      Rcpp::DoubleVector vectCriteria = computeCriteria_Rcpp(logLik, dim, N, entropy, Cte);
	  for(j = 0; j < NBER_CRITERIA; j++)
		criteria(i, j) = vectCriteria[j];
	
	  i++;
	}
	orderLine++;
  }
  
  
  fIn.close();
  
  // Find the best model
  Rcpp::IntegerVector bestModels = findBestModels(criteria);
  outPutList["criteria"] = criteria;
  outPutList["bestModelsIndexes"] = bestModels;
  
  
  return outPutList;
}
  
/**
 * @brief	Write a model in a file
 * @param	x ) list(P, K, S, logLik, dim, entropy)
 * @return 	void
 */
//FIXME Delete this file and replace by writeModelInFile_Rcpp
//	[[Rcpp::export()]]
void writeParInFile_Rcpp(Rcpp::List x, std::string xfile)
{
  std::ofstream out ;
  out.open(xfile.c_str(), std::ios_base::app);
  out.precision(PRECISION);
  int P = Rcpp::as<int>(x["P"]);
  Rcpp::LogicalVector S = x["S"];
  out << Rcpp::as<int>(x["N"]) << " " << P << " " << Rcpp::as<int>(x["K"]) << " " ;
  for(int i = 0; i < P; i++)
	if(S[i]) out << 1 << " ";
	else out << 0 << " ";
  out  << std::fixed << Rcpp::as<double>(x["logLik"]) << " " << Rcpp::as<int>(x["dim"]) << " " << Rcpp::as<double>(x["entropy"]) << "\n";
  
  out.close();
}

/**
 * @brief	Write a model in a file (Explored models file)
 * @warning	
 */
//	[[Rcpp::export()]]
void writeModelInFile_Rcpp(Rcpp::List x, std::string xfile)
{
  std::ofstream out ;
  out.open(xfile.c_str(), std::ios_base::app);
  out.precision(PRECISION);
  int P = Rcpp::as<int>(x["P"]);
  Rcpp::LogicalVector S = x["S"];
  out << Rcpp::as<int>(x["N"]) << " " << P << " " << Rcpp::as<int>(x["K"]) << " " ;
  for(int i = 0; i < P; i++)
	if(S[i]) out << 1 << " ";
	else out << 0 << " ";
  out  << std::fixed << Rcpp::as<double>(x["logLik"]) << " " << Rcpp::as<int>(x["dim"]) << " " << Rcpp::as<double>(x["entropy"]) << "\n";
  
  out.close();
}



/**
 * @brief writeCriteriaInFile_Rcpp
 * @param criteria
 * @param xfile
 */
//	[[Rcpp::export()]]
void writeCriteriaInFile_Rcpp(Rcpp::DoubleVector criteria, std::string xfile)
{
  int n = criteria.length();
  std::ofstream out;
  out.open(xfile.c_str(), std::ios_base::app);
  out.precision(PRECISION);
  for(int i = 0; i < n; i++)
	out << std::fixed << criteria[i] << " ";
  out << "\n";
  
  out.close();
}



/**
 * @brief readModelAt
 * @param xfile
 * @param orderLine
 * @param header
 * @return
 */
//	[[Rcpp::export()]]
Rcpp::List readModelAt_Rcpp(std::string xfile, int n, bool header)
{
  Rcpp::List outPutList;
  int nLines = header? nberOfLines(xfile) - 1 : nberOfLines(xfile);
  if(n >= nLines)
  {
    MyError("line out of range");
		return outPutList;
  }
  
  std::ifstream fIn(xfile.c_str(), std::ios_base::in);
  
  if(fIn.bad())
  {
    ErrorOpeningFile();
		return outPutList;
  }
/*  
  try 
  {
		fIn.open(xfile.c_str(), std::ios::in);
		fIn.exceptions(fIn.failbit);
  } catch (const std::ios_base::failure& e)
  {
		Rcpp::Rcout << "Caught an ios_base::failure.\n"
                  << "Explanatory string: " << e.what() << '\n';
    return outPutList;
  }*/
  
  std::string modelString;
  if(header)
	nextLine(fIn, modelString);

  for(int i = 0; i <= n; i++)
    if(!nextLine(fIn, modelString))
      {
        MyError("cannot read line some line ");
        return outPutList;
      }
  fIn.close();
  return readModelFromString_Rcpp(modelString);
}


/**
 * @brief selectedDimnensionRcpp
 * @param xfileExploredModels
 * @param constantGrid
 * @param vectLogLik
 * @param vectDim
 * @param header
 */
// [[Rcpp::export()]]
bool selectDimFromFile_Rcpp(std::string xfileExploredModels,
                                           Rcpp::DoubleVector constantGrid,
                                           Rcpp::DoubleVector vectLogLik,
                                           Rcpp::IntegerVector vectDim,
                                           bool header)
{
  std::ifstream fIn(xfileExploredModels.c_str(), std::ios_base::in);
  
  if(fIn.bad())
    {
      ErrorOpeningFile();
      return false;
    }
  
// 	std::ifstream fIn;
//   try 
//   {
// 		fIn.open(xfileExploredModels.c_str(), std::ios::in);
// 		fIn.exceptions(fIn.failbit);
//   } catch (const std::ios_base::failure& e)
//   {
// 		Rcpp::Rcout << "Caught an ios_base::failure.\n"
//                   << "Explanatory string: " << e.what() << '\n';
//     return false;
//   }

  std::string modelString;
  double logLik;
  int N, P, K, i, j, dim, nGrid = constantGrid.length();
  Rcpp::DoubleVector vectCrit(nGrid);
  vectCrit.fill(std::numeric_limits<double>::infinity());
  double crit_i;

  if(header) nextLine(fIn, modelString);

  while(nextLine(fIn, modelString))
    {
      std::istringstream tampon(modelString);

      if(!(tampon >> N))
        {
          fIn.close();
          MyError("Incorrect value in first column");
          return false;
        }

      if(!(tampon >> P))
        {
          fIn.close();
          MyError("Incorrect value in second column");
          return false;
        }

      if(howmanyWords(modelString) < P + NBER_MIN_COL_EXPLORED_MODELS_FILE)
        {
          fIn.close();
          MyError("Incorrect number of column in file");
          return false;
        }


      if(!(tampon >> K))
        {
          fIn.close();
          MyError("Incorrect number of population");
          return false;
        }

      for(j = 0; j < P; j++)
        {
          if(!(tampon >> i))
            {
              fIn.close();
              MyError("Incorrect value for S");
              return false;
            }
        }

      if(!(tampon >> logLik))
        {
          fIn.close();
          MyError("Incorrect value of logLik");
          return false;
        }

      if(!(tampon >> dim))
        {
          fIn.close();
          MyError("Incorrect value of dimension");
          return false;
        }

      for(i = 0; i < nGrid; i++)
        {
          crit_i = - logLik + constantGrid[i]*dim;
          if(crit_i < vectCrit[i])
            {
              vectCrit[i] = crit_i;
              vectLogLik[i] = logLik;
              vectDim[i] = dim;
            }
        }
    }

  fIn.close();

  return true;
}

/**
 * @brief selectDimFromData_Rcpp
 * @param xlogLik
 * @param xdim
 * @param xconstantGrid
 * @param outLogLik
 * @param outDim
 * @return
 */
// [[Rcpp::export()]]
bool selectDimFromData_Rcpp(Rcpp::DoubleVector xlogLik,
                            Rcpp::IntegerVector xdim,
                            Rcpp::DoubleVector xconstantGrid,
                            Rcpp::DoubleVector outLogLik,
                            Rcpp::IntegerVector outDim)
{
    int n_in = xlogLik.length(), n_out = xconstantGrid.length();
    if(n_in != xdim.length() ||
            outLogLik.length() != n_out ||
            outDim.length() != n_out)
    {
        MyError("Dimensions not compatibles");
        return false;
    }

    Rcpp::DoubleVector vectCrit(n_out);
    vectCrit.fill(std::numeric_limits<double>::infinity());

    int i, j;
    double crit_i;
    for(j = 0; j < n_out; j++)
    {
        for(i = 0; i < n_in; i++)
          {
            crit_i = - xlogLik[i] + xconstantGrid[j]*xdim[i];
            if(crit_i < vectCrit[j])
              {
                vectCrit[j] = crit_i;
                outLogLik[j] = xlogLik[i];
                outDim[j] = xdim[i];
              }
          }
    }

    return true;
}

/**
 * @brief detectJump
 * @param vectDim
 * @param pas
 * @param BeginEnd
 */
bool dimJumpRcpp_old(Rcpp::IntegerVector vectDim,
                int pas,
                Rcpp::IntegerVector BeginEnd1,
                Rcpp::IntegerVector BeginEnd2)
{
  int *dim = vectDim.begin();
  int *ptBeginEnd1 = BeginEnd1.begin();
  int *ptBeginEnd2 = BeginEnd2.begin();
  int lengthVectDim = vectDim.length();
  if(lengthVectDim < 1 || pas < 1)
  {
      MyError("Incorrect argument");
      return false;
  }
  int i, j, k;
  ptBeginEnd1[1] = 1;
  ptBeginEnd2[1] = 0;
  ptBeginEnd1[0] = 0;
  ptBeginEnd2[0] = 0;
  double saut, saut_1 = 0.0;

  for(i = 1; i < lengthVectDim; i++)
  {
    if(i < pas) j = 0;
    else j = i - pas;
    saut = std::abs(dim[j] - dim[i]);

    if(saut > saut_1)
    {
      // trouver le debut du saut
      for(k = j+1; k < i; k++)
      {
        if(saut == std::abs(dim[k] - dim[i]))
          j = k;
        else break;
      }

      // enregistrer le precedent comme 2e meilleur saut
      if(ptBeginEnd1[1] >= ptBeginEnd2[1] + pas && saut_1 > 0.0)
      {
            ptBeginEnd2[1] = ptBeginEnd1[1];
            ptBeginEnd2[0] = ptBeginEnd1[0];
      }

      // enregistrer le saut courrant comme meilleur saut
      saut_1 = saut;
      ptBeginEnd1[1] = i;
      ptBeginEnd1[0] = j;
    }
  }
  return true;
}

/**
 * @brief detectJumpRcpp
 * @param vectDim
 * @param pas
 * @param BeginEnd1
 * @param BeginEnd2
 * @return
 */
// [[Rcpp::export()]]
bool dimJump_Rcpp(Rcpp::IntegerVector vectDim,
                int pas,
                Rcpp::IntegerVector BeginEnd1,
                Rcpp::IntegerVector BeginEnd2)
{
  int *dim = vectDim.begin();
  int *ptBeginEnd1 = BeginEnd1.begin();
  int *ptBeginEnd2 = BeginEnd2.begin();
  int lengthVectDim = vectDim.length();
  if(lengthVectDim < 1 || pas < 1)
  {
      MyError("Incorrect argument");
      return false;
  }
  int i, j, k;
  ptBeginEnd1[1] = 1;
  ptBeginEnd2[1] = 0;
  ptBeginEnd1[0] = 0;
  ptBeginEnd2[0] = 0;
  double saut, saut1 = 0.0;
  int j0;

  for(i = 1; i < lengthVectDim; i++)
  {
    if(i < pas) j = 0;
    else j = i - pas;
    saut = std::abs(dim[j] - dim[i]);

    j0 = j;

    if(saut > saut1)
    {
      // trouver le debut du saut
      for(k = j0 + 1; k < i; k++)
      {
        if(saut == std::abs(dim[k] - dim[i]))
          j = k;
        else break;
      }

      // enregistrer le saut courrant comme meilleur saut
      saut1 = saut;
      ptBeginEnd1[0] = j;
      ptBeginEnd1[1] = i;
    }
  }

  saut1 = 0.0;
  for(i = 1; i < lengthVectDim; i++)
  {
    if(i < pas) j = 0;
    else j = i - pas;
    saut = std::abs(dim[j] - dim[i]);

    j0 = j;

    if(saut > saut1)
    {
      // trouver le debut du saut
      for(k = j0 + 1; k < i; k++)
      {
        if(saut == std::abs(dim[k] - dim[i]))
          j = k;
        else break;
      }

      // enregistrer le saut courrant comme meilleur saut
      //if(j > ptBeginEnd1[1])// old
      if((i < ptBeginEnd1[0]) | (j > ptBeginEnd1[1]))
        {
          saut1 = saut;
          ptBeginEnd2[0] = j;
          ptBeginEnd2[1] = i;
        }
    }
  }

  return true;
}

/**
 * @brief selectModelRcpp
 * @param xfile
 * @param N
 * @param cte
 * @return
 */
// FIXME
//  [[Rcpp::export()]]
void selectModelFromFile_Rcpp(std::string xfileExploredModels,
                              Rcpp::IntegerVector vectN,
                                      Rcpp::IntegerVector vectK,
                                      Rcpp::IntegerMatrix matS,
                                      Rcpp::DoubleVector vectLogLik,
                                      Rcpp::IntegerVector vectDim,
                                      Rcpp::DoubleVector vectEntropy,
                                      Rcpp::DoubleVector vectCriteria,
                                      double cte,
                                bool header,
                                Rcpp::IntegerVector lines)
{
  //vectN.length() = 1 for N
  // Compatibility of arguments
  if(vectK.length() != NBER_CRITERIA)
    {
      throw Rcpp::exception("Incorrect dimension of argument 'vectK'");
      return;
    }

  if((matS.nrow() != NBER_CRITERIA))
    {
      throw Rcpp::exception("Incorrect number of lines in 'matS");
      return;
    }

  if(vectLogLik.length() != NBER_CRITERIA)
    {
      throw Rcpp::exception("Incorrect dimension of argument 'vectLogLik'");
      return;
    }

  if(vectDim.length() != NBER_CRITERIA)
    {
      throw Rcpp::exception("Incorrect dimension of argument 'vectDim'");
      return;
    }

  if(vectEntropy.length() != NBER_CRITERIA)
    {
      throw Rcpp::exception("Incorrect dimension of argument 'vectEntropy'");
      return;
    }

  if(vectCriteria.length() != NBER_CRITERIA)
    {
      throw Rcpp::exception("Incorrect dimension of argument 'vectCriteria'");
      return;
    }

  std::ifstream fIn(xfileExploredModels.c_str(), std::ios_base::in);
  
  if(fIn.bad())
    {
      throw Rcpp::exception("Can not open file");
    }
    
// 	std::ifstream fIn;
//   try 
//   {
// 		fIn.open(xfileExploredModels.c_str(), std::ios::in);
// 		fIn.exceptions(fIn.failbit);
//   } catch (const std::ios_base::failure& e)
//   {
// 		Rcpp::Rcout << "Caught an ios_base::failure.\n"
//                   << "Explanatory string: " << e.what() << '\n';
// 		return;
//   }

  std::string modelString;
  double logLik, entropy;
  int xP, xK, i, j, dim, s, k;

  // Pointers
  int *ptVectK = vectK.begin(), *ptMatS = matS.begin();
  int *ptVectDim = vectDim.begin(), *ptVectN = vectN.begin();
  double *ptVectLogLik = vectLogLik.begin(), *ptVectEntropy = vectEntropy.begin();
  double *ptVectCriteria = vectCriteria.begin();

  // Initialise vectCriteria
  for(j = 0; j < NBER_CRITERIA; j++)
    *(ptVectCriteria + j) = std::numeric_limits<double>::infinity();

  int nberLines = header? nberOfLines(xfileExploredModels) - 1 : nberOfLines(xfileExploredModels);

  if(lines.length() == 0)
    {
      lines = Rcpp::IntegerVector(nberLines);
      for(i = 0; i < nberLines; i++)
        *(lines.begin() + i) = i;
    }

  // Remove header
  if(header)
    nextLine(fIn, modelString);


  i = 0;
  int ii = 0, N;
  while(nextLine(fIn, modelString))
    {
      if(i == (*(lines.begin() + ii)))
        {
          std::istringstream tampon(modelString);

          if(!(tampon >> N))
            {
              fIn.close();
              throw Rcpp::exception("Incorrect value in the first column");
            }

          if(!(tampon >> xP))
            {
              fIn.close();
              throw Rcpp::exception("Incorrect value in the second column");
            }

          if(matS.ncol() != xP)
            {
              fIn.close();
              throw Rcpp::exception("Incorrect number of columns in 'matS'");
              return;
            }

          if(howmanyWords(modelString) < xP + NBER_MIN_COL_EXPLORED_MODELS_FILE)
            {
              fIn.close();
              throw Rcpp::exception("Incorrect number of column in file");
            }


          if(!(tampon >> xK))
            {
              fIn.close();
              throw Rcpp::exception("Incorrect number of population");
            }

          std::vector<int> xS(xP);

          for(j = 0; j < xP; j++)
            {
              if(!(tampon >> s) && ((s != 0) & (s!= 1)))
                {
                  fIn.close();
                  throw Rcpp::exception("Incorrect value for S");
                }
              xS[j] = s;
            }

          if(!(tampon >> logLik))
            {
              fIn.close();
              throw Rcpp::exception("Incorrect value of logLik");
            }

          if(!(tampon >> dim))
            {
              fIn.close();
              throw Rcpp::exception("Incorrect value of dimension");
            }

          if(!(tampon >> entropy))
            {
              fIn.close();
              throw Rcpp::exception("Incorrect value of Entropy");
            }

          Rcpp::DoubleVector xCriteria = computeCriteria_Rcpp(logLik, dim, N, entropy, cte);
          for(j = 0; j < NBER_CRITERIA; j++)
            {
              if(*(xCriteria.begin() + j) < *(ptVectCriteria + j))
                {
                  *(ptVectCriteria + j) = *(xCriteria.begin() + j);
                  *(ptVectK + j) = xK;
                  for(k = 0; k < xP; k++)
                    *(ptMatS + NBER_CRITERIA*k + j) = xS[k];
                  *(ptVectLogLik + j) = logLik;
                  *(ptVectDim + j) = dim;
                  *(ptVectEntropy + j) = entropy;
                }
            }

          ii++;
        }
      i++;
    }

  *ptVectN = N;
  fIn.close();
}

//  [[Rcpp::export()]]
void selectModelFromData_Rcpp(Rcpp::DoubleVector vectLogLik,
                              Rcpp::IntegerVector vectDim,
                              Rcpp::DoubleVector vectEntropy,
                              int N,
                              double cte,
                              Rcpp::IntegerVector vectIndexes,
                              Rcpp::DoubleVector vectCriteria)
{
  int nVectLogLik = vectLogLik.length();
  if((vectDim.length() != nVectLogLik) ||
     (vectEntropy.length() != nVectLogLik))
    {
      MyError("incompatible lengths of arguments");
      return;
    }
  int i, j;
  int *ptVectIndexes = vectIndexes.begin();

  double *ptVectCriteria = vectCriteria.begin();
  for(i = 0; i < NBER_CRITERIA; i++)
      *(ptVectCriteria + i) = std::numeric_limits<double>::infinity();

  for(i = 0; i < nVectLogLik; i++)
    {
      Rcpp::DoubleVector criteria = computeCriteria_Rcpp(vectLogLik[i],
                                                    vectDim[i],
                                                    N,
                                                    vectEntropy[i],
                                                    cte);
      for(j = 0; j < NBER_CRITERIA; j++)
        {
          if(criteria[j] < *(ptVectCriteria + j))
            {
              *(ptVectCriteria + j) = criteria[j];
              *(ptVectIndexes + j) = i;
            }
        }
    }

}

/**
 * @brief	Read a set of parameters of a model from a file.
 */
//	[[Rcpp::export()]]
Rcpp::List readParKS_Rcpp(std::string xfile)
{
  std::string lu, lu_1;
  int n, k, P, K;
  double p;
  std::ifstream fp(xfile.c_str(), std::ios_base::in);
  Rcpp::List outList;
  
  if(fp.bad())
  {
    ErrorOpeningFile();
    return outList;
  }

/*  try 
  {
		fp.open(xfile.c_str(), std::ios_base::in);
		fp.exceptions(fp.failbit);
  } catch (const std::ios_base::failure& e)
  {
		MyError("Opening a file");
		Rcpp::Rcout << "Caught an ios_base::failure.\n"
                  << "Explanatory string: " << e.what() << '\n';
		return outList;
  }
 */ 
  //
  if(nextLine(fp, lu)) // L
  {
    std::istringstream tampon (lu);
    tampon >> lu_1;
	
    if(!(lu_1=="P" && tampon >> P))
    {
      MyError("Incorrect value of number of loci\n");
      fp.close();
      return outList;
    }
  }
  else
  {
    MyError("Incomplete parameter file ");
    fp.close();
    return outList;
  }
  Rcpp::IntegerVector N_Levels(P);
  if(nextLine(fp, lu)) // 
  {
    std::istringstream tampon (lu);
    tampon >> lu_1;
    if(lu_1!="N_Levels")
    {
      MyError("Incorrect values of number of allele states\n");
      fp.close();
      return outList;
    }
    for(k=0; k < P; k++)
    {
      if(tampon >> n)
        N_Levels[k] = n;
      else
      {
        MyError("Incorrect value in N_Levels\n");
        fp.close();
        return outList;
      }
    }
    lu_1.erase();
  }
  else
  {
    MyError("Incomplete parameter file ");
    fp.close();
    return outList;
  }
  
  if(nextLine(fp, lu)) // K
  {
    std::istringstream tampon (lu);
    tampon >> lu_1;
    if(!(lu_1=="K" && tampon >> K && K>0))
    {
      MyError("Incorrect value of K\n");
      fp.close();
      return outList;
    }
  }
  else
  {
    MyError("Incomplete parameter file ");
    fp.close();
    return outList;
  }
  //cout << "\n > OK K\n";
  Rcpp::LogicalVector S(P);
  if(nextLine(fp, lu)) // S
  {
    std::istringstream tampon (lu);
    tampon >> lu_1;
    if(lu_1=="S")
    {
      for(k=0; k<P; k++)
      {
        if(tampon >> n && (n==1 || n==0))
        {
          S[k]= (n==1);
        }
        else
        {
          MyError("Incorrect value S\n");
          fp.close();
          return outList;
        }
      }
    }
    else
    {
      MyError("Incorrect value of S\n");
      fp.close();
      return outList;
    }
  }
  else
  {
    MyError("Incomplete parameter file ");
    fp.close();
    return outList;
  }
  
  Rcpp::DoubleVector pi_k(K);
  if(nextLine(fp, lu)) // pi_k
  {
    std::istringstream tampon (lu);
    tampon >> lu_1;
    if(lu_1=="mixingProportions")
    {
      //  Lecture des pi_k
      for(k = 0; k < K; k++)
      {
        if(tampon >> p) pi_k[k] = p;
        else
        {
          MyError("Incorrect value of mixing proportions\n");
          fp.close();
          return outList;
        }
      }
    }
  }
  else
  {
    MyError("Incomplete parameter file ");
    fp.close();
    return outList;
  }
  
  int j, i;
  //std::vector<int> variables_lu(P); // pour verifier si tous les loci sont lus
  int maxNRow = myMax<int>(N_Levels.begin(), P);
  
  Rcpp::NumericMatrix xmatrix(maxNRow, K);
  Rcpp::CharacterVector xCharactVector(maxNRow);
  
  Rcpp::List listProba(P);
  Rcpp::List listLevels(P);
  //std::fill(variables_lu.data(), variables_lu.data() + P, 0);
  
  j = 0;
  while(nextLine(fp, lu))
  {
	std::istringstream tampon (lu);
    tampon >> lu_1;
	
	if(lu_1 == "Variable" || lu_1 == "END")
	{
	  i = 0;
	}
	else
	{
	  xCharactVector[i] = lu_1;
	  if(S[j])
		{
		  for(k = 0; k < K; k++)
		  {
			if(!(tampon >> p && p >= 0.0 && p <= 1.0))
			{
			  MyError("Incorrect value of probability");
			  fp.close();
			  return outList;
			}
			xmatrix(i, k) = p;
		  }
		}
		else
		{
		  if(!(tampon >> p && p >= 0.0 && p <= 1.0))
			{
			  Rcout << "\n >> j = " << j+1 << "\n";
			  fp.close();
			  return outList;
			}
		  for(k = 0; k < K; k++)
			xmatrix(i, k) = p;
		}

	  i++;
	  if(i == N_Levels[j])
	  {
		listProba[j] = xmatrix(Range(0, N_Levels[j] - 1), _);
		listLevels[j] = xCharactVector;
		j++;
	  }
	}
	
	lu.erase();
	if(j == P) break;
  }
  
  fp.close();
  
  int dim, s = 0, sc = 0;

  for(j = 0; j < P; j++)
    {
      if(S[j]) s += N_Levels[j] - 1;
      else sc += N_Levels[j] - 1;
    }

  dim = K - 1 + K*s + sc;
  
  outList["P"] = P;
  outList["K"] = K;
  outList["S"] = S;
  outList["NbersLevels"] = N_Levels;
  outList["levels"] = listLevels;
  outList["pi_K"] = pi_k;
  outList["prob"] = listProba;
  outList["dim"] = dim;
  
  return outList;
}

//	[[Rcpp::export()]]
void writeParKS_InFile_Rcpp(Rcpp::List modelList, std::string file)
{
  PAR_KS par(modelList);
  par.printInFile(file);
}
