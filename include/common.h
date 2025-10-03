#include <vector>
#include <cmath>
#include <iostream>
#include <ctime>   // localtime
#include <sstream> // stringstream
#include <iomanip> // put_time
#include <string>  // string



#define _SUCCESS_ 0 /**< integer returned after successful call of a function */
#define _FAIL_ 1 /**< integer returned after failed call of a function */


const std::string _Path_= "..";


const double _PI_       = 3.141592653589793238462643383279;
const double _E_        = 2.7182818284590452353602874713527;


const double _MPL_      = 7.354086028747625*pow(10.,21.);
const double _gstarrad_ = 10.75;
const double _Tref_     = 20.0;
const double _MW_       = 80.377*pow(10,3.);
const double _Me_       = 0.511;
const double _Mmu_      = 105.7;
const double _GF_       = 1.166*pow(10,-11.);
const double _Zeta3_    = 1.2020569031595942853997381615114;
const double _ThetaW_   = 0.491878;
const double _Mn_       = 939.56542194;
const double _Mp_       = 938.27208943;
const double _Q_        = 1.2933325099999138; /* Q = Mn - Mp*/
const double _hbar_     = 6.582119569*pow(10.,-16.);
const double _taunsec_  = 878.4;

// hierarchy independent
const double _dm21_     = 7.530e-17;
const double _theta12_  = 0.587;
const double _theta13_  = 0.148;
// hierarchy dependent
const double _dm31NH_   = 2.534e-15;
const double _theta23NH_= 0.831;
const double _dm31IH_   = -2.46e-15;
const double _theta23IH_= 0.824;


