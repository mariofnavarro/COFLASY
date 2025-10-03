/*
Set interpolation function for Weak Rates Correstions
*/

class WeakRatesCorrection{
private:
    /*
    Set variables used within the class
    */
    vector<double> myval_T;

    vector<double> myval_Gamma_np_Corr;
    vector<double> myval_Gamma_pn_Corr;



    /* --- */

    vector<double> vec_T;

    vector<double> vec_Gamma_np_Corr;
    vector<double> vec_Gamma_pn_Corr;



    /* --- */

    const double *double_T;

    const double *double_Gamma_np_Corr;
    const double *double_Gamma_pn_Corr;



    /* --- */
    int num_T;


    /* --- */

    gsl_spline *spline_Gamma_np_Corr;
    gsl_spline *spline_Gamma_pn_Corr;


    /* --- */

    const gsl_interp_type *InterType;
    gsl_interp_accel *acc;

public:
    /*
    Constructor of the class; read files of the rates and set the shared interpolation settings
    */
    WeakRatesCorrection(void){
        ifstream myfile_T;
        myfile_T.open(_Path_ + "/data/WeakRates/Corrections/Tg.dat");
        if ( myfile_T.is_open() ) {
            double tmp = 0.0;
            while (myfile_T >> tmp) {
                myval_T.push_back(tmp);
            }
        }
        myfile_T.close();
        for (int i = 0; i < int(myval_T.size()); i++){
            vec_T.push_back(myval_T[i]);
        }


        /* --- */
        ifstream myfile_Gamma_np_Corr;
        myfile_Gamma_np_Corr.open(_Path_ + "/data/WeakRates/Corrections/Gamma_np_Corr.dat");
        if ( myfile_Gamma_np_Corr.is_open() ) {
            double tmp = 0.0;
            while (myfile_Gamma_np_Corr >> tmp) {
                myval_Gamma_np_Corr.push_back(tmp);
            }
        }
        myfile_Gamma_np_Corr.close();
        for (int i = 0; i < int(myval_Gamma_np_Corr.size()); i++){
            vec_Gamma_np_Corr.push_back(myval_Gamma_np_Corr[i]);
        }

        /* --- */
        ifstream myfile_Gamma_pn_Corr;
        myfile_Gamma_pn_Corr.open(_Path_ + "/data/WeakRates/Corrections/Gamma_pn_Corr.dat");
        if ( myfile_Gamma_pn_Corr.is_open() ) {
            double tmp = 0.0;
            while (myfile_Gamma_pn_Corr >> tmp) {
                myval_Gamma_pn_Corr.push_back(tmp);
            }
        }
        myfile_Gamma_pn_Corr.close();
        for (int i = 0; i < int(myval_Gamma_pn_Corr.size()); i++){
            vec_Gamma_pn_Corr.push_back(myval_Gamma_pn_Corr[i]);
        }






        /* --- */
        double_T                = &vec_T[0];

        double_Gamma_np_Corr    = &vec_Gamma_np_Corr[0];
        double_Gamma_pn_Corr    = &vec_Gamma_pn_Corr[0];


        /* --- */
        num_T                   = vec_T.size()*sizeof(double_T) / sizeof(double); 


        /* --- */
        InterType               = gsl_interp_linear;
        acc                     = gsl_interp_accel_alloc();

        /* --- */
        spline_Gamma_np_Corr    = gsl_spline_alloc(InterType, num_T);
        spline_Gamma_pn_Corr    = gsl_spline_alloc(InterType, num_T);



    }


    /**************************************************
    initialize interpolation 
    ***************************************************/

    void set_Gamma_np_Corr(void){
        gsl_spline_init(spline_Gamma_np_Corr, double_T, double_Gamma_np_Corr, num_T);
    }

    void set_Gamma_pn_Corr(void){
        gsl_spline_init(spline_Gamma_pn_Corr, double_T, double_Gamma_pn_Corr, num_T);
    }



    /**************************************************
    Declare methods returning the interpolated value
    ***************************************************/

    double Gamma_np_Corr(double z){
        return gsl_spline_eval(spline_Gamma_np_Corr, z , acc);
    }

    double Gamma_pn_Corr(double z){
        return gsl_spline_eval(spline_Gamma_pn_Corr, z , acc);
    }


};


/*
set weak rates globally, accesible by various modules
*/
WeakRatesCorrection weakratescorrection;
/*
call this in the int main() (amiqs.cpp) to set the interpolation grid
*/
void set_weakratescorrection(void){
    weakratescorrection.set_Gamma_np_Corr();
    weakratescorrection.set_Gamma_pn_Corr();
}