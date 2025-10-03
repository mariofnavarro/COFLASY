/*
Set interpolation function for 1D Weak Rates 
*/

class WeakRates1D{
private:
    /*
    Set variables used within the class
    */
    vector<double> myval_T;

    vector<double> myval_Gamma_ne_pnu1;
    vector<double> myval_Gamma_n_penu2;
    vector<double> myval_Gamma_n_penu3;

    vector<double> myval_Gamma_pe_nnu1;


    /* --- */

    vector<double> vec_T;

    vector<double> vec_Gamma_ne_pnu1;
    vector<double> vec_Gamma_n_penu2;
    vector<double> vec_Gamma_n_penu3;

    vector<double> vec_Gamma_pe_nnu1;
    



    /* --- */

    const double *double_T;

    const double *double_Gamma_ne_pnu1;
    const double *double_Gamma_n_penu2;
    const double *double_Gamma_n_penu3;

    const double *double_Gamma_pe_nnu1;



    /* --- */
    int num_T;


    /* --- */

    gsl_spline *spline_Gamma_ne_pnu1;
    gsl_spline *spline_Gamma_n_penu2;
    gsl_spline *spline_Gamma_n_penu3;

    gsl_spline *spline_Gamma_pe_nnu1;


    /* --- */

    const gsl_interp_type *InterType;
    gsl_interp_accel *acc;

public:
    /*
    Constructor of the class; read files of the rates and set the shared interpolation settings
    */
    WeakRates1D(void){
        ifstream myfile_T;
        myfile_T.open(_Path_ + "/data/WeakRates/T.dat");
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
        ifstream myfile_Gamma_ne_pnu1;
        myfile_Gamma_ne_pnu1.open(_Path_ + "/data/WeakRates/1D/Gamma_ne_pnu1.dat");
        if ( myfile_Gamma_ne_pnu1.is_open() ) {
            double tmp = 0.0;
            while (myfile_Gamma_ne_pnu1 >> tmp) {
                myval_Gamma_ne_pnu1.push_back(tmp);
            }
        }
        myfile_Gamma_ne_pnu1.close();
        for (int i = 0; i < int(myval_Gamma_ne_pnu1.size()); i++){
            vec_Gamma_ne_pnu1.push_back(myval_Gamma_ne_pnu1[i]);
        }

        /* --- */
        ifstream myfile_Gamma_n_penu2;
        myfile_Gamma_n_penu2.open(_Path_ + "/data/WeakRates/1D/Gamma_n_penu2.dat");
        if ( myfile_Gamma_n_penu2.is_open() ) {
            double tmp = 0.0;
            while (myfile_Gamma_n_penu2 >> tmp) {
                myval_Gamma_n_penu2.push_back(tmp);
            }
        }
        myfile_Gamma_n_penu2.close();
        for (int i = 0; i < int(myval_Gamma_n_penu2.size()); i++){
            vec_Gamma_n_penu2.push_back(myval_Gamma_n_penu2[i]);
        }



        /* --- */
        ifstream myfile_Gamma_n_penu3;
        myfile_Gamma_n_penu3.open(_Path_ + "/data/WeakRates/1D/Gamma_n_penu3.dat");
        if ( myfile_Gamma_n_penu3.is_open() ) {
            double tmp = 0.0;
            while (myfile_Gamma_n_penu3 >> tmp) {
                myval_Gamma_n_penu3.push_back(tmp);
            }
        }
        myfile_Gamma_n_penu3.close();
        for (int i = 0; i < int(myval_Gamma_n_penu3.size()); i++){
            vec_Gamma_n_penu3.push_back(myval_Gamma_n_penu3[i]);
        }


        /* --- */
        ifstream myfile_Gamma_pe_nnu1;
        myfile_Gamma_pe_nnu1.open(_Path_ + "/data/WeakRates/1D/Gamma_pe_nnu1.dat");
        if ( myfile_Gamma_pe_nnu1.is_open() ) {
            double tmp = 0.0;
            while (myfile_Gamma_pe_nnu1 >> tmp) {
                myval_Gamma_pe_nnu1.push_back(tmp);
            }
        }
        myfile_Gamma_pe_nnu1.close();
        for (int i = 0; i < int(myval_Gamma_pe_nnu1.size()); i++){
            vec_Gamma_pe_nnu1.push_back(myval_Gamma_pe_nnu1[i]);
        }



        /* --- */
        double_T                = &vec_T[0];

        double_Gamma_ne_pnu1    = &vec_Gamma_ne_pnu1[0];
        double_Gamma_n_penu2    = &vec_Gamma_n_penu2[0];
        double_Gamma_n_penu3    = &vec_Gamma_n_penu3[0];

        double_Gamma_pe_nnu1    = &vec_Gamma_pe_nnu1[0];


        /* --- */
        num_T                   = vec_T.size()*sizeof(double_T) / sizeof(double); 


        /* --- */
        InterType               = gsl_interp_linear;
        acc                     = gsl_interp_accel_alloc();

        /* --- */
        spline_Gamma_ne_pnu1    = gsl_spline_alloc(InterType, num_T);
        spline_Gamma_n_penu2    = gsl_spline_alloc(InterType, num_T);
        spline_Gamma_n_penu3    = gsl_spline_alloc(InterType, num_T);

        spline_Gamma_pe_nnu1    = gsl_spline_alloc(InterType, num_T);


    }


    /**************************************************
    initialize interpolation 
    ***************************************************/

    void set_Gamma_ne_pnu1(void){
        gsl_spline_init(spline_Gamma_ne_pnu1, double_T, double_Gamma_ne_pnu1, num_T);
    }

    void set_Gamma_n_penu2(void){
        gsl_spline_init(spline_Gamma_n_penu2, double_T, double_Gamma_n_penu2, num_T);
    }

    void set_Gamma_n_penu3(void){
        gsl_spline_init(spline_Gamma_n_penu3, double_T, double_Gamma_n_penu3, num_T);
    }



    void set_Gamma_pe_nnu1(void){
        gsl_spline_init(spline_Gamma_pe_nnu1, double_T, double_Gamma_pe_nnu1, num_T);
    }


    /**************************************************
    Declare methods returning the interpolated value
    ***************************************************/

    double Gamma_ne_pnu1(double z){
        return gsl_spline_eval(spline_Gamma_ne_pnu1, z , acc);
    }


    double Gamma_n_penu2(double z){
        return gsl_spline_eval(spline_Gamma_n_penu2, z , acc);
    }

    double Gamma_n_penu3(double z){
        return gsl_spline_eval(spline_Gamma_n_penu3, z , acc);
    }


    double Gamma_pe_nnu1(double z){
        return gsl_spline_eval(spline_Gamma_pe_nnu1, z , acc);
    }
};


/*
set weak rates globally, accesible by various modules
*/
WeakRates1D weakrates1D;
/*
call this in the int main() (amiqs.cpp) to set the interpolation grid
*/
void set_weakrates1D(void){
    // n+e->p+nu
    weakrates1D.set_Gamma_ne_pnu1();
    // n -> p+e+nu
    weakrates1D.set_Gamma_n_penu2();
    weakrates1D.set_Gamma_n_penu3();

    // p+nu -> n+e
    weakrates1D.set_Gamma_pe_nnu1();
}