/*
Set the weak rates as function of Tnu and Tgamma (znu, zg)
*/

class WeakRates2D{
private:
    /*
    Set variables used within the class
    */
    vector<double> myval_Tnu;
    vector<double> myval_Tg;
    // rates
    vector<double> myval_Gamma_n_penu4;
    vector<double> myval_Gamma_ne_pnu2;
    vector<double> myval_Gamma_nnu_pe;
    vector<double> myval_Gamma_pe_nnu2;
    vector<double> myval_Gamma_penu_n;
    vector<double> myval_Gamma_pnu_ne;





    vector<double> vec_Tnu;
    vector<double> vec_Tg;
    // rates
    vector<double> vec_Gamma_n_penu4;
    vector<double> vec_Gamma_ne_pnu2;
    vector<double> vec_Gamma_nnu_pe;
    vector<double> vec_Gamma_pe_nnu2;
    vector<double> vec_Gamma_penu_n;
    vector<double> vec_Gamma_pnu_ne;




    const double *double_Tnu;
    const double *double_Tg;
    // rates
    const double *double_Gamma_n_penu4;
    const double *double_Gamma_ne_pnu2;
    const double *double_Gamma_nnu_pe;
    const double *double_Gamma_pe_nnu2;
    const double *double_Gamma_penu_n;
    const double *double_Gamma_pnu_ne;






    int num_Tnu;
    int num_Tg;
    void *za_alloc;
    double* za;

    // rates
    gsl_spline2d *spline_Gamma_n_penu4;
    gsl_spline2d *spline_Gamma_ne_pnu2;
    gsl_spline2d *spline_Gamma_nnu_pe;
    gsl_spline2d *spline_Gamma_pe_nnu2;
    gsl_spline2d *spline_Gamma_penu_n;
    gsl_spline2d *spline_Gamma_pnu_ne;




    const gsl_interp2d_type *InterType;
    gsl_interp_accel *xacc;
    gsl_interp_accel *yacc;





public:
    /*
    Constructor of the class; read files of the Transfer Rates and set the shared interpolation settings
    */
    WeakRates2D(void){
        ifstream myfile_Tnu;
        myfile_Tnu.open(_Path_ + "/data/WeakRates/T.dat");
        if ( myfile_Tnu.is_open() ) {
            double tmp = 0.0;
            while (myfile_Tnu >> tmp) {
                myval_Tnu.push_back(tmp);
            }
        }
        myfile_Tnu.close();
        for (int i = 0; i < int(myval_Tnu.size()); i++){
            vec_Tnu.push_back(myval_Tnu[i]);
        }

        ifstream myfile_Tg;
        myfile_Tg.open(_Path_ + "/data/WeakRates/T.dat");
        if ( myfile_Tg.is_open() ) {
            double tmp = 0.0;
            while (myfile_Tg >> tmp) {
                myval_Tg.push_back(tmp);
            }
        }
        myfile_Tg.close();
        for (int i = 0; i < int(myval_Tg.size()); i++){
            vec_Tg.push_back(myval_Tg[i]);
        }


        /* ---- rates*/

        // n -> p+e+nu
        ifstream myfile_Gamma_n_penu4;
        myfile_Gamma_n_penu4.open(_Path_ + "/data/WeakRates/2D/Gamma_n_penu4.dat");
        if ( myfile_Gamma_n_penu4.is_open() ) {
            double tmp = 0.0;
            while (myfile_Gamma_n_penu4 >> tmp) {
                myval_Gamma_n_penu4.push_back(tmp);
            }
        }
        myfile_Gamma_n_penu4.close();
        for (int i = 0; i < int(myval_Gamma_n_penu4.size()); i++){
            vec_Gamma_n_penu4.push_back(myval_Gamma_n_penu4[i]);
        }

        // n+e -> p+nu
        ifstream myfile_Gamma_ne_pnu2;
        myfile_Gamma_ne_pnu2.open(_Path_ + "/data/WeakRates/2D/Gamma_ne_pnu2.dat");
        if ( myfile_Gamma_ne_pnu2.is_open() ) {
            double tmp = 0.0;
            while (myfile_Gamma_ne_pnu2 >> tmp) {
                myval_Gamma_ne_pnu2.push_back(tmp);
            }
        }
        myfile_Gamma_ne_pnu2.close();
        for (int i = 0; i < int(myval_Gamma_ne_pnu2.size()); i++){
            vec_Gamma_ne_pnu2.push_back(myval_Gamma_ne_pnu2[i]);
        }

        // n+nu -> p+e
        ifstream myfile_Gamma_nnu_pe;
        myfile_Gamma_nnu_pe.open(_Path_ + "/data/WeakRates/2D/Gamma_nnu_pe.dat");
        if ( myfile_Gamma_nnu_pe.is_open() ) {
            double tmp = 0.0;
            while (myfile_Gamma_nnu_pe >> tmp) {
                myval_Gamma_nnu_pe.push_back(tmp);
            }
        }
        myfile_Gamma_nnu_pe.close();
        for (int i = 0; i < int(myval_Gamma_nnu_pe.size()); i++){
            vec_Gamma_nnu_pe.push_back(myval_Gamma_nnu_pe[i]);
        }

        // p+e->n+nu
        ifstream myfile_Gamma_pe_nnu2;
        myfile_Gamma_pe_nnu2.open(_Path_ + "/data/WeakRates/2D/Gamma_pe_nnu2.dat");
        if ( myfile_Gamma_pe_nnu2.is_open() ) {
            double tmp = 0.0;
            while (myfile_Gamma_pe_nnu2 >> tmp) {
                myval_Gamma_pe_nnu2.push_back(tmp);
            }
        }
        myfile_Gamma_pe_nnu2.close();
        for (int i = 0; i < int(myval_Gamma_pe_nnu2.size()); i++){
            vec_Gamma_pe_nnu2.push_back(myval_Gamma_pe_nnu2[i]);
        }


        // p+e+nu->n
        ifstream myfile_Gamma_penu_n;
        myfile_Gamma_penu_n.open(_Path_ + "/data/WeakRates/2D/Gamma_penu_n.dat");
        if ( myfile_Gamma_penu_n.is_open() ) {
            double tmp = 0.0;
            while (myfile_Gamma_penu_n >> tmp) {
                myval_Gamma_penu_n.push_back(tmp);
            }
        }
        myfile_Gamma_penu_n.close();
        for (int i = 0; i < int(myval_Gamma_penu_n.size()); i++){
            vec_Gamma_penu_n.push_back(myval_Gamma_penu_n[i]);
        }


        // p+nu->n+e
        ifstream myfile_Gamma_pnu_ne;
        myfile_Gamma_pnu_ne.open(_Path_ + "/data/WeakRates/2D/Gamma_pnu_ne.dat");
        if ( myfile_Gamma_pnu_ne.is_open() ) {
            double tmp = 0.0;
            while (myfile_Gamma_pnu_ne >> tmp) {
                myval_Gamma_pnu_ne.push_back(tmp);
            }
        }
        myfile_Gamma_pnu_ne.close();
        for (int i = 0; i < int(myval_Gamma_pnu_ne.size()); i++){
            vec_Gamma_pnu_ne.push_back(myval_Gamma_pnu_ne[i]);
        }





        double_Tnu      = &vec_Tnu[0];
        double_Tg       = &vec_Tg[0];
        num_Tnu         = vec_Tnu.size()*sizeof(double_Tnu) / sizeof(double); /* x grid points */
        num_Tg          = vec_Tg.size()*sizeof(double_Tg) / sizeof(double); /* y grid points */
        za_alloc        = malloc(num_Tnu * num_Tg * sizeof(double));
        za              = (double *)za_alloc;

        InterType       = gsl_interp2d_bilinear;
        xacc            = gsl_interp_accel_alloc();
        yacc            = gsl_interp_accel_alloc();


        // rates
        spline_Gamma_n_penu4      = gsl_spline2d_alloc(InterType, num_Tnu, num_Tg);
        spline_Gamma_ne_pnu2      = gsl_spline2d_alloc(InterType, num_Tnu, num_Tg);
        spline_Gamma_nnu_pe      = gsl_spline2d_alloc(InterType, num_Tnu, num_Tg);
        spline_Gamma_pe_nnu2      = gsl_spline2d_alloc(InterType, num_Tnu, num_Tg);
        spline_Gamma_penu_n      = gsl_spline2d_alloc(InterType, num_Tnu, num_Tg);
        spline_Gamma_pnu_ne      = gsl_spline2d_alloc(InterType, num_Tnu, num_Tg);

    }




    /**************************************************
    Set interpolation grid
    ***************************************************/

    /* ----- rates  ---- */


    /*
    set the interpolation grid for the Gamma_n_penu4 rate
    */
    void set_Gamma_n_penu4(void){
        double_Gamma_n_penu4    = &vec_Gamma_n_penu4[0];
        for (int i = 0; i < int(myval_Tnu.size()); i++){
            for (int j = 0; j < int(myval_Tg.size()); j++){
                gsl_spline2d_set(spline_Gamma_n_penu4, za, i, j, double_Gamma_n_penu4[int(myval_Tg.size())*i+j]);
            }
        }
        /* initialize interpolation */
        gsl_spline2d_init(spline_Gamma_n_penu4, double_Tnu, double_Tg, za, num_Tnu, num_Tg);
    }


    /*
    set the interpolation grid for the Gamma_ne_pnu2 rate
    */
    void set_Gamma_ne_pnu2(void){
        double_Gamma_ne_pnu2    = &vec_Gamma_ne_pnu2[0];
        for (int i = 0; i < int(myval_Tnu.size()); i++){
            for (int j = 0; j < int(myval_Tg.size()); j++){
                gsl_spline2d_set(spline_Gamma_ne_pnu2, za, i, j, double_Gamma_ne_pnu2[int(myval_Tg.size())*i+j]);
            }
        }
        /* initialize interpolation */
        gsl_spline2d_init(spline_Gamma_ne_pnu2, double_Tnu, double_Tg, za, num_Tnu, num_Tg);
    }


    /*
    set the interpolation grid for the Gamma_nnu_pe rate
    */
    void set_Gamma_nnu_pe(void){
        double_Gamma_nnu_pe    = &vec_Gamma_nnu_pe[0];
        for (int i = 0; i < int(myval_Tnu.size()); i++){
            for (int j = 0; j < int(myval_Tg.size()); j++){
                gsl_spline2d_set(spline_Gamma_nnu_pe, za, i, j, double_Gamma_nnu_pe[int(myval_Tg.size())*i+j]);
            }
        }
        /* initialize interpolation */
        gsl_spline2d_init(spline_Gamma_nnu_pe, double_Tnu, double_Tg, za, num_Tnu, num_Tg);
    }


    /*
    set the interpolation grid for the Gamma_pe_nnu2 rate
    */
    void set_Gamma_pe_nnu2(void){
        double_Gamma_pe_nnu2    = &vec_Gamma_pe_nnu2[0];
        for (int i = 0; i < int(myval_Tnu.size()); i++){
            for (int j = 0; j < int(myval_Tg.size()); j++){
                gsl_spline2d_set(spline_Gamma_pe_nnu2, za, i, j, double_Gamma_pe_nnu2[int(myval_Tg.size())*i+j]);
            }
        }
        /* initialize interpolation */
        gsl_spline2d_init(spline_Gamma_pe_nnu2, double_Tnu, double_Tg, za, num_Tnu, num_Tg);
    }



    /*
    set the interpolation grid for the Gamma_penu_n rate
    */
    void set_Gamma_penu_n(void){
        double_Gamma_penu_n    = &vec_Gamma_penu_n[0];
        for (int i = 0; i < int(myval_Tnu.size()); i++){
            for (int j = 0; j < int(myval_Tg.size()); j++){
                gsl_spline2d_set(spline_Gamma_penu_n, za, i, j, double_Gamma_penu_n[int(myval_Tg.size())*i+j]);
            }
        }
        /* initialize interpolation */
        gsl_spline2d_init(spline_Gamma_penu_n, double_Tnu, double_Tg, za, num_Tnu, num_Tg);
    }



    /*
    set the interpolation grid for the Gamma_pnu_ne rate
    */
    void set_Gamma_pnu_ne(void){
        double_Gamma_pnu_ne    = &vec_Gamma_pnu_ne[0];
        for (int i = 0; i < int(myval_Tnu.size()); i++){
            for (int j = 0; j < int(myval_Tg.size()); j++){
                gsl_spline2d_set(spline_Gamma_pnu_ne, za, i, j, double_Gamma_pnu_ne[int(myval_Tg.size())*i+j]);
            }
        }
        /* initialize interpolation */
        gsl_spline2d_init(spline_Gamma_pnu_ne, double_Tnu, double_Tg, za, num_Tnu, num_Tg);
    }





    /**************************************************
    Declare methods returning the interpolated value
    ***************************************************/

    /* ----- rates  ---- */

    /*
    Function returning the interpolated Gamma_n_penu4 rate
    */
    double Gamma_n_penu4(double znu, double zg){
        return gsl_spline2d_eval(spline_Gamma_n_penu4, znu , zg, xacc, yacc);
    }

    /*
    Function returning the interpolated Gamma_ne_pnu2 rate
    */
    double Gamma_ne_pnu2(double znu, double zg){
        return gsl_spline2d_eval(spline_Gamma_ne_pnu2, znu , zg, xacc, yacc);
    }


    /*
    Function returning the interpolated Gamma_nnu_pe rate
    */
    double Gamma_nnu_pe(double znu, double zg){
        return gsl_spline2d_eval(spline_Gamma_nnu_pe, znu , zg, xacc, yacc);
    }


    /*
    Function returning the interpolated Gamma_pe_nnu2 rate
    */
    double Gamma_pe_nnu2(double znu, double zg){
        return gsl_spline2d_eval(spline_Gamma_pe_nnu2, znu , zg, xacc, yacc);
    }


    /*
    Function returning the interpolated Gamma_penu_n rate
    */
    double Gamma_penu_n(double znu, double zg){
        return gsl_spline2d_eval(spline_Gamma_penu_n, znu , zg, xacc, yacc);
    }


    /*
    Function returning the interpolated Gamma_pnu_ne rate
    */
    double Gamma_pnu_ne(double znu, double zg){
        return gsl_spline2d_eval(spline_Gamma_pnu_ne, znu , zg, xacc, yacc);
    }




};


/*
set rates globally, accesible by various modules
*/
WeakRates2D weakrates2D;
/*
call this in the int main() (amiqs.cpp) to set the interpolation grid
*/
void set_weakrates2D(void){
    // rates
    weakrates2D.set_Gamma_n_penu4();
    weakrates2D.set_Gamma_ne_pnu2();
    weakrates2D.set_Gamma_nnu_pe();
    weakrates2D.set_Gamma_pe_nnu2();
    weakrates2D.set_Gamma_penu_n();
    weakrates2D.set_Gamma_pnu_ne();
}