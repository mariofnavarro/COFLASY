/*
Set the Transfer rate term as function of (znu, zg)
*/

class Transfer2D{
private:
    /*
    Set variables used within the class
    */
    vector<double> myval_znu;
    vector<double> myval_zg;
    // nu e <-> nu e
    vector<double> myval_A12_nue_nue_I1;
    vector<double> myval_A34_nue_nue_I1;






    vector<double> vec_znu;
    vector<double> vec_zg;
    // nu e <-> nu e
    vector<double> vec_A12_nue_nue_I1;
    vector<double> vec_A34_nue_nue_I1;



    const double *double_znu;
    const double *double_zg;
    // nu e <-> nu e
    const double *double_A12_nue_nue_I1;
    const double *double_A34_nue_nue_I1;


    int num_znu;
    int num_zg;
    void *za_alloc;
    double* za;

    // nu e <-> nu e
    gsl_spline2d *spline_A12_nue_nue_I1;
    gsl_spline2d *spline_A34_nue_nue_I1;



    const gsl_interp2d_type *InterType;
    gsl_interp_accel *xacc;
    gsl_interp_accel *yacc;





public:
    /*
    Constructor of the class; read files of the Transfer Rates and set the shared interpolation settings
    */
    Transfer2D(void){
        ifstream myfile_znu;
        myfile_znu.open(_Path_ + "/data/Transfer_Rate/z.dat");
        if ( myfile_znu.is_open() ) {
            double tmp = 0.0;
            while (myfile_znu >> tmp) {
                myval_znu.push_back(tmp);
            }
        }
        myfile_znu.close();
        for (int i = 0; i < int(myval_znu.size()); i++){
            vec_znu.push_back(myval_znu[i]);
        }

        ifstream myfile_zg;
        myfile_zg.open(_Path_ + "/data/Transfer_Rate/z.dat");
        if ( myfile_zg.is_open() ) {
            double tmp = 0.0;
            while (myfile_zg >> tmp) {
                myval_zg.push_back(tmp);
            }
        }
        myfile_zg.close();
        for (int i = 0; i < int(myval_zg.size()); i++){
            vec_zg.push_back(myval_zg[i]);
        }


        /* ---- nu e <-> nu e*/

        // A12nuenue
        ifstream myfile_A12_nue_nue_I1;
        myfile_A12_nue_nue_I1.open(_Path_ + "/data/Transfer_Rate/A12_nue_nue_I1.dat");
        if ( myfile_A12_nue_nue_I1.is_open() ) {
            double tmp = 0.0;
            while (myfile_A12_nue_nue_I1 >> tmp) {
                myval_A12_nue_nue_I1.push_back(tmp);
            }
        }
        myfile_A12_nue_nue_I1.close();
        for (int i = 0; i < int(myval_A12_nue_nue_I1.size()); i++){
            vec_A12_nue_nue_I1.push_back(myval_A12_nue_nue_I1[i]);
        }

        // A34nuenue
        ifstream myfile_A34_nue_nue_I1;
        myfile_A34_nue_nue_I1.open(_Path_ + "/data/Transfer_Rate/A34_nue_nue_I1.dat");
        if ( myfile_A34_nue_nue_I1.is_open() ) {
            double tmp = 0.0;
            while (myfile_A34_nue_nue_I1 >> tmp) {
                myval_A34_nue_nue_I1.push_back(tmp);
            }
        }
        myfile_A34_nue_nue_I1.close();
        for (int i = 0; i < int(myval_A34_nue_nue_I1.size()); i++){
            vec_A34_nue_nue_I1.push_back(myval_A34_nue_nue_I1[i]);
        }




        double_znu      = &vec_znu[0];
        double_zg       = &vec_zg[0];
        num_znu         = vec_znu.size()*sizeof(double_znu) / sizeof(double); /* x grid points */
        num_zg          = vec_zg.size()*sizeof(double_zg) / sizeof(double); /* y grid points */
        za_alloc        = malloc(num_znu * num_zg * sizeof(double));
        za              = (double *)za_alloc;

        InterType       = gsl_interp2d_bilinear;
        xacc            = gsl_interp_accel_alloc();
        yacc            = gsl_interp_accel_alloc();


        // nu e <-> nu e
        spline_A12_nue_nue_I1      = gsl_spline2d_alloc(InterType, num_znu, num_zg);
        spline_A34_nue_nue_I1      = gsl_spline2d_alloc(InterType, num_znu, num_zg);

    }




    /**************************************************
    Set interpolation grid
    ***************************************************/

    /* ----- nu e <-> nu e  ---- */


    /*
    set the interpolation grid for the A12_nue_nue_I1 rate
    */
    void set_A12_nue_nue_I1(void){
        double_A12_nue_nue_I1    = &vec_A12_nue_nue_I1[0];
        for (int i = 0; i < int(myval_znu.size()); i++){
            for (int j = 0; j < int(myval_zg.size()); j++){
                gsl_spline2d_set(spline_A12_nue_nue_I1, za, i, j, double_A12_nue_nue_I1[int(myval_zg.size())*i+j]);
            }
        }
        /* initialize interpolation */
        gsl_spline2d_init(spline_A12_nue_nue_I1, double_znu, double_zg, za, num_znu, num_zg);
    }


    /*
    set the interpolation grid for the A34_nue_nue_I1 rate
    */
    void set_A34_nue_nue_I1(void){
        double_A34_nue_nue_I1    = &vec_A34_nue_nue_I1[0];
        for (int i = 0; i < int(myval_znu.size()); i++){
            for (int j = 0; j < int(myval_zg.size()); j++){
                gsl_spline2d_set(spline_A34_nue_nue_I1, za, i, j, double_A34_nue_nue_I1[int(myval_zg.size())*i+j]);
            }
        }
        /* initialize interpolation */
        gsl_spline2d_init(spline_A34_nue_nue_I1, double_znu, double_zg, za, num_znu, num_zg);
    }

    


    /**************************************************
    Declare methods returning the interpolated value
    ***************************************************/

    /* ----- nu e <-> nu e  ---- */

    /*
    Function returning the interpolated A12_nue_nue_I1 rate
    */
    double A12_nue_nue_I1(double znu, double zg){
        return gsl_spline2d_eval(spline_A12_nue_nue_I1, znu , zg, xacc, yacc);
    }

    /*
    Function returning the interpolated A34_nue_nue_I1 rate
    */
    double A34_nue_nue_I1(double znu, double zg){
        return gsl_spline2d_eval(spline_A34_nue_nue_I1, znu , zg, xacc, yacc);
    }

    
};


/*
set rates globally, accesible by various modules
*/
Transfer2D transfer2D;
/*
call this in the int main() (amiqs.cpp) to set the interpolation grid
*/
void set_transfer2D(void){
    // nu e <-> nu e
    transfer2D.set_A12_nue_nue_I1();
    transfer2D.set_A34_nue_nue_I1();
}