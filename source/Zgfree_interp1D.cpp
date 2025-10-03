/*
Set interpolation function for 1D zg free evolution 
*/

class Zgfree{
private:
    /*
    Set variables used within the class
    */
    vector<double> myval_x;

    vector<double> myval_zg;


    /* --- */

    vector<double> vec_x;

    vector<double> vec_zg;
    



    /* --- */

    const double *double_x;

    const double *double_zg;



    /* --- */
    int num_x;


    /* --- */

    gsl_spline *spline_zg;



    /* --- */

    const gsl_interp_type *InterType;
    gsl_interp_accel *acc;

public:
    /*
    Constructor of the class; read files of the rates and set the shared interpolation settings
    */
    Zgfree(void){
        ifstream myfile_x;
        myfile_x.open(_Path_ + "/data/zg_free/z.dat");
        if ( myfile_x.is_open() ) {
            double tmp = 0.0;
            while (myfile_x >> tmp) {
                myval_x.push_back(tmp);
            }
        }
        myfile_x.close();
        for (int i = 0; i < int(myval_x.size()); i++){
            vec_x.push_back(myval_x[i]);
        }


        /* --- */
        ifstream myfile_zg;
        myfile_zg.open(_Path_ + "/data/zg_free/zg.dat");
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






        /* --- */
        double_x                = &vec_x[0];

        double_zg      = &vec_zg[0];


        /* --- */
        num_x                   = vec_x.size()*sizeof(double_x) / sizeof(double); 


        /* --- */
        InterType               = gsl_interp_linear;
        acc                     = gsl_interp_accel_alloc();

        /* --- */
        spline_zg      = gsl_spline_alloc(InterType, num_x);

    }


    /**************************************************
    initialize interpolation 
    ***************************************************/

    void set_zg(void){
        gsl_spline_init(spline_zg, double_x, double_zg, num_x);
    }


    /**************************************************
    Declare methods returning the interpolated value
    ***************************************************/

    double zg(double z){
        return gsl_spline_eval(spline_zg, z , acc);
    }



};


/*
set collision terms globally, accesible by various modules
*/
Zgfree zgfree;
/*
call this in the int main() (amiqs.cpp) to set the interpolation grid
*/
void set_zgfree(void){
    zgfree.set_zg();
}