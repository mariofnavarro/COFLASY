/*
Set interpolation function for 1D Transfer rates 
*/

class Transfer1D{
private:
    /*
    Set variables used within the class
    */
    vector<double> myval_z;

    vector<double> myval_Aij_nunu_ee_I1;
    /* --- */

    vector<double> vec_z;

    vector<double> vec_Aij_nunu_ee_I1;

    /* --- */

    const double *double_z;

    const double *double_Aij_nunu_ee_I1;


    /* --- */
    int num_z;


    /* --- */

    gsl_spline *spline_Aij_nunu_ee_I1;

    /* --- */

    const gsl_interp_type *InterType;
    gsl_interp_accel *acc;


public:
    /*
    Constructor of the class; read files of the rates and set the shared interpolation settings
    */
    Transfer1D(void){
        ifstream myfile_z;
        myfile_z.open(_Path_ + "/data/Transfer_Rate/z.dat");
        if ( myfile_z.is_open() ) {
            double tmp = 0.0;
            while (myfile_z >> tmp) {
                myval_z.push_back(tmp);
            }
        }
        myfile_z.close();
        for (int i = 0; i < int(myval_z.size()); i++){
            vec_z.push_back(myval_z[i]);
        }


        /* --- */
        ifstream myfile_Aij_nunu_ee_I1;
        myfile_Aij_nunu_ee_I1.open(_Path_ + "/data/Transfer_Rate/Aij_nunu_ee_I1.dat");
        if ( myfile_Aij_nunu_ee_I1.is_open() ) {
            double tmp = 0.0;
            while (myfile_Aij_nunu_ee_I1 >> tmp) {
                myval_Aij_nunu_ee_I1.push_back(tmp);
            }
        }
        myfile_Aij_nunu_ee_I1.close();
        for (int i = 0; i < int(myval_Aij_nunu_ee_I1.size()); i++){
            vec_Aij_nunu_ee_I1.push_back(myval_Aij_nunu_ee_I1[i]);
        }



        /* --- */
        double_z                = &vec_z[0];

        double_Aij_nunu_ee_I1      = &vec_Aij_nunu_ee_I1[0];


        /* --- */
        num_z                   = vec_z.size()*sizeof(double_z) / sizeof(double); 


        /* --- */
        InterType               = gsl_interp_linear;
        acc                     = gsl_interp_accel_alloc();

        /* --- */
        spline_Aij_nunu_ee_I1      = gsl_spline_alloc(InterType, num_z);

    }


    /**************************************************
    initialize interpolation 
    ***************************************************/

    void set_Aij_nunu_ee_I1(void){
        gsl_spline_init(spline_Aij_nunu_ee_I1, double_z, double_Aij_nunu_ee_I1, num_z);
    }


    /**************************************************
    Declare methods returning the interpolated value
    ***************************************************/

    double Aij_nunu_ee_I1(double z){
        return gsl_spline_eval(spline_Aij_nunu_ee_I1, z , acc);
    }


};


/*
set collision terms globally, accesible by various modules
*/
Transfer1D transfer1D;
/*
call this in the int main() (amiqs.cpp) to set the interpolation grid
*/
void set_transfer1D(void){
    // nu nu -> e e
    transfer1D.set_Aij_nunu_ee_I1();
}