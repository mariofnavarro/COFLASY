class Interpolator1D {
public:
    Interpolator1D(const std::vector<double>& x_vals, const std::vector<double>& y_vals) {
        if (x_vals.size() != y_vals.size() || x_vals.size() < 2) {
            throw std::invalid_argument("Vectors must be the same size and contain at least 2 points.");
        }

        size = x_vals.size();
        acc = gsl_interp_accel_alloc();
        spline = gsl_spline_alloc(gsl_interp_cspline, size);
        gsl_spline_init(spline, x_vals.data(), y_vals.data(), size);
        xmin = x_vals.front();
        xmax = x_vals.back();
    }

    ~Interpolator1D() {
        gsl_spline_free(spline);
        gsl_interp_accel_free(acc);
    }

    double operator()(double x) const {
        if (x < xmin || x > xmax) {
            throw std::out_of_range("Interpolation query x is out of bounds.");
        }
        return gsl_spline_eval(spline, x, acc);
    }

private:
    size_t size;
    gsl_interp_accel* acc = nullptr;
    gsl_spline* spline = nullptr;
    double xmin, xmax;
};


struct Xn_ODE_Params {
    Interpolator1D* interp_Tnu;
    Interpolator1D* interp_Tg;
    Interpolator1D* interp_r11;
    Interpolator1D* interp_rb11;
    Interpolator1D* interp_H;
};
