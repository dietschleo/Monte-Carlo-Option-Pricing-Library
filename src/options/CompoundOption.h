#ifndef COMPOUNDOPTION_H
#define COMPOUNDOPTION_H

#include<vector>
#include<string>
#include<map>
#include "RandomNumber.h"



class CompoundOption {
private:
    double S_hist, sigma_hist, r_hist, T_hist;
    bool check_update();

//compound_option(double S, double K1, double K2, double T1, double T2, double sigma, double r, const std::vector<double>& z1, const std::vector<double>& z2)
public:
    RandomNumber rand;
    double S, K, K2, T, T2, t, sigma, r, s_h, v_h;
    std::map<std::string, double> prices, greeks;

    CompoundOption(RandomNumber rand, double S, double K, double K2, double T, double T2, double t, double sigma, double r,double v_h,double s_h);

    std::map<std::string, double> Price();

    std::map<std::string, double> ExtractGreeks();
    
};



#endif