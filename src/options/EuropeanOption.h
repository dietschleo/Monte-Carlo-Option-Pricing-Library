#ifndef EUROPEANOPTION_H
#define EUROPEANOPTION_H

#include<vector>
#include<string>
#include<map>
#include "RandomNumber.h"

class EuropeanOption {
private:
    double S_hist, sigma_hist, r_hist, T_hist;
    bool MCdefault = false;
    std::map<std::string, double> MCprice; //cache for MC calculated price

    bool check_update();

public:
    RandomNumber rand;
    double S, K, T, t, sigma, r, s_h, v_h;
    std::map<std::string, double> prices, greeks;

    //constructor
    EuropeanOption(RandomNumber rand, double S, double K, double T, double t, double sigma, double r,double v_h,double s_h);

    //function for prices:

    std::map<std::string, double> Price();

    std::map<std::string, double> AnalyticalPrice();

    std::map<std::string, double> MonteCarloPrice();

    void Analytical(bool analytical);

    std::map<std::string, double> ExtractGreeks();


};


#endif