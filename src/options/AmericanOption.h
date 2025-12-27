#ifndef AMERICANOPTION_H
#define AMERICANOPTION_H

#include<vector>
#include<string>
#include<map>
#include "RandomNumber.h"



class AmericanOption{
private:
    double price;
    int n_threads;
    double S_hist, sigma_hist, r_hist, T_hist;

    bool check_update();

public:
    RandomNumber rand;
    double S, K, T, t, sigma, r, v_h, s_h;
    std::map<std::string, double> prices, greeks;

    AmericanOption(RandomNumber rand, double S, double K, double T, double t, double sigma, double r, double v_h, double s_h, int n_threads);

    double CallPrice();

    double PutPrice();

    std::map<std::string, double> Price();


    std::map<std::string, double> ExtractGreeks();

};



#endif