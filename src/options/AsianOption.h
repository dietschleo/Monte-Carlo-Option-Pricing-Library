#ifndef ASIANOPTION_H
#define ASIANOPTION_H

#include<vector>
#include<string>
#include<map>
#include "RandomNumber.h"




class AsianOption {
private:
    double S_hist, sigma_hist, r_hist, T_hist;
    
    bool check_update();
    
public:
    RandomNumber rand;
    int n_threads;
    double S, K, T, t, sigma, r, v_h, s_h;
    std::map<std::string, double> prices, greeks;

    AsianOption(RandomNumber rand, double S, double K, double T, double t, double sigma, double r, double v_h, double s_h, int n_threads);


    std::map<std::string, double> Price();

    void Clear();

    std::map<std::string, double> ExtractGreeks();
};



#endif