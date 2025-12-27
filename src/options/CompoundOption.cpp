
#include<vector>
#include<string>
#include<map>

#include "functions.h"
#include "CompoundOption.h"



    bool CompoundOption::check_update(){//function that returns true if a parameter has been updated, propagate parameter update to rand instance
        bool result = true; //true until proven false
        if(S_hist != S || sigma_hist != sigma || r_hist != r || T_hist != T){
            rand.S = S;
            rand.r = r;
            rand.sigma = sigma;
            S_hist = S;
            r_hist = r;
            sigma_hist = sigma;
            result = false;
            rand.Reset();
        }
        return result;
    }


    CompoundOption::CompoundOption(RandomNumber rand, double S, double K, double K2, double T, double T2, double t, double sigma, double r,double v_h,double s_h)
    : rand(rand), S(S), K(K), K2(K2), T(T), T2(T2), t(t), sigma(sigma), r(r), v_h(v_h), s_h(s_h) {}

    std::map<std::string, double> CompoundOption::Price(){
        check_update();
        prices.clear();
        rand.Reset();
        rand.CreateRandomSeries();
        prices = compound_option(S, K, K2, T, T2, sigma, r, rand.Zvector, rand.Zvector2);
        return prices;
    }

    std::map<std::string, double> CompoundOption::ExtractGreeks(){
        if (prices.empty() || !check_update()){
            prices = Price();
        }
        //placeholders for current variables and prices
        std::map<std::string, double> temp_prices = prices;
        double temp_S = S;
        double temp_sigma = sigma;
        S = S + s_h;
        check_update();
		std::map<std::string, double> p_plus = Price();
        S = temp_S - s_h;
		check_update();
        std::map<std::string, double> p_minus = Price();
        sigma = sigma + v_h;
		check_update();
        std::map<std::string, double>  p_plus_v = Price();
        sigma = temp_sigma - v_h;
		check_update();
        std::map<std::string, double>  p_minus_v = Price();

        prices = temp_prices;
        std::map<std::string, double> results;
        for (auto& pair : prices) { //for each type of contract:
            results[pair.first + "_price"] = pair.second;
            std::map<std::string, double> greeks = delta_gamma_extraction(prices[pair.first], p_plus[pair.first], p_minus[pair.first], s_h);
            greeks["vega"] = vega_vomma_extraction(prices[pair.first], p_plus_v[pair.first], p_minus_v[pair.first], v_h)["vega"];
            greeks["vomma"] = vega_vomma_extraction(prices[pair.first], p_plus_v[pair.first], p_minus_v[pair.first], v_h)["vomma"];
            for (auto& duo : greeks) { //for each greek:
                std::string contract_greek = pair.first + "_" + duo.first;
                results[contract_greek] = duo.second;
                //std::cout << "Value for " << contract_greek << " is : " << results[contract_greek] << "\n";
            }
        }
        //reset original price and variable
        
        sigma = temp_sigma;
        S = temp_S;
        greeks = results;
        return results;
    }
    
