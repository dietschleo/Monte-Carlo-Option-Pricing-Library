
#include<vector>
#include<string>
#include<map>

#include "functions.h"
#include "AmericanOption.h"



    bool AmericanOption::check_update(){//function that returns true if a parameter has been updated, propagate parameter update to rand instance
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

    AmericanOption::AmericanOption(RandomNumber rand, double S, double K, double T, double t, double sigma, double r, double v_h, double s_h, int n_threads)
        : rand(rand), S(S), K(K), T(T), t(t), sigma(sigma), r(r), v_h(v_h), s_h(s_h), n_threads(n_threads){}

    double AmericanOption::CallPrice(){
        check_update();
        rand.CreateBrownianMotion();
        prices["call"] = american_options(S, r, T, K, rand.Znestedvect, true);
        return prices["call"];
    }

    double AmericanOption::PutPrice(){
        check_update();
        rand.CreateBrownianMotion();
        prices["put"] = american_options(S, r, T, K, rand.Znestedvect, false);
        return prices["put"];
    }

    std::map<std::string, double> AmericanOption::Price() {
        if (n_threads == 1) {
            prices["call"] = CallPrice();
            prices["put"] = PutPrice();
        } else {
            // Define lambda function
            auto lambda_american_pricer = [this](double S, double r, double T, double K, const std::vector<std::vector<double>>& Znestedvect) -> std::map<std::string, double> {
                RandomNumber local_rand = rand;  // Create local copy of rand
                local_rand.CreateBrownianMotion(); // Generate Brownian motion
                std::map<std::string, double> temp;
                temp["call"] = american_options(S, r, T, K, local_rand.Znestedvect, true);
                temp["put"] = american_options(S, r, T, K, local_rand.Znestedvect, false);
                return temp;
            };

            // Multi-threading using parallel_run
            auto results = parallel_run(
                [this, lambda_american_pricer](double S, double r, double T, double K, const std::vector<std::vector<double>>& Znestedvect, int thread_index) {
                    RandomNumber local_rand = rand;  // Create a thread-local RandomNumber
                    local_rand.SetSeed(rand.GetSeed() + thread_index); // Assign unique seed
                    local_rand.CreateBrownianMotion(); // Generate Brownian motion
                    return lambda_american_pricer(S, r, T, K, local_rand.Znestedvect);
                },
                n_threads,
                S, r, T, K, rand.Znestedvect, 0 // Pass required arguments
            );

            // Average the results
            prices.clear();
            for (const auto& [key, value] : results) {
                prices[key] += value / static_cast<double>(n_threads);
            }
        }

        return prices;
    }


    std::map<std::string, double> AmericanOption::ExtractGreeks(){
        //placeholders for current variables and prices
        if (prices.empty() || !check_update()){
            prices = Price();
        }
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
                //std::cout << S << " Value for " << contract_greek << " is : " << results[contract_greek] << "\n";
            }
        }
        //reset original price and variable

        sigma = temp_sigma;
        S = temp_S;
        greeks = results;
        return results;
    }
