
#include<vector>
#include<string>
#include<map>

#include "functions.h"
#include "AsianOption.h"

    
    bool AsianOption::check_update(){//function that returns true if a parameter has been updated, propagate parameter update to rand instance
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
    

    AsianOption::AsianOption(RandomNumber rand, double S, double K, double T, double t, double sigma, double r, double v_h, double s_h, int n_threads)
        : rand(rand), S(S), K(K), T(T), t(t), sigma(sigma), r(r), v_h(v_h), s_h(s_h), n_threads(n_threads) {}
    


    std::map<std::string, double> AsianOption::Price() {
        check_update();
        if(n_threads==1){
            rand.CreateBrownianMotion();
            prices = asian_options(S, r, T, K, rand.Znestedvect);
            return prices;

        } else {
            // Define a lambda function to wrap the call to asian_options
			auto asian_option_pricer = [this](double S, double r, double T, double K, const std::vector<std::vector<double>>& Znestedvect, int seed) {
			RandomNumber local_rand = rand;   // Create a local copy of rand
			local_rand.SetSeed(seed);     // Assign a unique seed for this thread
			local_rand.CreateBrownianMotion(); // Generate Brownian motion with the new seed
			return asian_options(S, r, T, K, local_rand.Znestedvect);
			};

			// Use parallel_run to compute results across multiple threads
			
			auto results = parallel_run(
			[this, asian_option_pricer](double S, double r, double T, double K, const std::vector<std::vector<double>>& Znestedvect, int thread_index) {
				int unique_seed = rand.GetSeed() + thread_index; // Unique seed for each thread
				return asian_option_pricer(S, r, T, K, Znestedvect, unique_seed);
			},
			n_threads,
			S, r, T, K, rand.Znestedvect, 0 // Thread index passed as the last argument
			);

			// Compute the average of the results
			std::map<std::string, double> averaged_results;
			for (const auto& [key, value] : results) {
			averaged_results[key] = value / static_cast<double>(n_threads); // Divide by the number of threads
			}
            prices = averaged_results;
			return averaged_results;
        }
    }

    void AsianOption::Clear(){
        prices.clear();
        greeks.clear();
    }
    
    std::map<std::string, double> AsianOption::ExtractGreeks(){
        //placeholders for current variables and prices
        if (prices.empty()){
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
