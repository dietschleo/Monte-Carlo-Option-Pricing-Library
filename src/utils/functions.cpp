// defining functions to be called in another file

#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <string>
#include <optional>
#include <map>
#include <type_traits>
#include <tuple>

//// Transverse functions and random numbers ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double discount(double v, double t, double r){
    return v * std::exp(-r * t);
}

std::vector<double> SBM(double T, int N, double sigma, double S, double r, int seed) {
    /* a function that return the vector of double for a single trajectory of a standard
    brownian motion of dimension N, maturity T and standard deviation sigma and starting value S*/

    double dt = T / static_cast<double>(N);
    double sqrt_dt = std::sqrt(dt);
    std::vector<double> trajectory(N, 0.0); //empty vector

    std::mt19937 generator(seed); //Mersenne Twister engine w/ random seed
    std::normal_distribution<double> normal_dist(0, sqrt_dt);

    trajectory[0] = S;
    for (int i=1; i<N; i++) {
        double z = normal_dist(generator);
        trajectory[i]=trajectory[i-1] * std::exp((r - 0.5 * sigma * sigma) * dt + sigma * z);
    //std::cout << trajectory[0] << " i : " << trajectory[i] << "\n";
    }

    return trajectory;
}

std::vector<std::vector<double>> multi_simmulations_SBM(double T, int N, double sigma, double S, double r, int Nmc, int seed) {
    /* a function that generates Nmc trajectories of SBM of dimension N, maturity T and std sigma */
    std::vector<std::vector<double>> trajectories; // nested vector to store Nmc trajectories
    //memory allocation
    trajectories.reserve(Nmc);

    //populating nested vectors with SBM
    for (int i = 0; i < Nmc; ++i) {
        trajectories.emplace_back(SBM(T, N, sigma, S, r, seed + i)); //move the vector from the function to the outer vector
    }

    return trajectories;
}

std::vector<double> vector_std_dist(double mean, double sigma, int Nmc, int seed) {
    /* function that generates Nmc number of normally distributed numbers
    */
    std::mt19937 generator(seed); // if seed is null use random seed
    std::normal_distribution<double> std_dist(mean, sigma);

    std::vector<double> z_vector(Nmc);

    for (int i = 0; i < Nmc; ++i) {
        z_vector[i] = std_dist(generator);
    }
    return z_vector;
}

double norm_cdf(double x) {
    /*  function that takes x a value of the standard deviation and returns the corresponding centile
    */
    return 0.5 * (1 + std::erf(x / std::sqrt(2.0)));
}

std::map<std::string, double> delta_gamma_extraction(double p, double p_plus, double p_minus, double h){
    /*Function that approximates delta and gamma with final central difference*/
    std::map <std::string, double> results;

    results["delta"] = (p_plus - p_minus) / (2 * h);
    results["gamma"] = (p_plus - 2 * p + p_minus) / (h * h);

    return results;
}

std::map<std::string, double> vega_vomma_extraction(double p, double p_plus, double p_minus, double h){
    /*Function that approximates vega and vomma with final central difference*/
    std::map <std::string, double> results;

    results["vega"] = (p_plus - p_minus) / (2 * h);
    results["vomma"] = (p_plus - 2 * p + p_minus) / (h * h);

    return results;
}

//// Option pricing functions //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::map<std::string, double> monte_carlo_european(double S, double K, double T, double t, double sigma, double r, const std::vector<double>& z_vector) {
    /* Function that calculates the average pay-off of a European call option
       given Nmc random numbers*/

    std::map <std::string, double> results;
    double A = 0;
    double B = 0;
    int Nmc = z_vector.size(); // Use the size of the input vector as Nmc
    double drift = (r - 0.5 * (sigma * sigma)) * (T - t); //drift component

    for (int i = 0; i < Nmc; i++) { // Iterate over all random numbers
        double diffusion = z_vector[i] * std::sqrt(T - t); //diffusion component
        //std::cout << z_vector[i] << std::endl;
        double St = S * std::exp(drift + diffusion);
        A += std::max(St - K, 0.0); // Add pay-off to total pay-off
        B += std::max(K - St, 0.0);
    }
    
    A = discount(A / Nmc, (T - t), r);
    B = discount(B / Nmc, (T - t), r);
    results["call"] = A;
    results["put"] = B;

    return results;
}

std::map<std::string, double> black_scholes_european(double S, double K, double T, double t, double sigma, double r) {
    /* function that takes as input <S> (spot price), <K> (strike price), the type of <contract>, <T> (time at maturity),
    <t> (time of valuation), <sigma> (implied volatility), <r> (the risk-free rate) and returns the price of the specified
    contract using the Black-Scholes analytical condition.  */
    
    double d1 = (std::log((S / K)) + (r + 0.5 * sigma * sigma) * (T - t)) / (sigma * std::sqrt((T - t))) ;
    double d2 = d1 - std::sqrt(T-t) * sigma ;
    
    std::map <std::string, double> results;
    results["call"] = S * norm_cdf(d1) - K * std::exp(-r *(T - t)) * norm_cdf(d2);
    results["put"] = K * std::exp(-r * (T - t)) * norm_cdf(-d2) - S * norm_cdf(-d1);

    return results;
}

std::map<std::string, double> compound_option(double S, double K1, double K2, double T1, double T2, double sigma, double r, const std::vector<double>& z1, const std::vector<double>& z2) {
    /* function that prices a compound option ie. call on call where K1/T1 are the strike and maturity of the outer option*/

    // Initialize map and keys
    std::map<std::string, double> results;
    results["call_on_call"] = 0.0;
    results["put_on_call"] = 0.0;
    results["put_on_put"] = 0.0;
    results["call_on_put"] = 0.0;

    int Nmc = z1.size();
    if(Nmc != z2.size()) {
        std::cerr << "Error: Random vector dimension mismatch" << std::endl;
        return {};
    }

    //Variable def
    double payoff = 0.0, sum_payoff = 0.0;
    std::vector<double> ST1(Nmc, 0.0), ST2(Nmc, 0.0);

    for(int i = 0; i<Nmc; ++i){
        //simulate undiscounted underlying prices at T1 and T2
        ST1[i] = S * std::exp((r - 0.5 * (sigma * sigma)) * T1 + sigma * std::sqrt(T1) * z1[i]);
        ST2[i] = ST1[i] * std::exp((r - 0.5*(sigma * sigma))*(T2 - T1) + sigma * std::sqrt(T2 - T1) * z2[i]);
        
        //compute inner payoff, discount to T1, computer outer payoff
        results["call_on_call"] += std::max(discount(std::max(ST2[i] - K2, 0.0), (T2 - T1), r) - K1, 0.0);
        results["put_on_call"] += std::max(K1 - discount(std::max(ST2[i] - K2, 0.0), (T2 - T1), r), 0.0);      
        results["put_on_put"] += std::max(K1 - discount(std::max(K2 - ST2[i], 0.0), (T2 - T1), r), 0.0);
        results["call_on_put"] += std::max(discount(std::max(K2 - ST2[i], 0.0), (T2 - T1), r) - K1, 0.0);
    }

    for (auto& pair : results) { //discount sum of payoff to t=0 using average linearity property
        pair.second = discount(pair.second  / Nmc, T1, r);
    }
    return results;
}    

std::map<std::string, double> asian_options(double S, double r, double T, double K, std::vector<std::vector<double>>& underlying) {
    /* function that returns the price of an asian option based on underlying price, risk free rate, strike price for each contract type
    (ie. fixed strike call, fixed strike put, floating strike call, floating strike put)*/

    std::map<std::string, double> results;
    double Nmc = underlying.size(); //nb of paths
    double N = underlying[0].size();//length of paths

    double fi_call = 0.0, fi_put = 0.0, fl_call = 0.0, fl_put = 0.0;

    if (underlying.empty() || underlying[0].empty()) {
        std::cerr << "Error: 'underlying' data is empty.\n";
        return results;
    }


    for(int sim = 0; sim < Nmc; ++sim){
        double A = 0.0;  //average underlying price
        double ST = 0.0; //underlying price at maturity
        for(int i = 0; i < N; ++i){
            A += underlying[sim][i];
        }
        ST = underlying[sim][N-1];
        A=A/N;
        fi_call += std::max((A - K), 0.0);
        fi_put += std::max(K - A, 0.0);
        fl_call += std::max(A - ST, 0.0);
        fl_put += std::max(ST - A, 0.0);
    }

    results["fixed_strike_call"] = discount(fi_call/Nmc, T, r);
    results["fixed_strike_put"] = discount(fi_put/Nmc, T, r);
    results["floating_strike_call"] = discount(fl_call/Nmc, T, r);
    results["floating_strike_put"] = discount(fl_put/Nmc, T, r);

    return results;
}

//// Functions related to american options /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double mean(std::vector<double> series){
    //I can't believe I need to code this.
    double mean = 0.0, N = series.size();
    for(int i = 0; i < N; i++){
        mean += series[i];
    }
    return mean / N;
}

std::map<std::string, double> OLSregression(const std::vector<double> X, std::vector<double> Y){
    //function that returns slope and intercept of an OLS SLR between X and Y series

    std::map<std::string, double> results;
    double N = X.size();

    if(N != Y.size()){
        std::cout << "Error: expected two series of the same size.";
        return results;
    }
    double Xmean = mean(X);
    double Ymean = mean(Y);
    //using: slope = Cov(X,Y)/Var(X)
    double Cov = 0.0;
    double Var = 0.0;

    for(int i=0; i < N; ++i){
        Cov += (X[i] - Xmean) * (Y[i] - Ymean);
            double temp = (X[i] - Xmean); //placeholder to prevent computing twice
            Var += temp * temp;
    }
    results["slope"] = Cov / Var;
    results["intercept"] = Ymean - results["slope"] * Xmean;
    return results;
}

std::map<std::string, std::vector<double>> is_in_the_money(std::vector<double>& underlying, const std::vector<int>& exercise_time, double K, bool iscall){
    // function that returns the prices and their index when they're in the money and haven't been exercised yet
    int Nmc = underlying.size();
    std::map<std::string, std::vector<double>> results;
    std::vector<double> ITM(Nmc), index(Nmc);
    int count_itm = 0;    
  
    for (int i = 0; i < Nmc; ++i) {
        //check if the option is ITM and hasn't been exercised yet
        if(((underlying[i] > K && iscall) || (underlying[i] < K && !iscall)) && exercise_time[i] == -1){
            ITM[count_itm] = underlying[i];
            index[count_itm] = (double)i;
            count_itm += 1;
        }
    }

    index.resize(count_itm);
    ITM.resize(count_itm);
    results["index"] = index;
    results["price"] = ITM;
    return results;
}

double american_options(double S, double r, double T, double K, std::vector<std::vector<double>>& underlying, bool iscall) {
    //function that computes an american option price using the Longstaff-Schwartz method. if iscall is true, the pricer will return the call price. 
    //returns a double instead of a map since it's not computationally interesting to compute put and call prices all at once.

    
    double payoff = 0.0;
    int Nmc = underlying.size(); //nb of paths
    double N = underlying[0].size();//length of paths
    double dt = T / N; //time step size
    double inv_Nmc = 1.0 / Nmc; //weight of a single sim

    std::vector<int> exercise_time(Nmc, -1); //a vector to store the timestep options are exercised.

    if (underlying.empty() || underlying[0].empty()) {
        std::cerr << "Error: 'underlying' data is empty.\n";
        return -1;
    }
    
    //LSM method
    for(int t = N - 2; t >= 0; t--){ //loop backward
        
        //note for later: the cross section vectors could probably be pointers 
        std::vector<double> cross_section(Nmc);
        std::vector<double> CS_payoff;

        //build cross-sectional vector at t
        for (int j = 0; j < Nmc; ++j){cross_section[j] = underlying[j][t];}
        //identify in the money paths 
        std::map<std::string, std::vector<double>> ITM = is_in_the_money(cross_section, exercise_time, K, iscall);

        int N_ITM = ITM["index"].size();
        if (N_ITM == 0){break;}

        //cross-sectional payoff (at t+1) not exercised 
        if (iscall){
            for(int i = 0; i < N_ITM; ++i){//iterate on ITM
                CS_payoff.emplace_back( //discounted ITM cf in t+1
                    discount(
                        std::max(underlying[ITM["index"][i]][t + 1] - K, 0.0),
                        dt, r));
            }
        } else {for(int i = 0; i < N_ITM; ++i){//iterate on ITM
                CS_payoff.emplace_back( //discounted ITM cf in t+1
                    discount(
                        std::max(K - underlying[ITM["index"][i]][t + 1], 0.0),
                        dt, r));
            }
        }
        std::map<std::string, double> reg = OLSregression(ITM["price"], CS_payoff);
        std::vector<double> continuation_val(N_ITM);
        //estimate continuation value
        for(int i = 0; i < N_ITM; ++i){
            continuation_val[i] = reg["intercept"] + reg["slope"] * ITM["price"][i];
            //check if it's more profitable to exercise now
            if((iscall && std::max(ITM["price"][i] - K, 0.0) > continuation_val[i]) //call
            ||(!iscall && std::max(K - ITM["price"][i], 0.0) > continuation_val[i])){ //put
                exercise_time[ITM["index"][i]] = t; //exercised
            }
        }
    }
    
    for (int i = 0; i< Nmc; ++i){
        if(exercise_time[i] >= 0){ //don't sum unexercized options
            if(iscall){payoff += inv_Nmc * discount(
                    std::max(underlying[i][exercise_time[i]] - K, 0.0), 
                    exercise_time[i] * dt, r); 
            } else {payoff += inv_Nmc * discount(
                    std::max(K - underlying[i][exercise_time[i]], 0.0), 
                    exercise_time[i] * dt, r);
            }
        }
    }
    return payoff;
}

