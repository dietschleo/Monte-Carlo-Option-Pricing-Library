#ifndef FUNCTIONS_H  // Check if the header is not already included
#define FUNCTIONS_H  // Define a unique identifier for this header
#include <future>
#include <functional>
#include <random>
#include <map>

// Function declaration
template <typename Function, typename... Args>
auto parallel_run(Function pricing, int n_thread, Args... args) -> decltype(pricing(args...)){
    /*function that taks the pricing function and its arguments and computes the average of its results on multiple threads*/
    using ReturnType = decltype(pricing(args...));
    std::vector<std::future<ReturnType>> futures;

    for (int i = 0; i < n_thread; i++){ //launch the threads
        futures.emplace_back(std::async(std::launch::async, pricing, args...));
    }
    std::map<std::string, double> results;
    for (auto& future : futures){
        auto result = future.get();
        for (const auto& [key, value] : result){
            results[key] += value;
        }
    }
    return results;
}


double discount(double v, double t, double r);
std::vector<double> SBM(double T, int N, double sigma);
std::vector<std::vector<double>> multi_simmulations_SBM(double T, int N, double sigma, double S, double r, int Nmc, int seed);
std::vector<double> vector_std_dist(double mean, double sigma, int Nmc, int seed);
std::map<std::string, double> monte_carlo_european(double S, double K, double T, double t, double sigma, double r, const std::vector<double>& z_vector);
std::map<std::string, double> black_scholes_european(double S, double K, double T, double t, double sigma, double r);
double norm_cdf(double x);
std::map<std::string, double> compound_option(double S, double K1, double K2, double T1, double T2, double sigma, double r, const std::vector<double>& z1, const std::vector<double>& z2);
std::map<std::string, double> delta_gamma_extraction(double p, double p_plus, double p_minus, double h);
std::map<std::string, double> vega_vomma_extraction(double p, double p_plus, double p_minus, double h);
std::map<std::string, double> asian_options(double S, double K, double T, double r, std::vector<std::vector<double>>& underlying);

double mean(std::vector<double> series);
std::map<std::string, double> OLSregression(const std::vector<double> X, std::vector<double> Y);
std::map<std::string, std::vector<double>> is_in_the_money(std::vector<double>& underlying, const std::vector<int>& exercise_time, double K, bool iscall);
double american_options(double S, double r, double T, double K, std::vector<std::vector<double>>& underlying, bool iscall);



#endif