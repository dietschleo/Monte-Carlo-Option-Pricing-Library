
#include<vector>
#include<string>
#include<map>

#include "functions.h"
#include "RandomNumber.h"


    std::vector<std::vector<double>> RandomNumber::CreateRandomSeries(){
        Zvector = vector_std_dist(0.0, sigma, SimulationNumber, seed);
        Zvector2 = vector_std_dist(0.0, sigma, SimulationNumber, seed + 1);
        return {Zvector, Zvector2};
    }

    std::vector<std::vector<double>> RandomNumber::CreateBrownianMotion(){
        //check if the vector has never been computed or variables have been updated since last computation
        if (Znestedvect.size() == 0 || S != S_hist || sigma != sigma_hist || r != r_hist){
            Znestedvect = multi_simmulations_SBM(T, DayNumber, sigma, S, r, SimulationNumber, seed);
            S_hist = S;
            sigma_hist= sigma;
            r_hist=r;
        }
        return Znestedvect;
    }

    void RandomNumber::Reset(){
        //method to clear memory of instance
        Zvector.clear();
        Zvector2.clear();
        Znestedvect.clear();
    }

    void RandomNumber::SetSeed(int NewSeed){
        seed = NewSeed;
    }

    int RandomNumber::GetSeed(){//getter func for seed
        return seed;
    }
    RandomNumber::RandomNumber(int seed, int SimulationNumber, int DayNumber, double S, double r, double sigma, double T)
        : seed(seed), SimulationNumber(SimulationNumber), DayNumber(DayNumber), S(S), sigma(sigma), r(r), T(T){}

