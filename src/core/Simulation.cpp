
#include "Simulation.h"

#include <map>


void Simulation::update_attribute(){ //function updating the rand instance with the current attributes
    if (S_hist != S){rand.S = S; S_hist = S;}
    if (r_hist != r){rand.r = r; r_hist = r;}
    if (sigma_hist != sigma){rand.sigma = sigma; sigma_hist = sigma;}
    if (T_hist != T){rand.T = T; T_hist = T;}
    if (DayNumber_hist != DayNumber){rand.DayNumber = DayNumber; DayNumber_hist = DayNumber;}
    if (SimulationNumber_hist != SimulationNumber){rand.SimulationNumber = SimulationNumber; SimulationNumber_hist = SimulationNumber;}
}


    Simulation::Simulation()
        : rand(seed, 1000, 252, S, r, sigma, T){} //constructor automatically creates instance rand w/ default values

    
    //int seed, int SimulationNumber, int DayNumber, double S, double r, double sigma, double T)
    RandomNumber Simulation::CreateRandomNumber(int seed){
        update_attribute();
        RandomNumber newrand(seed, SimulationNumber, DayNumber, S, r, sigma, T);
        seed = seed + 1; //change seed for future instance of RandomNumber
        return newrand;
    }
    
    CompoundOption Simulation::CreateCompoundOption(){
        update_attribute();
        CompoundOption compound(rand, S, K, K2, T, T2, t, sigma, r, v_h, s_h);
        return compound;
    }

    AsianOption Simulation::CreateAsianOption(){
        update_attribute();
        AsianOption asian(rand, S, K, T, t, sigma, r, v_h, s_h, n_threads);
        return asian;
    }

    AmericanOption Simulation::CreateAmericanOption(){
        update_attribute();
        AmericanOption american(rand, S, K, T, t, sigma, r, v_h, s_h, n_threads); //initiate an instance of the EuropeanOption object
        return american;
    }

    EuropeanOption Simulation::CreateEuropeanOption(){
        EuropeanOption euro(rand, S, K, T, t, sigma, r, v_h, s_h); //initiate an instance of the EuropeanOption object
        return euro;
    }
