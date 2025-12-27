#ifndef OBJECTS_H
#define OBJECTS_H

//import dependencies
#include "functions.h"
#include "RandomNumber.h"
#include "EuropeanOption.h"
#include "CompoundOption.h"
#include "AmericanOption.h"
#include "AsianOption.h"

#include <map>


class Simulation{ //Factory class
private:
//control variables for change in attributes:
double S_hist = 100, r_hist = 0.05, sigma_hist = 0.2, T_hist = 1.0;
int DayNumber_hist = 252, SimulationNumber_hist = 1000;

    void update_attribute();

public:
    int seed = 42; //instance of random number

    //define default simmulation values 
    double S = 100, K = 105, T = 1.0, t = 0.0, sigma = 0.2;
    double r = 0.05, K2 = 5, T2 = 2.0, v_h = 0.01, s_h = 0.2 ;

    RandomNumber rand;


    int DayNumber = 252, SimulationNumber = 1000, n_threads = 1;

    Simulation();

    RandomNumber CreateRandomNumber(int seed);
    
    CompoundOption CreateCompoundOption();

    AsianOption CreateAsianOption();

    AmericanOption CreateAmericanOption();

    EuropeanOption CreateEuropeanOption();
};


#endif