#ifndef RANDOMNUMBER_H  // Check if the header is not already included
#define RANDOMNUMBER_H 

#include<vector>
#include<string>
#include<map>

class RandomNumber{

private:    
    double S_hist, sigma_hist, r_hist, T_hist; //control variables for updates
        int seed;
public:
    int SimulationNumber, DayNumber;
    std::vector<double> Zvector, Zvector2;
    std::vector<std::vector<double>> Znestedvect;
    double S, sigma, r, T;

    std::vector<std::vector<double>> CreateRandomSeries();

    std::vector<std::vector<double>> CreateBrownianMotion();

    void Reset();

    void SetSeed(int NewSeed);

    int GetSeed();
    
    RandomNumber(int seed, int SimulationNumber, int DayNumber, double S, double r, double sigma, double T);


};

#endif