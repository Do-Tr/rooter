#ifndef Mcmc_H
#define Mcmc_H

#include <fstream>
#include <string>
#include <vector>
class BipartitionSamples;
class Model;
class RandomVariable;
class Settings;



class Mcmc {

    public:
                            Mcmc(void) = delete;
                            Mcmc(RandomVariable* r, Model* m, Settings* s);
                           ~Mcmc(void);
        void                run(void);
    
    protected:
        double              acceptanceProb(double lnX);
        std::vector<double> calculatePowers(int numStones, double alpha, double beta);
        std::string         formatTimeRemaining(double numSeconds);
        void                printSummaryFile(double lnL);
        void                sample(bool isFirstLine);
        void                sample(int gen, double lnL);
        
        Model*              model;
        RandomVariable*     rv;
        Settings*           settings;
        BipartitionSamples* rootSamples;
        std::ofstream       parmOut;
        std::ofstream       treeOut;
};

#endif
