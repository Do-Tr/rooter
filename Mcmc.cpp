#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include "BipartitionSamples.h"
#include "EigenSystem.h"
#include "Mcmc.h"
#include "Model.h"
#include "Move.h"
#include "ParameterTree.h"
#include "RandomVariable.h"
#include "Settings.h"
#include "SteppingStones.h"



Mcmc::Mcmc(RandomVariable* r, Model* m, Settings* s) {

    model = m;
    rv = r;
    settings = s;
    rootSamples = new BipartitionSamples();
}

Mcmc::~Mcmc(void) {

    delete rootSamples;
}

double Mcmc::acceptanceProb(double lnX) {

    double prob = 0.0;
    if (lnX < -300.0)
        prob = 0.0;
    else if (lnX > 0.0)
        prob = 1.0;
    else
        prob = exp(lnX);
    return prob;
}

std::vector<double> Mcmc::calculatePowers(int numStones, double alpha, double beta) {

    int ns = numStones - 1;
    std::vector<double> pwrs;
    double intervalProb = (double)1.0 / ns;
    pwrs.push_back(1.0);
    for (int i=ns-1; i>0; i--)
        pwrs.push_back( rv->betaQuantile(alpha, beta, i * intervalProb) );
    pwrs.push_back(0.0);
#   if 0
    for (int i=0; i<pwrs.size(); i++)
        std::cout << i+1 << " -- " << std::fixed << std::setprecision(15) << pwrs[i] << std::endl;
#   endif
    return pwrs;
}

std::string Mcmc::formatTimeRemaining(double numSeconds) {

    int numDays = numSeconds / 86400;
    numSeconds -= 86400 * numDays;
    
    int numHours = numSeconds / 3600;
    numSeconds -= 3600 * numHours;
    
    int numMinutes = numSeconds / 60;
    numSeconds -= 60 * numMinutes;
    
    int ns = numSeconds;
    
    std::string str = "";
    if (numDays > 0)
        str += std::to_string(numDays) + "d ";
    if ( numHours > 0 || (numDays > 0 && numHours == 0) )
        str += std::to_string(numHours) + ":";
    if (numMinutes < 10)
        str += "0";
    str += std::to_string(numMinutes) + ":";
    if (ns < 10)
        str += "0";
    str += std::to_string(ns);
    
    return str;
}

void Mcmc::printSummaryFile(double lnL) {

    std::string sumFile  = settings->getOutPutFileName() + ".sum";
    std::ofstream sumOut;
    sumOut.open( sumFile.c_str(), std::ios::out );
    
    sumOut << "Input file name                     = " << settings->getInputFileName() << std::endl;
    if (settings->getIsReversible() == true)
        sumOut << "Time reversible                     = YES" << std::endl;
    else
        sumOut << "Time reversible                     = NO" << std::endl;
    sumOut << "Number of gamma rate categories     = " << settings->getNumGammaCats() << std::endl;
    sumOut << "Branch length exponential parameter = " << settings->getBrlenLambda() << std::endl;
    sumOut << "Gamma shape exponential parameter   = " << settings->getAsrvLambda() << std::endl;
    sumOut << "Log of the marginal likelihood      = " << lnL << std::endl;

    sumOut.close();
}

void Mcmc::run(void) {

    std::cout << "   * Running Markov chain Monte Carlo analysis" << std::endl;

    // parameters of run
    int chainLength      = settings->getChainLength();
    int printFrequency   = settings->getPrintFrequency();
    int sampleFrequency  = settings->getSampleFrequency();
    std::string parmFile = settings->getOutPutFileName() + ".p";
    std::string treeFile = settings->getOutPutFileName() + ".t";

    // output to screen
    std::cout << "   * Number of MCMC interations = " << chainLength << std::endl;

    // open files for logging
    parmOut.open( parmFile.c_str(), std::ios::out );
    treeOut.open( treeFile.c_str(), std::ios::out );

    // do some initialization
    model->getTree(1)->updateAllConditionalLikelihoods(true);
    double curLnL = model->lnLikelihood();
    double curLnPrior = model->lnPrior();
    std::cout << "   * Initial log likelihood = " << curLnL << std::endl;

    // get information on moves
    std::vector<Move*>& moves = model->getMoves();
    std::set<Move*> uniqueMoves;
    for (int i=0; i<moves.size(); i++)
        uniqueMoves.insert(moves[i]);

    sample(true);
    auto t1 = std::chrono::high_resolution_clock::now();
    
            for (Move* m : uniqueMoves)
                m->clear();

    for (int n=1; n<=chainLength; n++)
        {
        // cycle through all the moves
        for (Move* m : moves)
            {
            // update the parameter
            double lnProposalRatio = m->update();
            
            // calculate the acceptance probability
            double newLnL = model->lnLikelihood();
            double newLnPrior = model->lnPrior();
            double R = acceptanceProb( (newLnL-curLnL) + (newLnPrior-curLnPrior) + lnProposalRatio );
            
            // accept or reject the proposal
            if (rv->uniformRv() < R)
                {
                m->accept();
                curLnL = newLnL;
                curLnPrior = newLnPrior;
                }
            else
                {
                m->reject();
                }
            }

        if (n % printFrequency == 0)
            {
            auto tNow = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = tNow - t1;
            double timePerIteration = elapsed.count() / n;
            double timeRemaining = (chainLength - n) * timePerIteration;
            std::string tr = formatTimeRemaining(timeRemaining);

            std::string rootString = model->getRootBipartition();
            std::cout << "   * " << std::setw(6) << n << " -- " << std::fixed << std::setprecision(2) << curLnL << " -- " << rootString << " (Estimated time remaining: " << tr << ")" << std::endl;
            }
            
         if (n % sampleFrequency == 0)
            {
            rootSamples->addBipartition( model->getRootBipartition() );
            sample(n, curLnL);
            }
            
        if (n % 1000 == 0)
            rootSamples->print();

        }
    sample(false);

    // print acceptance rates for stone
    for (Move* m : uniqueMoves)
        std::cout << "   * Accepted " << m->acceptanceProbability() * 100.0 << "\% of proposals for updates to " << m->getName() << std::endl;
    
    // result time
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = t2 - t1;
    std::cout << "   * MCMC took " << formatTimeRemaining(elapsed.count()) << std::endl;

    // close files
    parmOut.close();
    treeOut.close();
}

void Mcmc::sample(bool isFirstLine) {

    if (isFirstLine == true)
        {
        // parameter file
        std::vector<Move*> moves = model->getMoves();
        parmOut << "Generation" << '\t';
        for (int i=0; i<moves.size(); i++)
            {
            std::string mn = moves[i]->getName();
            parmOut << mn << '\t';
            }
        parmOut << std::endl;

        // tree file
        ParameterTree* t = model->getTree(0);
        treeOut << "#NEXUS" << std::endl << std::endl;
        treeOut << "begin trees;" << std::endl;
        treeOut << "   translate" << std::endl;
        std::vector<std::string> taxonNames = t->getTaxonNames();
        for (int i=0; i<taxonNames.size(); i++)
            {
            treeOut << "   " << i+1 << " " << taxonNames[i];
            if (i + 1 < taxonNames.size())
                treeOut << "," << std::endl;
            else
                treeOut << ";" << std::endl;
            }
        }
    else
        {
        treeOut << "end;" << std::endl;
        }
}

void Mcmc::sample(int gen, double lnL) {

    // parameter file
    std::vector<Move*> moves = model->getMoves();
    parmOut << gen << '\t';
    for (int i=0; i<moves.size(); i++)
        {
        std::vector<double> vals = moves[i]->values();
        for (int j=0; j<vals.size(); j++)
            parmOut << std::fixed << std::setprecision(6) << vals[j] << '\t';
        }
    
    EigenSystem& eigs = EigenSystem::eigenSystem();
    std::vector<double> pi = eigs.getStationaryFrequencies();
    for (int i=0; i<4; i++)
        parmOut << std::fixed << std::setprecision(6) << pi[i] << '\t';
    parmOut << std::endl;

    // tree file
    ParameterTree* t = model->getTree(0);
    std::string newickString = t->getNewick();
    treeOut << "   tree " << "sample_" << gen << " = " << newickString << ";" << std::endl;
}
