#include <iomanip>
#include <iostream>
#include "Alignment.h"
#include "ConditionalLikelihoods.h"
#include "Msg.h"


ConditionalLikelihoods::ConditionalLikelihoods(Alignment* a, int nc) {

    // initialize some instance variables
    numGammaCats   = nc;
    numPatterns    = a->getNumPatterns();
    
    // dynamically allocate the information for the conditional likelihoods
    int clsSize = numPatterns * 4 * numGammaCats;
    cls = new double[clsSize];
    for (int i=0; i<clsSize; i++)
        cls[i] = 0.0;

    // allocate the log scaler vector
    lnScaler = new double[numPatterns];
    lnScalerDp = new double[numPatterns];
    for (int i=0; i<numPatterns; i++)
        {
        lnScaler[i] = 0.0;
        lnScalerDp[i] = 0.0;
        }
}

ConditionalLikelihoods::~ConditionalLikelihoods(void) {

    delete [] cls;
    delete [] lnScaler;
    delete [] lnScalerDp;
}

void ConditionalLikelihoods::initializeTipConditonalLikelihoods(Alignment* a, std::string tName) {

    int idx = a->getTaxonIndex(tName);
    if (idx == -1)
        Msg::error("Couldn't find taxon " + tName + " in alignment");
    
    double* p = &cls[0];
    for (int c=0; c<numPatterns; c++)
        {
        int nucCode = a->getNucleotideForPattern(idx, c);
        int possibleNucs[4];
        a->getPossibleNucs(nucCode, possibleNucs);
        for (int k=0; k<numGammaCats; k++)
            {
            for (int s=0; s<4; s++)
                p[s] = (double)possibleNucs[s];
            p += 4;
            }
        }
}

void ConditionalLikelihoods::print(void) {

    double* p = &cls[0];
    for (int c=0; c<numPatterns; c++)
        {
        std::cout << std::fixed << std::setprecision(0);
        for (int s=0; s<4; s++)
            std::cout << p[s];
        std::cout << " ";
        p += 4 * numGammaCats;
        }
    std::cout << std::endl;
}
