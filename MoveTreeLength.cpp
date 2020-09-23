#include <iomanip>
#include <iostream>
#include "Branch.h"
#include "EigenSystem.h"
#include "Model.h"
#include "MoveTreeLength.h"
#include "Node.h"
#include "ParameterBaseFrequencies.h"
#include "ParameterTree.h"
#include "RandomVariable.h"



MoveTreeLength::MoveTreeLength(RandomVariable* r, Model* m, std::string nm, ParameterTree* t[2]) : Move(r, m, nm) {

    myParameter[0] = t[0];
    myParameter[1] = t[1];
    tuning = log(2.0);
}

void MoveTreeLength::accept(void) {

    numTries++;
    numAccepted++;

    (*myParameter[0]) = (*myParameter[1]);
}

void MoveTreeLength::reject(void) {

    numTries++;

    (*myParameter[1]) = (*myParameter[0]);
    myParameter[1]->updateAllTransitionProbabilities(true);
    model->updateTransitionProbabilities();
}

void MoveTreeLength::restore(void) {

}

void MoveTreeLength::tune(void) {

    double targetRate = 0.22;
    double currentRate = acceptanceProbability();
    
    if ( currentRate > targetRate )
        tuning /= (1.0 + ((currentRate-targetRate)/(1.0 - targetRate)) );
    else
        tuning *= (2.0 - currentRate/targetRate);
    
    if (tuning < 4.0)
        tuning = 4.0;
    else if (tuning > 3000.0)
        tuning = 3000.0;
}

double MoveTreeLength::update(void) {

    ParameterTree* t = myParameter[1];
    
    // get the current tree length
    double oldValue = t->getTreeLength();
    
    // propose a new tree length
    double randomFactor = exp( tuning * (rv->uniformRv()-0.5) );
    double newTreeLength = oldValue * randomFactor;
    
    //std::cout << std::fixed << std::setprecision(8) << oldValue << " -> " << newTreeLength << std::endl;

    // update all of the branch lengths
    std::vector<Node*> dpSeq = t->getDownPassSequence();
    for (int i=0; i<dpSeq.size(); i++)
        {
        Node* p = dpSeq[i];
        Branch* b = p->getMyBranch();
        if (b != NULL)
            {
            double v = b->getLength();
            b->setLength( v * randomFactor );
            }
        }

    // update the Eigen system
    EigenSystem& eigs = EigenSystem::eigenSystem();
    eigs.flipActiveValues();
    model->updateRateMatrix();

    // update the flags for conditional likelihoods
    t->flipAllActiveConditionalLikelihoods();
    t->updateAllConditionalLikelihoods(true);

    // update the transition probability for the branch
    t->flipAllActiveTransitionProbabilities();
    t->updateAllTransitionProbabilities(true);
    model->updateTransitionProbabilities();

    // return the log of the Hastings ratio
    return (dpSeq.size()-1) * log(randomFactor);
}

std::vector<double> MoveTreeLength::values(void) {

    std::vector<double> vals;
    vals.push_back( myParameter[1]->getTreeLength() );
    return vals;
}
