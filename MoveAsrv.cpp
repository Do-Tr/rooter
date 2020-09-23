#include <cmath>
#include "Model.h"
#include "MoveAsrv.h"
#include "ParameterAsrv.h"
#include "ParameterTree.h"
#include "RandomVariable.h"



MoveAsrv::MoveAsrv(RandomVariable* r, Model* m, std::string nm, ParameterAsrv* gs[2]) : Move(r, m, nm) {

    myParameter[0] = gs[0];
    myParameter[1] = gs[1];
    tuning = log(2.0);
}

void MoveAsrv::accept(void) {

    numTries++;
    numAccepted++;

    (*myParameter[0]) = (*myParameter[1]);
    ParameterTree* t0 = model->getTree(0);
    ParameterTree* t1 = model->getTree(1);
    (*t0) = (*t1);
}

void MoveAsrv::reject(void) {

    numTries++;

    (*myParameter[1]) = (*myParameter[0]);
    ParameterTree* t0 = model->getTree(0);
    ParameterTree* t1 = model->getTree(1);
    (*t1) = (*t0);
}

void MoveAsrv::restore(void) {

}

void MoveAsrv::tune(void) {

    double targetRate = 0.44;
    double currentRate = acceptanceProbability();
    
    if ( currentRate > targetRate )
        tuning *= (1.0 + ((currentRate-targetRate)/(1.0 - targetRate)) );
    else
        tuning /= (2.0 - currentRate/targetRate);

    if (tuning < log(1.01) )
        tuning = log(1.01);
    else if (tuning > log(1000.0))
        tuning = log(1000.0);
}

double MoveAsrv::update(void) {

    // modify the gamma shape parameter
    oldValue = myParameter[1]->getAsrv();
    double randomFactor = exp( tuning * (rv->uniformRv()-0.5) );
    double newValue = oldValue * randomFactor;
    myParameter[1]->setAsrv(newValue);

    // update the flags for conditional likelihoods
    ParameterTree* t = model->getTree(1);
    t->flipAllActiveConditionalLikelihoods();
    t->updateAllConditionalLikelihoods(true);

    // update the transition probability for the branch
    t->updateAllTransitionProbabilities(true);
    t->flipAllActiveTransitionProbabilities();
    model->updateTransitionProbabilities();

    // return the log of the Hastings ratio
    return log(newValue) - log(oldValue);
}

std::vector<double> MoveAsrv::values(void) {

    std::vector<double> val;
    val.push_back(myParameter[1]->getAsrv());
    return val;
}
