#include "EigenSystem.h"
#include "Model.h"
#include "MoveExchangabilityRates.h"
#include "ParameterExchangabilityRates.h"
#include "ParameterTree.h"
#include "RandomVariable.h"



MoveExchangabilityRates::MoveExchangabilityRates(RandomVariable* r, Model* m, std::string nm, ParameterExchangabilityRates* f[2], bool isTimeReversible) : Move(r, m, nm) {

    if (isTimeReversible == true)
        numRates = 6;
    else
        numRates = 12;
    myParameter[0] = f[0];
    myParameter[1] = f[1];
    oldValue.resize(numRates);
    tuning = 1000.0;
}

void MoveExchangabilityRates::accept(void) {

    numTries++;
    numAccepted++;

    (*myParameter[0]) = (*myParameter[1]);
    ParameterTree* t0 = model->getTree(0);
    ParameterTree* t1 = model->getTree(1);
    (*t0) = (*t1);
}

void MoveExchangabilityRates::reject(void) {

    numTries++;

    (*myParameter[1]) = (*myParameter[0]);
    ParameterTree* t0 = model->getTree(0);
    ParameterTree* t1 = model->getTree(1);
    (*t1) = (*t0);
    EigenSystem& eigs = EigenSystem::eigenSystem();
    eigs.flipActiveValues();
}

void MoveExchangabilityRates::restore(void) {

}

void MoveExchangabilityRates::tune(void) {

    double targetRate = 0.22;
    double currentRate = acceptanceProbability();
    
    if ( currentRate > targetRate )
        tuning /= (1.0 + ((currentRate-targetRate)/(1.0 - targetRate)) );
    else
        tuning *= (2.0 - currentRate/targetRate);

    if (tuning < numRates)
        tuning = numRates;
    else if (tuning > 5000.0)
        tuning = 5000.0;
}

double MoveExchangabilityRates::update(void) {

    // modify the exchangability rate parameters
    oldValue = myParameter[1]->getExchangabilityRates();
    std::vector<double> aForward(numRates);
    std::vector<double> aReverse(numRates);
    std::vector<double> oldFreqs(numRates);
    std::vector<double> newFreqs(numRates);
    for (int i=0; i<numRates; i++)
        {
        oldFreqs[i] = oldValue[i];
        aForward[i] = oldValue[i] * tuning;
        if (aForward[i] < 1E-100)
            return -(1E-310);
        }
    rv->dirichletRv(aForward, newFreqs);
    myParameter[1]->setExchangabilityRates(newFreqs);
    for (int i=0; i<numRates; i++)
        {
        aReverse[i] = newFreqs[i] * tuning;
        if (aReverse[i] < 1E-100)
            return -(1E-310);
        }

    // update the Eigen system
    EigenSystem& eigs = EigenSystem::eigenSystem();
    eigs.flipActiveValues();
    model->updateRateMatrix();

    // update the flags for conditional likelihoods
    ParameterTree* t = model->getTree(1);
    t->flipAllActiveConditionalLikelihoods();
    t->updateAllConditionalLikelihoods(true);

    // update the transition probability for the branch
    t->updateAllTransitionProbabilities(true);
    t->flipAllActiveTransitionProbabilities();
    model->updateTransitionProbabilities();

    // return the log of the Hastings ratio
    return rv->lnDirichletPdf(aReverse, oldFreqs) - rv->lnDirichletPdf(aForward, newFreqs);
}

std::vector<double> MoveExchangabilityRates::values(void) {

    std::vector<double> val = myParameter[1]->getExchangabilityRates();
    return val;
}
