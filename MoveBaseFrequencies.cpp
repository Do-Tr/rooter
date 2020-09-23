#include "EigenSystem.h"
#include "Model.h"
#include "MoveBaseFrequencies.h"
#include "ParameterBaseFrequencies.h"
#include "ParameterTree.h"
#include "RandomVariable.h"



MoveBaseFrequencies::MoveBaseFrequencies(RandomVariable* r, Model* m, std::string nm, ParameterBaseFrequencies* f[2]) : Move(r, m, nm) {

    myParameter[0] = f[0];
    myParameter[1] = f[1];
    oldValue.resize(4);
    tuning = 100.0;
}

void MoveBaseFrequencies::accept(void) {

    numTries++;
    numAccepted++;

    (*myParameter[0]) = (*myParameter[1]);
    ParameterTree* t0 = model->getTree(0);
    ParameterTree* t1 = model->getTree(1);
    (*t0) = (*t1);
}

void MoveBaseFrequencies::reject(void) {

    numTries++;

    (*myParameter[1]) = (*myParameter[0]);
    ParameterTree* t0 = model->getTree(0);
    ParameterTree* t1 = model->getTree(1);
    (*t1) = (*t0);
    EigenSystem& eigs = EigenSystem::eigenSystem();
    eigs.flipActiveValues();
}

void MoveBaseFrequencies::restore(void) {

}

void MoveBaseFrequencies::tune(void) {

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

double MoveBaseFrequencies::update(void) {

    // modify the base frequencies
    oldValue = myParameter[1]->getBaseFrequencies();
    std::vector<double> aForward(4);
    std::vector<double> aReverse(4);
    std::vector<double> oldFreqs(4);
    std::vector<double> newFreqs(4);
    for (int i=0; i<4; i++)
        {
        oldFreqs[i] = oldValue[i];
        aForward[i] = oldValue[i] * tuning;
        if (aForward[i] < 1E-100)
            return -(1E-310);
        }
    rv->dirichletRv(aForward, newFreqs);
    myParameter[1]->setBaseFrequencies(newFreqs);
    for (int i=0; i<4; i++)
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
    t->flipAllActiveTransitionProbabilities();
    t->updateAllTransitionProbabilities(true);
    model->updateTransitionProbabilities();

    // return the log of the Hastings ratio
    return rv->lnDirichletPdf(aReverse, oldFreqs) - rv->lnDirichletPdf(aForward, newFreqs);
}

std::vector<double> MoveBaseFrequencies::values(void) {

    std::vector<double> val = myParameter[1]->getBaseFrequencies();
    return val;
}
