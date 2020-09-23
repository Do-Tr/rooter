#include <cmath>
#include "Branch.h"
#include "EigenSystem.h"
#include "Model.h"
#include "MoveBranchLength.h"
#include "Msg.h"
#include "Node.h"
#include "ParameterTree.h"
#include "RandomVariable.h"
#include "TransitionProbabilities.h"



MoveBranchLength::MoveBranchLength(RandomVariable* r, Model* m, std::string nm, ParameterTree* t[2]) : Move(r, m, nm) {

    myParameter[0] = t[0];
    myParameter[1] = t[1];
    tuning = log(2.0);
}

void MoveBranchLength::accept(void) {

    numTries++;
    numAccepted++;

    //(*myParameter[0]) = (*myParameter[1]);
    // update the transition probability for the branch
    Branch* oldBranch0 = myParameter[0]->findBranch( myParameter[0]->findNodeIndexed(oldBranch1->getEnd1()->getIndex()), myParameter[0]->findNodeIndexed(oldBranch1->getEnd2()->getIndex()) );
    oldBranch0->flipActiveTi();
    
    oldBranch0->setLength( oldBranch1->getLength() );
    
    // update the flags for conditional likelihoods
    Node* p = oldBranch0->getAncestralNode();
    while (p != NULL)
        {
        p->flipActiveCl();
        p = p->getAncestor();
        }
}

void MoveBranchLength::reject(void) {

    numTries++;

    //(*myParameter[1]) = (*myParameter[0]);
    
    // update the transition probability for the branch
    oldBranch1->flipActiveTi();

    oldBranch1->setLength(oldValue);

    // update the flags for conditional likelihoods
    Node* p = oldBranch1->getAncestralNode();
    while (p != NULL)
        {
        p->flipActiveCl();
        p = p->getAncestor();
        }
}

void MoveBranchLength::restore(void) {

}

void MoveBranchLength::tune(void) {

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

double MoveBranchLength::update(void) {

    // find the branch
    ParameterTree* t = myParameter[1];
    Branch* b = t->randomBranch();
    if (b == NULL)
        Msg::error("Could not find branch");
    oldBranch1 = b;
    
    // modify the branch length
    oldValue = b->getLength();
    double randomFactor = exp( tuning * (rv->uniformRv()-0.5) );
    double newValue = oldValue * randomFactor;
    b->setLength(newValue);
    
    // update the transition probability for the branch
    t->updateAllTransitionProbabilities(false);
    b->setNeedsUpdate(true);
    b->flipActiveTi();
    model->updateTransitionProbabilities();
    
    // update the flags for conditional likelihoods
    Node* p = b->getAncestralNode();
    while (p != NULL)
        {
        p->flipActiveCl();
        p->setNeedsUpdate(true);
        p = p->getAncestor();
        }
    
    return log(randomFactor);
}

std::vector<double> MoveBranchLength::values(void) {

    std::vector<double> val;
    //val.push_back(myNode->getLength());
    return val;
}
