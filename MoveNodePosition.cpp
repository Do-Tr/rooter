#include <cmath>
#include <iostream>
#include "Branch.h"
#include "EigenSystem.h"
#include "Model.h"
#include "MoveNodePosition.h"
#include "Msg.h"
#include "Node.h"
#include "ParameterTree.h"
#include "RandomVariable.h"
#include "TransitionProbabilities.h"



MoveNodePosition::MoveNodePosition(RandomVariable* r, Model* m, std::string nm, ParameterTree* t[2]) : Move(r, m, nm) {

    myParameter[0] = t[0];
    myParameter[1] = t[1];
    tuning = 0.0;
}

void MoveNodePosition::accept(void) {

    numTries++;
    numAccepted++;

    //(*myParameter[0]) = (*myParameter[1]);
    // update the transition probability for the branch
    for (int i=0; i<10; i++)
        {
        Branch* b1 = oldBranch[i];
        if (b1 != NULL)
            {
            Branch* b0 = myParameter[0]->findBranch( myParameter[0]->findNodeIndexed(b1->getEnd1()->getIndex()), myParameter[0]->findNodeIndexed(b1->getEnd2()->getIndex()) );
            b0->flipActiveTi();
            b0->setLength( b1->getLength() );
            }
        }
    
    
    // update the flags for conditional likelihoods
    Branch* b0 = myParameter[0]->findBranch( myParameter[0]->findNodeIndexed(oldBranch[0]->getEnd1()->getIndex()), myParameter[0]->findNodeIndexed(oldBranch[0]->getEnd2()->getIndex()) );
    Node* p = b0->getDescendantNode();
    while (p != NULL)
        {
        p->flipActiveCl();
        p = p->getAncestor();
        }
}

void MoveNodePosition::reject(void) {

    numTries++;

    //(*myParameter[1]) = (*myParameter[0]);
    
    // update the transition probability for the branch
    for (int i=0; i<10; i++)
        {
        Branch* b1 = oldBranch[i];
        if (b1 != NULL)
            {
            b1->flipActiveTi();
            b1->setLength( oldValue[i] );
            }
        }

    // update the flags for conditional likelihoods
    Node* p = oldBranch[0]->getDescendantNode();
    while (p != NULL)
        {
        p->flipActiveCl();
        p = p->getAncestor();
        }
}

void MoveNodePosition::restore(void) {

}

void MoveNodePosition::tune(void) {

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

double MoveNodePosition::update(void) {

    // choose a branch
    ParameterTree* t = myParameter[1];
    Branch* b = t->randomInteriorBranch();
    if (b == NULL)
        Msg::error("Could not find branch");
    
    // remember the touched branches
    int numTouchedBranches = 0;
    for (int i=0; i<10; i++)
        oldBranch[i] = NULL;
    oldBranch[0] = b;
    oldValue[0] = b->getLength();
    numTouchedBranches++;
    
    // find the descendant branches, add them to the list of touched branches, find the smallest branch length
    std::vector<Branch*> db;
    Node* n = b->getDescendantNode();
    std::vector<Node*> descendants = n->getDescendants();
    if (descendants.size() == 0)
        {
        t->print();
        std::cout << "Problem node = " << b->getDescendantNode()->getIndex() << std::endl;
        Msg::error("No descendants of brnach chosen");
        }
    double minDescendantLength = -1.0;
    for (Node* p : descendants)
        {
        Branch* pb = p->getMyBranch();
        double v = pb->getLength();
        db.push_back( pb );
        oldBranch[numTouchedBranches] = pb;
        oldValue[numTouchedBranches++] = v;
        if (minDescendantLength < 0.0)
            minDescendantLength = v;
        else if (v < minDescendantLength)
            minDescendantLength = v;
        }

    // find the range of the change
    double h = minDescendantLength + b->getLength();
    double newV = rv->uniformRv() * h;
    double diff = oldBranch[0]->getLength() - newV;
   
    // modify the lengths of the branches
    b->setLength( b->getLength() - diff );
    for (Branch* pb : db)
        pb->setLength( pb->getLength() + diff );
    
    // check the branch lengths
    for (int i=0; i<numTouchedBranches; i++)
        {
        if (oldBranch[i]->getLength() < 0.0)
            Msg::error("Negative branch length produced");
        }
    
    // update the transition probability for the branch
    t->updateAllTransitionProbabilities(false);
    for (int i=0; i<numTouchedBranches; i++)
        {
        oldBranch[i]->setNeedsUpdate(true);
        oldBranch[i]->flipActiveTi();
        }
    model->updateTransitionProbabilities();
    
    // update the flags for conditional likelihoods
    Node* p = b->getDescendantNode();
    while (p != NULL)
        {
        p->flipActiveCl();
        p->setNeedsUpdate(true);
        p = p->getAncestor();
        }
    
#   if 0
    for (int i=0; i<10; i++)
        {
        Branch* b1 = oldBranch[i];
        if (b1 != NULL)
            std::cout << i << " -- " << b1->getLength() << " " << oldValue[i] << std::endl;
        }
#   endif
    
    return 0.0;
}

std::vector<double> MoveNodePosition::values(void) {

    std::vector<double> val;
    //val.push_back(myNode->getLength());
    return val;
}
