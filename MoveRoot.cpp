#include <cmath>
#include <iostream>
#include "Model.h"
#include "MoveRoot.h"
#include "Msg.h"
#include "Node.h"
#include "ParameterTree.h"
#include "RandomVariable.h"
#include "TransitionProbabilities.h"

#undef DEBUG_ROOT_MOVE


MoveRoot::MoveRoot(RandomVariable* r, Model* m, std::string nm, ParameterTree* t[2]) : Move(r, m, nm) {

    myParameter[0] = t[0];
    myParameter[1] = t[1];
    tuning = log(4.0);
}

void MoveRoot::accept(void) {
    
    numTries++;
    numAccepted++;

    (*myParameter[0]) = (*myParameter[1]);
}

void MoveRoot::reject(void) {

    numTries++;

    (*myParameter[1]) = (*myParameter[0]);
    myParameter[1]->updateAllTransitionProbabilities(true);
    model->updateTransitionProbabilities();
}

void MoveRoot::restore(void) {

}

void MoveRoot::tune(void) {

}

double MoveRoot::update(void) {

    ParameterTree* t = myParameter[1];
    Node* r = t->getRoot();
    std::vector<Node*> originalRootDescendants = r->getDescendants();
    if (originalRootDescendants.size() != 2)
        Msg::error("Root should have two descendants");

    // find the area of arrangement at the root of the tree
    Node* potentialU[2];
    int numPotentialU = 0;
    if (originalRootDescendants[0]->getIsLeaf() == false)
        potentialU[numPotentialU++] = originalRootDescendants[0];
    if (originalRootDescendants[1]->getIsLeaf() == false)
        potentialU[numPotentialU++] = originalRootDescendants[1];
    Node* u = potentialU[(int)(numPotentialU * rv->uniformRv())];
    Node* v = u->getAncestor();
    std::vector<Node*> uDescendants = u->getDescendants();
    if (uDescendants.size() != 2)
        {
        t->print();
        for (int i=0; i<originalRootDescendants.size(); i++)
            std::cout << "originalRootDescendants[" << i << "] = " << originalRootDescendants[i]->getIndex() << " (" << originalRootDescendants[i] << ")" << std::endl;
        std::cout << "potentialU[0] = " << potentialU[0] << std::endl;
        std::cout << "potentialU[1] = " << potentialU[1] << std::endl;
        std::cout << "numPotentialU = " << numPotentialU << std::endl;
        std::cout << "u = " << u->getIndex() << std::endl;
        Msg::error("Node u should have two descendants");
        }
    Node* a = uDescendants[0];
    Node* b = uDescendants[1];
    Node* c = NULL;
    if (u == originalRootDescendants[0])
        c = originalRootDescendants[1];
    else
        c = originalRootDescendants[0];
    
    // get lengths of branches
    double aLength = a->getMyBranch()->getLength();
    double bLength = b->getMyBranch()->getLength();
    double cLength = c->getMyBranch()->getLength();
    double uLength = u->getMyBranch()->getLength();

    // is reverse move forced
    bool isForcedBackMove = false;
    if (uLength < cLength)
        isForcedBackMove = false;
    else
        isForcedBackMove = true;

    double dAV = aLength + uLength;
    double dBV = bLength + uLength;
    double dCV = cLength;
    double oldH1 = 0.0;
    if (dAV < dBV && dAV < dCV)
        oldH1 = dAV;
    else if (dBV < dAV && dBV < dCV)
        oldH1 = dBV;
    else
        oldH1 = dCV;
    double newH1 = oldH1 * exp(tuning*(rv->uniformRv() - 0.5));
    dAV = dAV + (newH1 - oldH1);
    dBV = dBV + (newH1 - oldH1);
    dCV = dCV + (newH1 - oldH1);
    double h[3];
    h[0] = dAV;
    h[1] = dBV;
    h[2] = dCV;
    int dID[3];
    dID[0] = 0;
    dID[1] = 1;
    dID[2] = 2;
    if (h[0] > h[1])
        {
        double temp = h[0];
        h[0] = h[1];
        h[1] = temp;
        int dTemp = dID[0];
        dID[0] = dID[1];
        dID[1] = dTemp;
        }
    if (h[0] > h[2])
        {
        double temp = h[0];
        h[0] = h[2];
        h[2] = temp;
        int dTemp = dID[0];
        dID[0] = dID[2];
        dID[2] = dTemp;
        }
    if (h[1] > h[2])
        {
        double temp = h[1];
        h[1] = h[2];
        h[2] = temp;
        int dTemp = dID[1];
        dID[1] = dID[2];
        dID[2] = dTemp;
        }

    double newDAV = 0.0, newDBV = 0.0, newDCV = 0.0;
    for (int i=0; i<3; i++)
        {
        if (dID[i] == 0)
            newDAV = h[i];
        else if (dID[i] == 1)
            newDBV = h[i];
        else
            newDCV = h[i];
        }

    double x = rv->uniformRv() * h[1];
    
#   if defined (DEBUG_ROOT_MOVE)
    printf ("h[0] = %lf h[1] = %lf h[2] = %lf (%d %d %d)\n", h[0], h[1], h[2], dID[0], dID[1], dID[2]);
    printf ("x = %lf\n", x);
    t->print("before");
#   endif

    double lnProposalRatio = log(newH1) - log(oldH1);
    
    bool isForcedForwardMove = false;
    if (x > h[0])
        {
        isForcedForwardMove = true;
        if (dID[0] == 0)
            {
            /* a is shortest, (b,c) forced */
            u->removeNeighbor(a);
            u->addNeighbor(c);
            v->removeNeighbor(c);
            v->addNeighbor(a);
            a->removeNeighbor(u);
            a->addNeighbor(v);
            a->setAncestor(v);
            c->removeNeighbor(v);
            c->addNeighbor(u);
            c->setAncestor(u);
            
            a->getMyBranch()->setEnds( a, v );
            c->getMyBranch()->setEnds( c, u );
            
            /*u->left = b;
            u->right = c;
            u->anc = v;
            a->anc = v;
            b->anc = u;
            c->anc = u;
            v->left = a;
            v->right = u;*/
            
            u->getMyBranch()->setLength(x);
            a->getMyBranch()->setLength(dAV);
            b->getMyBranch()->setLength(dBV - x);
            c->getMyBranch()->setLength(dCV - x);
            }
        else if (dID[0] == 1)
            {
            /* b is shortest, (a,c) forced */
            u->removeNeighbor(b);
            u->addNeighbor(c);
            v->removeNeighbor(c);
            v->addNeighbor(b);
            b->removeNeighbor(u);
            b->addNeighbor(v);
            b->setAncestor(v);
            c->removeNeighbor(v);
            c->addNeighbor(u);
            c->setAncestor(u);

            b->getMyBranch()->setEnds( b, v );
            c->getMyBranch()->setEnds( c, u );

            /*u->left = a;
            u->right = c;
            u->anc = v;
            a->anc = u;
            b->anc = v;
            c->anc = u;
            v->left = b;
            v->right = u;*/
            
            u->getMyBranch()->setLength(x);
            a->getMyBranch()->setLength(dAV - x);
            b->getMyBranch()->setLength(dBV);
            c->getMyBranch()->setLength(dCV - x);
            }
        else
            {
            /* c is shortest, (a,b) forced */
            
            // no need to update neighbors as topology doesn't change
            
            /*u->left = a;
            u->right = b;
            u->anc = v;
            a->anc = u;
            b->anc = u;
            c->anc = v;
            v->left = c;
            v->right = u;*/
            
            u->getMyBranch()->setLength(x);
            a->getMyBranch()->setLength(dAV - x);
            b->getMyBranch()->setLength(dBV - x);
            c->getMyBranch()->setLength(dCV);
            }
        }
    else
        {
        isForcedForwardMove = false;
        double ran = rv->uniformRv();
        if (ran <= 0.33333)
            {
            /* a with c */
            u->removeNeighbor(b);
            u->addNeighbor(c);
            v->removeNeighbor(c);
            v->addNeighbor(b);
            b->removeNeighbor(u);
            b->addNeighbor(v);
            b->setAncestor(v);
            c->removeNeighbor(v);
            c->addNeighbor(u);
            c->setAncestor(u);

            b->getMyBranch()->setEnds( b, v );
            c->getMyBranch()->setEnds( c, u );

            /*u->right = c;
            b->anc = v;
            c->anc = u;
            if (v->left == u)
                v->right = b;
            else
                v->left = b;*/
                
            u->getMyBranch()->setLength(x);
            a->getMyBranch()->setLength(newDAV - x);
            b->getMyBranch()->setLength(newDBV);
            c->getMyBranch()->setLength(newDCV - x);
            }
        else if (ran > 0.33333 && ran <= 0.66666)
            {
            /* b with c */
            u->removeNeighbor(a);
            u->addNeighbor(c);
            v->removeNeighbor(c);
            v->addNeighbor(a);
            a->removeNeighbor(u);
            a->addNeighbor(v);
            a->setAncestor(v);
            c->removeNeighbor(v);
            c->addNeighbor(u);
            c->setAncestor(u);

            a->getMyBranch()->setEnds( a, v );
            c->getMyBranch()->setEnds( c, u );

            /*u->left = c;
            a->anc = v;
            c->anc = u;
            if (v->left == u)
                v->right = a;
            else
                v->left = a;*/
                
            //a->length += u->length;  // check this!
            //c->length -= u->length;
            u->getMyBranch()->setLength(x);
            a->getMyBranch()->setLength(newDAV);
            b->getMyBranch()->setLength(newDBV - x);
            c->getMyBranch()->setLength(newDCV - x);
            }
        else
            {
            /* no change */
            u->getMyBranch()->setLength(x);
            a->getMyBranch()->setLength(newDAV - x);
            b->getMyBranch()->setLength(newDBV - x);
            c->getMyBranch()->setLength(newDCV);
            }
        }

#   if defined (DEBUG_ROOT_MOVE)
    t->print("after");
#   endif

    /* get hastings ratio here */
    if (isForcedForwardMove == false && isForcedBackMove == true)
        lnProposalRatio += log(3.0);
    else if (isForcedForwardMove == true && isForcedBackMove == false)
        lnProposalRatio += log(1.0/3.0);
    
    // get down pass sequence
    t->initializeDownPassSequence();
    
    // update the flags for conditional likelihoods
    t->flipAllActiveConditionalLikelihoods();
    t->updateAllConditionalLikelihoods(true);

    // update the transition probability for the branch
    t->updateAllTransitionProbabilities(true);
    model->updateTransitionProbabilities();

    return lnProposalRatio;
    
#   if 0
    // get information on original root
    ParameterTree* t = myParameter[1];
    Node* r = t->getRoot();
    std::vector<Node*> des = r->getDescendants();
    if (des.size() != 2)
        Msg::error("Root should have two descendants");
    Node* oldP0 = des[0];
    Node* oldP1 = des[1];
    Branch* b0 = t->findBranch(oldP0, r);
    Branch* b1 = t->findBranch(oldP1, r);
    if (b0 == NULL || b1 == NULL)
        Msg::error("Could not find branches at root");
    double oldV = b0->getLength() + b1->getLength();

    // find the new root branch
    Branch* newB = NULL;
    do
        {
        newB = t->randomBranch();
        } while(newB == b0 || newB == b1);
    Node* newP0 = newB->getEnd1();
    Node* newP1 = newB->getEnd2();
    double newV = newB->getLength();

    // remove the old root
    t->removeBranch(b0);
    t->removeBranch(b1);
    oldP0->removeNeighbor(r);
    oldP1->removeNeighbor(r);
    r->removeNeighbors();
    oldP0->addNeighbor(oldP1);
    oldP1->addNeighbor(oldP0);
    Branch* oldB = t->addBranch(oldP0, oldP1);
    oldB->setLength(oldV);

    // add the new root
    t->removeBranch(newB);
    newP0->removeNeighbor(newP1);
    newP1->removeNeighbor(newP0);
    newP0->addNeighbor(r);
    newP1->addNeighbor(r);
    r->addNeighbor(newP0);
    r->addNeighbor(newP1);
    Branch* newB0 = t->addBranch(newP0, r);
    Branch* newB1 = t->addBranch(newP1, r);
    double v0 = newV * rv->uniformRv();
    double v1 = newV - v0;
    newB0->setLength(v0);
    newB1->setLength(v1);
    
    // get down pass sequence
    t->initializeDownPassSequence();
    
    // update the flags for conditional likelihoods
    t->flipAllActiveConditionalLikelihoods();
    t->updateAllConditionalLikelihoods(true);

    // update the transition probability for the branch
    t->updateAllTransitionProbabilities(true);
    model->updateTransitionProbabilities();

    return log(1.0/oldV) - log(1.0/newV);
#   endif
}

std::vector<double> MoveRoot::values(void) {

    std::vector<double> vals;
    vals.push_back( 0.0 );
    return vals;
}
