#include <iomanip>
#include <iostream>
#include <map>
#include "Alignment.h"
#include "ConditionalLikelihoods.h"
#include "EigenSystem.h"
#include "Model.h"
#include "MoveAsrv.h"
#include "MoveBaseFrequencies.h"
#include "MoveNodePosition.h"
#include "MoveExchangabilityRates.h"
#include "MoveRoot.h"
#include "MoveTreeLength.h"
#include "Msg.h"
#include "Node.h"
#include "Parameter.h"
#include "ParameterAsrv.h"
#include "ParameterBaseFrequencies.h"
#include "ParameterExchangabilityRates.h"
#include "ParameterTree.h"
#include "RandomVariable.h"
#include "RootBipartition.h"
#include "Settings.h"
#include "TipTimes.h"
#include "TransitionProbabilities.h"



Model::Model(RandomVariable* r, Settings* s, Alignment* a) {

    // remember important objects
    rv = r;
    settingsPtr = s;
    alignmentPtr = a;
    
    // initialize variables
    numPatterns = alignmentPtr->getNumPatterns();
    numGammaCats = settingsPtr->getNumGammaCats();
    isTimeReversible = settingsPtr->getIsReversible();
    
    // get tip times
    tipTimes = new TipTimes(alignmentPtr, settingsPtr);
    
    // initialize parameters
    initializeParameters();
    
    // initialize the conditional likelihoods
    initializeConditionalLikelihoods();
    
    // initialize the eigen system and transition probabilities
    updateRateMatrix();
    initiializeTransitionProbabilities();
        
    // initialize the moves
    initializeMoves();

}

Model::~Model(void) {

    for (int i=0; i<2; i++)
        {
        delete tree[i];
        delete subrates[i];
        if (asrv[i] != NULL)
            delete asrv[i];
        if (pi[i] != NULL)
            delete pi[i];

        for (int n=0; n<condLikes[i].size(); n++)
            delete condLikes[i][n];
        for (int n=0; n<transitionProbs[i].size(); n++)
            delete transitionProbs[i][n];
        }
    delete [] numSitesOfPattern;
    
    std::set<Move*> deletedMoves;
    for (int i=0; i<moves.size(); i++)
        {
        Move* m = moves[i];
        std::set<Move*>::iterator it = deletedMoves.find(m);
        if (it == deletedMoves.end())
            {
            deletedMoves.insert(m);
            delete m;
            }
        }
    
    delete tipTimes;
    delete rootPartition;
}

std::string Model::getRootBipartition(void) {

    std::string s = rootPartition->getRootBipartition(tree[0]);
    if (s == "blank")
        {
        tree[0]->print();
        rootPartition->print();
        }
    return s;
}

void Model::initializeConditionalLikelihoods(void) {

    int numNodes = (int)tree[0]->getNumNodes();
    std::vector<std::string> tn = alignmentPtr->getTaxonNames();

    for (int i=0; i<2; i++)
        {
        for (int n=0; n<numNodes; n++)
            condLikes[i].push_back( new ConditionalLikelihoods(alignmentPtr, numGammaCats) );
        
        for (int n=0; n<tn.size(); n++)
            condLikes[i][n]->initializeTipConditonalLikelihoods(alignmentPtr, tn[n]);
        }

    numSitesOfPattern = new int[numPatterns];
    for (int i=0; i<numPatterns; i++)
        numSitesOfPattern[i] = alignmentPtr->getNumSitesOfPattern(i);

    //printConditionalLikelihoods();
}

void Model::initializeMoves(void) {

    std::cout << "   * Initializing proposal for exchangability rates" << std::endl;
    moves.push_back( new MoveExchangabilityRates(rv, this, "exchangability rates", subrates, isTimeReversible) );
    
    if (numGammaCats > 1)
        {
        std::cout << "   * Initializing proposal for gamma shape parameter" << std::endl;
        moves.push_back( new MoveAsrv(rv, this, "gamma shape parameter", asrv) );
        }
    
    if (isTimeReversible == true)
        {
        std::cout << "   * Initializing proposal for base frequencies" << std::endl;
        moves.push_back( new MoveBaseFrequencies(rv, this, "base frequencies", pi) );
        }

    std::cout << "   * Initializing proposals for branch lengths" << std::endl;
    std::vector<Node*>& dp = tree[1]->getDownPassSequence();
    MoveNodePosition* npm = new MoveNodePosition(rv, this, "node position", tree);
    for (int i=0; i<dp.size(); i++)
        moves.push_back( npm );
    
    moves.push_back( new MoveTreeLength(rv, this, "tree length", tree) );

    std::cout << "   * Initializing proposal for root" << std::endl;
    moves.push_back( new MoveRoot(rv, this, "root", tree) );
}

void Model::initializeParameters(void) {

    std::cout << "   * Initializing tree from file \"" << settingsPtr->getTreeFileName() << "\"" << std::endl;
    tree[0] = new ParameterTree(rv, this, "Tree", alignmentPtr->getTaxonNames(), settingsPtr->getTreeFileName(), settingsPtr->getBrlenLambda());
    rootPartition = new RootBipartition();
    initializeTreeWithTipDateInformation(tree[0], tipTimes);
    tree[1] = new ParameterTree(*tree[0]);
    
    std::cout << "   * Initializing exchangability rate parameter" << std::endl;
    subrates[0] = new ParameterExchangabilityRates(rv, this, "Exchangability Rates", settingsPtr->getIsReversible());
    subrates[1] = new ParameterExchangabilityRates(*subrates[0]);
    
    if (settingsPtr->getNumGammaCats() > 1)
        {
        std::cout << "   * Initializing gamma shape parameter" << std::endl;
        asrv[0] = new ParameterAsrv(rv, this, "Gamma Shape Parameter", settingsPtr->getNumGammaCats(), settingsPtr->getAsrvLambda());
        asrv[1] = new ParameterAsrv(*asrv[0]);
        }
    else
        {
        std::cout << "   * Equal rates among sites" << std::endl;
        asrv[0] = NULL;
        asrv[1] = NULL;
        }
    
    if (isTimeReversible == true)
        {
        std::cout << "   * Initializing base frequencies parameter" << std::endl;
        pi[0] = new ParameterBaseFrequencies(rv, this, "Base Frequencies");
        pi[1] = new ParameterBaseFrequencies(*pi[0]);
        }
    else
        {
        std::cout << "   * Base frequencies determined from rate matrix" << std::endl;
        pi[0] = NULL;
        pi[1] = NULL;
        }

    //printParameters();
}

void Model::initiializeTransitionProbabilities(void) {

    int numNodes = (int)tree[0]->getNumNodes();
    for (int i=0; i<2; i++)
        {
        for (int n=0; n<numNodes; n++)
            transitionProbs[i].push_back( new TransitionProbabilities(numGammaCats, this) );
        }
        
    tree[1]->updateAllTransitionProbabilities(true);
    updateTransitionProbabilities();
    
    //printTransitionProbabilities();
}

void Model::initializeTreeWithTipDateInformation(ParameterTree* t, TipTimes* tt) {

    if (t == NULL || tt == NULL)
        {
        Msg::warning("   * Warning: Cannot initialize tree using tip date information");
        return;
        }
    std::cout << "   * Initializing tree branch lenghts using tip date information" << std::endl;
    
    // initialize the tip times
    double maxTipDepth = 0.0;
    std::map<std::string, double> tipDates = tt->getTipTimesMap();
    for (std::map<std::string, double>::iterator it=tipDates.begin(); it != tipDates.end(); it++)
        {
        Node* tipNode = t->findNodeWithTaxonName(it->first);
        if (tipNode == NULL)
            Msg::error("Could not find node " + it->first + " when initializing tip dates");
        tipNode->setScratchVariable(it->second);
        if (it->second > maxTipDepth)
            maxTipDepth = it->second;
        }

    // get the depth increment
    double increment = 0.0;
    if (maxTipDepth < 0.000000001)
        increment = 1.0;
    else
        increment = maxTipDepth / 10.0;
    
    // initialize interior node times
    std::vector<Node*> dpSeq = t->getDownPassSequence();
    for (int i=0; i<dpSeq.size(); i++)
        {
        Node* p = dpSeq[i];
        if (p->getIsLeaf() == false)
            {
            // get maximum time of descendant
            double maxTime = 0.0;
            std::vector<Node*> desc = p->getDescendants();
            for (int j=0; j<desc.size(); j++)
                {
                if (desc[j]->getScratchVariable() > maxTime)
                    maxTime = desc[j]->getScratchVariable();
                }
                
            // set the time for this node
            p->setScratchVariable(maxTime + increment);
            }
        }

    // initialize the tree height
    double treeHeight = rv->exponentialRv( settingsPtr->getBrlenLambda() );
    double factor = treeHeight / t->getRoot()->getScratchVariable();
    for (int i=0; i<dpSeq.size(); i++)
        {
        Node* p = dpSeq[i];
        if (p->getAncestor() != NULL)
            {
            Branch* b = p->getMyBranch();
            double v = factor * (p->getAncestor()->getScratchVariable() - p->getScratchVariable());
            b->setLength(v);
            }
        }

    //t->print();
}

double Model::lnLikelihood(void) {

    ParameterTree* t = tree[1];

    // conditional likelihoods down the tree
    std::vector<Node*>& dp = t->getDownPassSequence();
    for (int n=0; n<dp.size(); n++)
        {
        Node* p = dp[n];
        if (p->getIsLeaf() == false && p->needsUpdate() == true)
            {            
            std::vector<Node*> des = p->getDescendants();
            Node* ndeL = des[0];
            Node* ndeR = des[1];
            
            ConditionalLikelihoods* clInfoL = condLikes[ ndeL->getActiveCl() ][ ndeL->getIndex() ];
            ConditionalLikelihoods* clInfoR = condLikes[ ndeR->getActiveCl() ][ ndeR->getIndex() ];
            ConditionalLikelihoods* clInfoP = condLikes[    p->getActiveCl() ][    p->getIndex() ];
            double* clL = clInfoL->getCondLike();
            double* clR = clInfoR->getCondLike();
            double* clP = clInfoP->getCondLike();
            double* scL = clInfoL->getLnScalerDp();
            double* scR = clInfoR->getLnScalerDp();
            double* scP = clInfoP->getLnScalerDp();
            double* lnS = clInfoP->getLnScaler();
            
            Branch* bL = t->findBranch(ndeL, p);
            Branch* bR = t->findBranch(ndeR, p);
            double*** pL = transitionProbs[ bL->getActiveTi() ][ ndeL->getIndex() ]->getTiProbs();
            double*** pR = transitionProbs[ bR->getActiveTi() ][ ndeR->getIndex() ]->getTiProbs();

            for (int c=0; c<numPatterns; c++)
                {
                double* clPStart = clP;
                double maxCl = 0.0;
                for (int k=0; k<numGammaCats; k++)
                    {
                    for (int i=0; i<4; i++)
                        {
                        double sumL = 0.0, sumR = 0.0;
                        for (int j=0; j<4; j++)
                            {
                            sumL += pL[k][i][j] * clL[j];
                            sumR += pR[k][i][j] * clR[j];
                            }
                        clP[i] = sumL * sumR;
                        if (clP[i] > maxCl)
                            maxCl = clP[i];
                        }
                    clP += 4;
                    clL += 4;
                    clR += 4;
                    }
                double scaleFactor = 1.0 / maxCl;
                for (int i=0; i<numGammaCats*4; i++)
                    clPStart[i] *= scaleFactor;
                lnS[c] = log(maxCl);
                scP[c] = scL[c] + scR[c] + lnS[c];
                }
            
            p->setNeedsUpdate(false);
            }
        }
    
    // calculate the likelihood, using the conditional likelihoods at the root of the tree
    double gammaWeight = 1.0 / numGammaCats;
    ConditionalLikelihoods* clInfoR = condLikes[ t->getRoot()->getActiveCl() ][ t->getRoot()->getIndex() ];
    double* clR = clInfoR->getCondLike();
    double* scR = clInfoR->getLnScalerDp();
    double lnL = 0.0;
    EigenSystem& eigs = EigenSystem::eigenSystem();
    std::vector<double>& pi = eigs.getStationaryFrequencies();
    for (int c=0; c<numPatterns; c++)
        {
        double like = 0.0;
        for (int k=0; k<numGammaCats; k++)
            {
            for (int i=0; i<4; i++)
                like += clR[i] * pi[i] * gammaWeight;
            clR += 4;
            }
        lnL += numSitesOfPattern[c] * (log(like) + scR[c]);
        //lnL += numSitesOfPattern[c] * log(like);
        }
    
    return lnL;
}

double Model::lnPrior(void) {

    double lnP = 0.0;
    
    lnP += tree[1]->lnPriorProb();
    lnP += subrates[1]->lnPriorProb();
    lnP += asrv[1]->lnPriorProb();
    if (isTimeReversible == true)
        lnP += pi[1]->lnPriorProb();

    return lnP;
}

void Model::printParameters(void) {

    tree[0]->print("   * Tree[0]:");
    tree[1]->print("   * Tree[1]:");
    subrates[0]->print();
    subrates[1]->print();
    if (asrv[0] != NULL)
        {
        asrv[0]->print();
        asrv[1]->print();
        }
    if (pi[0] != NULL)
        {
        pi[0]->print();
        pi[1]->print();
        }
}

void Model::printConditionalLikelihoods(void) {

    int whichSpace = 0;
    int ntax = alignmentPtr->getNumTaxa();
    int npat = alignmentPtr->getNumPatterns();
    int stride = 7;
    int numStrides = npat / stride;
    if (npat % stride != 0)
        numStrides++;
    
    for (int s=0; s<numStrides; s++)
        {
        int beg = s * stride;
        int end = beg + stride;
        if (end >= npat)
            end = npat;
            
        std::cout << "         ";
        for (int j=beg; j<end; j++)
            {
            std::cout << j+1;
            int nspaces = 5*numGammaCats+2;
            int ndigits = floor(log10(j+1)+1);
            nspaces -= ndigits;
            for (int k=0; k<nspaces; k++)
                std::cout << " ";
            }
        std::cout << std::endl;
        
        for (int i=0; i<condLikes[whichSpace].size(); i++)
            {
            std::cout << std::setw(5) << i << " -- ";
            double* x = condLikes[whichSpace][i]->getCondLike();
            x += beg * 4 * numGammaCats;
            for (int j=beg; j<end; j++)
                {
                char c;
                if (i < ntax)
                    c = alignmentPtr->getNucleotideCharForPattern(i, j);
                else c = '*';
                std::cout << c << " ";
                for (int k=0; k<numGammaCats; k++)
                    {
                    for (int m=0; m<4; m++)
                        std::cout << std::fixed << std::setprecision(0) << x[m];
                    std::cout << " ";
                    x += 4;
                    }
                }
            std::cout << std::endl;
            }
        }
}

void Model::printTransitionProbabilities(void) {

    int whichSpace = 0;
    for (int i=0; i<transitionProbs[whichSpace].size(); i++)
        {
        std::cout << "Transition Probabilities for index " << i << ":" << std::endl;
        transitionProbs[whichSpace][i]->print();
        }
}

void Model::updateRateMatrix(void) {

    //NucleotideSquareMatrix_t Q;
    NucleotideSquareMatrix_t Q;
    
    // fill in off diagonal components of rate matrix
    std::vector<double>& sr = subrates[1]->getExchangabilityRates();
    if (isTimeReversible == false)
        {
        // non-reversible model
        int k = 0;
        for (int i=0; i<4; i++)
            {
            for (int j=0; j<4; j++)
                {
                if (i != j)
                    Q(i,j) = sr[k++];
                }
            }
        }
    else
        {
        // gtr model
        std::vector<double>& bf = pi[1]->getBaseFrequencies();
        int k = 0;
        for (int i=0; i<4; i++)
            {
            for (int j=i+1; j<4; j++)
                {
                Q(i,j) = sr[k] * bf[j];
                Q(j,i) = sr[k] * bf[i];
                k++;
                }
            }
        }
    
    // fill in the diagonal elements of the rate matrix
    for (int i=0; i<4; i++)
        {
        double sum = 0.0;
        for (int j=0; j<4; j++)
            {
            if (i != j)
                sum += Q(i,j);
            }
        Q(i,i) = -sum;
        }
    
    // calculate the stationary frequencies
    EigenSystem& eigs = EigenSystem::eigenSystem();
    if (isTimeReversible == false)
        eigs.setStationaryFrequencies( eigs.calulateStationaryFrequencies(Q) );
    else
        eigs.setStationaryFrequencies( pi[1]->getBaseFrequencies() );

    // rescale the rate matrix
    std::vector<double>& pi = eigs.getStationaryFrequencies();
    double averageRate = 0.0;
    for (int i=0; i<4; i++)
        averageRate += -pi[i] * Q(i,i);
    double scaleFactor = 1.0 / averageRate;
    Q *= scaleFactor;
    
    // calculate the Eigenvalues and Eigenvectors and do some precomputation
    eigs.calculateEigenSystem(Q);
}

void Model::updateTransitionProbabilities(void) {

    ParameterTree* t = tree[1];
    std::vector<Node*>& dp = t->getDownPassSequence();
    for (int n=0; n<dp.size(); n++)
        {
        Node* p = dp[n];
        if (p->getAncestor() != NULL)
            {
            Branch* b = t->findBranch(p, p->getAncestor());
            if (b->needsUpdate() == true)
                {
                TransitionProbabilities* tp = transitionProbs[b->getActiveTi()][p->getIndex()];
                tp->tiProbs( b->getLength() );
                b->setNeedsUpdate(false);
                }
            }
        }
}
