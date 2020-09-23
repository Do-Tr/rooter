#ifndef Model_H
#define Model_H

#include <string>
#include <vector>
class Alignment;
class ConditionalLikelihoods;
class Move;
class ParameterAsrv;
class ParameterBaseFrequencies;
class ParameterExchangabilityRates;
class ParameterTree;
class RandomVariable;
class RootBipartition;
class Settings;
class TipTimes;
class TransitionProbabilities;



class Model {

    public:
                                                Model(void) = delete;
                                                Model(RandomVariable* r, Settings* s, Alignment* a);
                                               ~Model(void);
        ParameterAsrv*                          getAsrv(void) { return asrv[1]; }
        std::vector<Move*>&                     getMoves(void) { return moves; }
        ParameterTree*                          getTree(int idx) { return tree[idx]; }
        std::string                             getRootBipartition(void);
        double                                  lnLikelihood(void);
        double                                  lnPrior(void);
        void                                    updateRateMatrix(void);
        void                                    updateTransitionProbabilities(void);

    protected:
        void                                    initializeConditionalLikelihoods(void);
        void                                    initializeMoves(void);
        void                                    initializeParameters(void);
        void                                    initiializeTransitionProbabilities(void);
        void                                    initializeTreeWithTipDateInformation(ParameterTree* t, TipTimes* tt);
        void                                    printParameters(void);
        void                                    printConditionalLikelihoods(void);
        void                                    printTransitionProbabilities(void);
    
        int                                     numGammaCats;
        int                                     numPatterns;
        bool                                    isTimeReversible;
        int*                                    numSitesOfPattern;
        Alignment*                              alignmentPtr;
        RandomVariable*                         rv;
        Settings*                               settingsPtr;
        TipTimes*                               tipTimes;
        RootBipartition*                        rootPartition;
    
        ParameterAsrv*                          asrv[2];
        ParameterBaseFrequencies*               pi[2];
        ParameterExchangabilityRates*           subrates[2];
        ParameterTree*                          tree[2];
        std::vector<Move*>                      moves;
    
        std::vector<ConditionalLikelihoods*>    condLikes[2];
        std::vector<TransitionProbabilities*>   transitionProbs[2];
};

#endif
