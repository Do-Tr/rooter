#ifndef TransitionProbabilities_H
#define TransitionProbabilities_H

class Model;


class TransitionProbabilities {

    public:
                                    TransitionProbabilities(void) = delete;
                                    TransitionProbabilities(int k, Model* m);
                                   ~TransitionProbabilities(void);
        void                        print(void);
        void                        tiProbs(double v);
        double***                   getTiProbs(void) { return p; }

    protected:
        Model*                      myModel;
        double***                   p;
        int                         numStates;
        int                         numGammaCats;
};

#endif
