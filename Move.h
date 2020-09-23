#ifndef Move_H
#define Move_H

#include <string>
#include <vector>
class Model;
class RandomVariable;



class Move {

    public:
                                    Move(void) = delete;
                                    Move(RandomVariable* r, Model* m, std::string nm);
        virtual                    ~Move(void) { }
        virtual void                accept(void) = 0;
        double                      acceptanceProbability(void) { return (double)numAccepted / numTries; }
        void                        clear(void);
        std::string                 getName(void) { return moveName; }
        int                         getNumTries(void) { return numTries; }
        virtual void                reject(void) = 0;
        virtual void                restore(void) = 0;
        virtual double              update(void) = 0;
        virtual void                tune(void) = 0;
        virtual std::vector<double> values(void) = 0;
    
    protected:
        Model*                      model;
        RandomVariable*             rv;
        std::string                 moveName;
        int                         numTries;
        int                         numAccepted;
};

#endif
