#ifndef BipartitionSamples_H
#define BipartitionSamples_H

#include <vector>


class BipartitionSamples {

    public:
                                    BipartitionSamples(void);
        void                        addBipartition(std::string bp);
        void                        print(void);

    protected:
        int                         nunSamples;
        std::vector<std::string>    bitsets;
        std::vector<int>            numVisits;
};

#endif
