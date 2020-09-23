#ifndef RootBipartition_H
#define RootBipartition_H

#include <string>
#include <vector>
class ParameterTree;
class RbBitSet;



class RootBipartition {

    public:
                                RootBipartition(void);
                               ~RootBipartition(void);
        std::string             getRootBipartition(ParameterTree* t);
        void                    print(void);

    private:
        void                    deleteBipartitions(void);
        std::vector<RbBitSet*>  partitions;
};

#endif
