#include <iostream>
#include "Msg.h"
#include "Node.h"
#include "ParameterTree.h"
#include "RbBitSet.h"
#include "RootBipartition.h"



RootBipartition::RootBipartition(void) {

}

RootBipartition::~RootBipartition(void) {

    deleteBipartitions();
}

void RootBipartition::deleteBipartitions(void) {

    for (int i=0; i<partitions.size(); i++)
        delete partitions[i];
    partitions.clear();
}

std::string RootBipartition::getRootBipartition(ParameterTree* t) {

    int n = (int)t->getNumTaxa();

    // check that we have enough partitions alreay allocated
    if (partitions.size() != t->getNumNodes())
        {
        deleteBipartitions();
        for (int i=0; i<t->getNumNodes(); i++)
            partitions.push_back( new RbBitSet(n, false) );
        }
    
    // check that the bipartitions are of the right size
    if (partitions[0]->size() != n)
        {
        deleteBipartitions();
        for (int i=0; i<t->getNumNodes(); i++)
            partitions.push_back( new RbBitSet(n, false) );
        }
    
    // fill in the bipartitions
    std::vector<Node*> dpSeq = t->getDownPassSequence();
    for (int i=0; i<dpSeq.size(); i++)
        {
        Node* p = dpSeq[i];
        if (p->getIsLeaf() == true)
            {
            RbBitSet* bs = partitions[p->getIndex()];
            bs->unset();
            bs->set(p->getIndex());
            }
        else
            {
            RbBitSet* bs = partitions[p->getIndex()];
            bs->unset();
            std::vector<Node*> desc = p->getDescendants();
            for (int j=0; j<desc.size(); j++)
                {
                (*bs) |= (*partitions[desc[j]->getIndex()]);
                }
            }
        }
    
    // report the root parition
    std::vector<Node*> rootDescendants = t->getRoot()->getDescendants();
    if (rootDescendants.size() != 2)
        Msg::error("Root node should have only two descendants, but instead " + std::to_string(rootDescendants.size()) + " were found");
    
    std::string s = "";
    for (int i=0; i<2; i++)
        {
        Node* p = rootDescendants[i];
        if (partitions[p->getIndex()]->isSet(0) == false)
            {
            s = partitions[p->getIndex()]->bitString();
            return s;
            }
        }
    return "blank";
}

void RootBipartition::print(void) {

    for (int i=0; i<partitions.size(); i++)
        {
        std::cout << i << " -- " << partitions[i]->bitString() << std::endl;
        }
}
