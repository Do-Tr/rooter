#ifndef Node_H
#define Node_H

#include <set>
#include <string>
#include <vector>
class Branch;



class Node {

    public:
                            Node(void);
        void                addNeighbor(Node* p) { neighbors.insert(p); }
        void                clean(void);
        int                 degree(void) { return (int)neighbors.size(); }
        void                flipActiveCl(void);
        int                 getActiveCl(void) { return activeCl; }
        Node*               getAncestor(void) { return ancestor; }
        std::vector<Node*>  getDescendants(void);
        Branch*             getMyBranch(void) { return myBranch; }
        std::set<Node*>&    getNeighbors(void) { return neighbors; }
        std::set<Node*>     getNeighborsCopy(void) { return neighbors; }
        std::vector<Node*>  getNeighborsAsVector(void);
        int                 getIndex(void) { return index; }
        bool                getIsLeaf(void) { return isLeaf; }
        std::string         getName(void) { return name; }
        int                 getOffset(void) { return offset; }
        double              getScratchVariable(void) { return scratchVariable; }
        bool                isDescendant(Node* p);
        bool                needsUpdate(void) { return updateCl; }
        size_t              numNeighbors(void) { return neighbors.size(); }
        void                print(void);
        void                removeNeighbor(Node* p) { neighbors.erase(p); }
        void                removeNeighbors(void) { neighbors.clear(); }
        void                setActiveCl(int x) { activeCl = x; }
        void                setAncestor(Node* p) { ancestor = p; }
        void                setIndex(int x) { index = x; }
        void                setIsLeaf(bool tf) { isLeaf = tf; }
        void                setOffset(int x) { offset = x; }
        void                setMyBranch(Branch* b) { myBranch = b; }
        void                setName(std::string s) { name = s; }
        void                setNeedsUpdate(bool tf) { updateCl = tf; }
        void                setScratchVariable(double x) { scratchVariable = x; }

    protected:
        std::set<Node*>     neighbors;
        Node*               ancestor;
        int                 index;
        int                 offset;
        bool                isLeaf;
        std::string         name;
        Branch*             myBranch;
        int                 activeCl;
        bool                updateCl;
        double              scratchVariable;
};


#endif
