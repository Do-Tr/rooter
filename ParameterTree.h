#ifndef ParameterTree_H
#define ParameterTree_H

#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include "Branch.h"
#include "Parameter.h"
class Node;


class ParameterTree : public Parameter {

    public:
                                        ParameterTree(RandomVariable* rp, Model* mp, std::string nm, std::vector<std::string> tn, std::string fn, double lam);
                                        ParameterTree(ParameterTree& t);
                                       ~ParameterTree(void);
        Parameter&                      operator=(Parameter& t);
        ParameterTree&                  operator=(ParameterTree& t);
        Branch*                         addBranch(Node* e1, Node* e2);
        Branch*                         findBranch(Node* e1, Node* e2);
        Node*                           findNodeIndexed(int idx);
        Node*                           findNodeWithTaxonName(std::string n);
        void                            flipAllActiveConditionalLikelihoods(void);
        void                            flipAllActiveTransitionProbabilities(void);
        std::set<Branch*>               getBranchesAround(Node* p);
        std::vector<Node*>&             getDownPassSequence(void) { return downPassSequence; }
        int                             getNumBranches(void) { return (int)branches.size(); }
        size_t                          getNumNodes(void) { return nodes.size(); }
        size_t                          getNumTaxa(void);
        std::string                     getNewick(void);
        std::string                     getParmHeader(int n);
        std::string                     getParmString(int n);
        Node*                           getRoot(void) { return root; }
        std::vector<std::string>&       getTaxonNames(void) { return taxonNames; }
        double                          getTreeLength(void);
        void                            initializeDownPassSequence(void);
        bool                            isBinary(void);
        double                          lnPriorProb(void);
        void                            print(void);
        void                            print(std::string header);
        void                            printNodes(void);
        Branch*                         randomBranch(void);
        Branch*                         randomInteriorBranch(void);
        Node*                           randomInternalNode(void);
        void                            removeBranch(Branch* b);
        double                          treeLength(void);
        void                            updateAllConditionalLikelihoods(bool tf);
        void                            updateAllTransitionProbabilities(bool tf);

    protected:
        Node*                           addNode(void);
        Node*                           addNode(int idx);
        Node*                           addNode(Node* nodeToCopy);
        double                          branchLengthSum(void);
        void                            buildTreeFromNewickString(std::vector<std::string> tn, std::string fn);
        void                            clone(ParameterTree& b);
        Node*                           copyNode(std::map<Node*, Node*>& visited, Node* orig);
        int                             getIndexForTaxon(std::string n);
        Node*                           getNodeForTaxon(std::string n);
        void                            listNodes(Node* p, Node* anc, size_t indent);
        std::vector<std::string>        parseNewickString(std::string ns);
        void                            passDown(Node* p, Node* anc);
        std::string                     readNewickTree(std::string fn);
        void                            removeBranches(void);
        void                            writeTree(Node* p, std::stringstream& ss);

        Node*                           root;
        std::vector<Node*>              nodes;
        std::set<Branch*,CompBranch>    branches;
        std::vector<std::string>        taxonNames;
        std::vector<Node*>              downPassSequence;
        double                          lambda;
};

#endif
