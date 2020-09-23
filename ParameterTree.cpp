#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include "Branch.h"
#include "BranchFactory.h"
#include "Msg.h"
#include "Node.h"
#include "ParameterTree.h"
#include "RandomVariable.h"



ParameterTree::ParameterTree(RandomVariable* rp, Model* mp, std::string nm, std::vector<std::string> tn, std::string fn, double lam) : Parameter(rp, mp, nm) {

    lambda = lam;
    buildTreeFromNewickString(tn, fn);
}

ParameterTree::ParameterTree(ParameterTree& t) : Parameter(t.rv, t.modelPtr, t.name) {

    clone(t);
}

ParameterTree::~ParameterTree(void) {

    for (Node* p : nodes)
        delete p;
    removeBranches();
}

Parameter& ParameterTree::operator=(Parameter& t) {

    if (this != &t)
        {
        ParameterTree* dc = dynamic_cast<ParameterTree*>(&t);
        clone(*dc);
        }
    return *this;
}

ParameterTree& ParameterTree::operator=(ParameterTree& t) {

    if (this != &t)
        {
        clone(t);
        }
    return *this;
}

Branch* ParameterTree::addBranch(Node* e1, Node* e2) {

    BranchFactory& bf = BranchFactory::branchFactoryInstance();
    Branch* b = bf.getBranch();
    b->setEnds(e1, e2);
    branches.insert(b);
    return b;
}

Node* ParameterTree::addNode(void) {

    Node* p = new Node;
    nodes.push_back(p);
    return p;
}

Node* ParameterTree::addNode(int idx) {

    Node* p = addNode();
    p->setIndex(idx);
    return p;
}

Node* ParameterTree::addNode(Node* nodeToCopy) {

    Node* p = addNode();
    p->setIndex( nodeToCopy->getIndex() ); // copy everything but pointers (neighbors & myBranch)
    p->setIsLeaf( nodeToCopy->getIsLeaf() );
    p->setName( nodeToCopy->getName() );
    p->setNeedsUpdate( nodeToCopy->needsUpdate() );
    p->setActiveCl( nodeToCopy->getActiveCl() );
    return p;
}

double ParameterTree::branchLengthSum(void) {

    double sum = 0.0;
    for (Branch* b : branches)
        sum += b->getLength();
    return sum;
}

void ParameterTree::buildTreeFromNewickString(std::vector<std::string> tn, std::string fn) {

    // copy the taxon names
    taxonNames = tn;
    
    std::string newickStr = readNewickTree(fn);

    std::vector<std::string> newickTokens;
    newickTokens = parseNewickString(newickStr);
    
    Node* p = NULL;
    int intIdx = (int)tn.size();
    bool readingBrlen = false;
    for (int i=0; i<newickTokens.size(); i++)
        {
        std::string token = newickTokens[i];
        if (token == "(")
            {
            if (p == NULL)
                {
                p = addNode(intIdx++);
                root = p;
                }
            else
                {
                Node* q = addNode(intIdx++);
                p->addNeighbor(q);
                q->addNeighbor(p);
                q->setAncestor(p);
                Branch* b = addBranch(p, q);
                b->setLength(1.0/lambda);
                p = q;
                p->setMyBranch(b);
                }
            p->setIsLeaf(false);
            readingBrlen = false;
            }
        else if (token == ")" || token == ",")
            {
            if (p->getAncestor() != NULL)
                p = p->getAncestor();
            else
                Msg::error("Tried to move down tree");
            readingBrlen = false;
            }
        else if (token == ":")
            {
            readingBrlen = true;
            }
        else if (token == ";")
            {
            
            }
        else
            {
            if (readingBrlen == false)
                {
                Node* q = addNode();
                p->addNeighbor(q);
                q->addNeighbor(p);
                q->setAncestor(p);
                Branch* b = addBranch(p, q);
                b->setLength(1.0/lambda);
                p = q;
                int idx = getIndexForTaxon(token);
                p->setIndex(idx);
                p->setName(token);
                p->setIsLeaf(true);
                p->setMyBranch(b);
                }
            else
                {
                Branch* b = findBranch(p, p->getAncestor());
                if (b != NULL)
                    {
                    double x = std::stod(token);
                    b->setLength(x);
                    }
                }
            readingBrlen = false;
            }
        }
    
    // randomly root the tree
    if (root->numNeighbors() > 2)
        {
        Branch* b = randomBranch();
        double v = b->getLength();
        Node* n1 = b->getEnd1();
        Node* n2 = b->getEnd2();
        n1->removeNeighbor(n2);
        n2->removeNeighbor(n1);
        removeBranch(b);
        
        Node* p = addNode(intIdx++);
        root = p;
        p->addNeighbor(n1);
        p->addNeighbor(n2);
        n1->addNeighbor(p);
        n2->addNeighbor(p);
        Branch* b1 = addBranch(n1, p);
        Branch* b2 = addBranch(n2, p);
        double v1 = v * rv->uniformRv();
        double v2 = v-v1;
        if (v1 < 0.0)
            v1 = -v1;
        if (v2 < 0.0)
            v2 = -v2;
        b1->setLength(v1);
        b2->setLength(v2);
        }
    
    // set the offset
    for (int i=0; i<nodes.size(); i++)
        nodes[i]->setOffset(i);

    // initialize the downpass sequence
    initializeDownPassSequence();
    
    // check that the tree is strictly bifurcating
    if (isBinary() == false)
        {
        std::cout << "   * Warning: tree is not binary" << std::endl;
        print();
        }
    else
        {
        std::cout << "   * Input tree is binary (good!)" << std::endl;
        }
}

void ParameterTree::clone(ParameterTree& t) {

    // copy prior info
    lambda = t.lambda;
    
    // copy the taxon names
    taxonNames = t.taxonNames;
    
    // set up the nodes
    if (nodes.size() != t.nodes.size())
        {
        for (size_t i=nodes.size(); i<t.nodes.size(); i++)
            nodes.push_back(new Node);
        }
    
    // set the root
    root = nodes[t.root->getOffset()];

    // copy node information
    for (int n=0; n<nodes.size(); n++)
        {
        Node* lftNode = nodes[n];
        Node* rhtNode = t.nodes[n];
        lftNode->setOffset(n);
        lftNode->setIndex( rhtNode->getIndex() );
        lftNode->setIsLeaf( rhtNode->getIsLeaf() );
        lftNode->setName( rhtNode->getName() );
        lftNode->setActiveCl( rhtNode->getActiveCl() );
        lftNode->setNeedsUpdate( rhtNode->needsUpdate() );
        
        lftNode->removeNeighbors();
        std::set<Node*>& rhtNeighbors = rhtNode->getNeighbors();
        for (Node* rn : rhtNeighbors)
            lftNode->addNeighbor( nodes[rn->getOffset()] );

        if (rhtNode->getAncestor() == NULL)
            lftNode->setAncestor(NULL);
        else
            lftNode->setAncestor( nodes[rhtNode->getAncestor()->getOffset()] );
        }
    
    // copy the downpass sequence
    if (downPassSequence.size() != t.downPassSequence.size())
        downPassSequence.resize(t.downPassSequence.size());
    for (int n=0; n<downPassSequence.size(); n++)
        downPassSequence[n] = nodes[ t.downPassSequence[n]->getOffset() ];

    // copy branch information
    removeBranches();
    for (int n=0; n<nodes.size(); n++)
        {
        if ( t.nodes[n]->getAncestor() != NULL )
            {
            Branch* rhtBranch = t.findBranch(t.nodes[n], t.nodes[n]->getAncestor());
            Branch* lftBranch = findBranch(nodes[n], nodes[n]->getAncestor());
            if (lftBranch == NULL)
                lftBranch = addBranch(nodes[n], nodes[n]->getAncestor());

            lftBranch->setLength( rhtBranch->getLength() );
            lftBranch->setActiveTi(  rhtBranch->getActiveTi() );
            lftBranch->setNeedsUpdate( rhtBranch->needsUpdate() );
            nodes[n]->setMyBranch(lftBranch);
            }
        }
}

Node* ParameterTree::copyNode(std::map<Node*, Node*>& visited, Node* orig) {

    if (orig == NULL)
        return NULL;
    
    std::map<Node*, Node*>::iterator it = visited.find(orig) ;
    if ( it != visited.end() )
        {
        return it->second;
        }

    Node* copy = addNode(orig);
    visited.insert( std::make_pair(orig, copy) );

    std::set<Node*> oNeighbors = orig->getNeighbors();
    for (std::set<Node*>::iterator it = oNeighbors.begin(); it != oNeighbors.end(); it++)
        {
        copy->addNeighbor( copyNode(visited, *it) );
        }
    copy->setAncestor( copyNode(visited, orig->getAncestor()) );
    return copy;
}

Branch* ParameterTree::findBranch(Node* e1, Node* e2) {

    if (e1 > e2)
        {
        Node* temp = e1;
        e1 = e2;
        e2 = temp;
        }

    for (Branch* b : branches)
        {
        if ( b->getEnd1() == e1 )
            {
            if ( b->getEnd2() == e2 )
                {
                return b;
                }
            }
        }
    
    return NULL;
}

Node* ParameterTree::findNodeIndexed(int idx) {

    for (Node* p : nodes)
        {
        if (p->getIndex() == idx)
            return p;
        }
    return NULL;
}

Node* ParameterTree::findNodeWithTaxonName(std::string n) {

    for (Node* nde : nodes)
        {
        if (nde->getName() == n)
            return nde;
        }
    return NULL;
}

void ParameterTree::flipAllActiveConditionalLikelihoods(void) {

    for  (Node* p : nodes)
        p->flipActiveCl();
}

void ParameterTree::flipAllActiveTransitionProbabilities(void) {

    for (Branch* b : branches)
        b->flipActiveTi();
}

std::set<Branch*> ParameterTree::getBranchesAround(Node* p) {

    std::set<Branch*> br;
    for (Branch* b : branches)
        {
        if (b->getEnd1() == p || b->getEnd2() == p)
            br.insert(b);
        }
    return br;
}

int ParameterTree::getIndexForTaxon(std::string n) {

    for (int i=0; i<taxonNames.size(); i++)
        {
        if (n == taxonNames[i])
            return i;
        }
    Msg::error("Cannot find taxon with name \"" + n + "\"");
    return 0;
}

std::string ParameterTree::getNewick(void) {

    std::stringstream ss;
    if (root->getIsLeaf() == true)
        {
        Node* oldRoot = root;
        std::vector<Node*> nbs = root->getNeighborsAsVector();
        if (nbs.size() > 1)
            Msg::error("Expecting only a single neighbor at the root of the tree");
        Node* newRoot = nbs[0];
        root->setAncestor(newRoot);
        oldRoot->setAncestor(newRoot);
        newRoot->setAncestor(NULL);
        root = newRoot;

        writeTree(root, ss);

        newRoot->setAncestor(oldRoot);
        oldRoot->setAncestor(NULL);
        root = oldRoot;
        }
    else
        {
        writeTree(root, ss);
        }
    std::string newick = ss.str();
    return newick;
}

Node* ParameterTree::getNodeForTaxon(std::string n) {

    for (std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); it++)
        {
        if ((*it)->getName() == n)
            return (*it);
        }
    return NULL;
}

size_t ParameterTree::getNumTaxa(void) {

    size_t x = 0;
    for (std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); it++)
        {
        if ( (*it)->getIsLeaf() == true )
            x++;
        }
    return x;
}

double ParameterTree::getTreeLength(void) {

    double sum = 0.0;
    for (Branch* b : branches)
        sum += b->getLength();
    return sum;
}

void ParameterTree::initializeDownPassSequence(void) {

    downPassSequence.clear();
    
    for (Node* n : nodes)
        n->setAncestor(NULL);
    passDown(root, root);
    if (root != NULL)
        root->setAncestor(NULL);
}

bool ParameterTree::isBinary(void) {

    for (Node* n : nodes)
        {
        if (n->getIsLeaf() == true)
            {
            
            }
        else if (n == root)
            {
            if (n->numNeighbors() != 2)
                return false;
            }
        else
            {
            if (n->numNeighbors() != 3)
                return false;
            }
        }
    return true;
}

void ParameterTree::listNodes(Node* p, Node* anc, size_t indent) {

    if (p != NULL)
        {
        std::set<Node*> neighbors = p->getNeighbors();
        
        for (int i=0; i<indent; i++)
            std::cout << " ";
        std::cout << p->getIndex() << " ( ";
        for (Node* n : neighbors)
            {
            if (n == p->getAncestor())
                std::cout << "a.";
            std::cout << n->getIndex() << " ";
            }
        std::cout << ") ";
        if (p->getMyBranch() != NULL)
            std::cout << std::fixed << std::setprecision(6) << p->getMyBranch()->getLength();
        else
            std::cout << "NULL";
        if (p->getIsLeaf() == true)
            std::cout << " (" << p->getName() << ")";
            
        //std::cout << " " << p->getIsLeaf();
        //std::cout << std::fixed << std::setprecision(5) << " " << p->getScratchVariable();
    
        std::cout << std::endl;

        for (std::set<Node*>::iterator it = neighbors.begin(); it != neighbors.end(); it++)
            {
            if ( (*it) != anc )
                listNodes( (*it), p, indent+3 );
            }
        }
}

double ParameterTree::lnPriorProb(void) {

    // find maximum depth to root
    double maxDepth = 0.0;
    for (int i=0; i<nodes.size(); i++)
        {
        Node* p = nodes[i];
        if (p->getIsLeaf() == true)
            {
            double depth = 0.0;
            Node* q = p;
            while (q != root)
                {
                depth += q->getMyBranch()->getLength();
                q = q->getAncestor();
                }
            if (depth > maxDepth)
                maxDepth = depth;
            }
        }
        
    return log(lambda) - lambda * maxDepth;
    
#   if 0
    double nBranches = (double)getNumBranches();
    double T = treeLength();
    double lnP = nBranches * log(lambda) - lambda * T;
    return lnP;
#   endif
}

std::vector<std::string> ParameterTree::parseNewickString(std::string ns) {

    std::vector<std::string> tks;
    for (int i=0; i<ns.size(); i++)
        {
        char c = ns[i];
        if (c == '(' || c == ')' || c == ',' || c == ':' || c == ';')
            {
            std::string tempStr;
            tempStr = c;
            tks.push_back(tempStr);
            }
        else
            {
            int j = i;
            std::string tempStr = "";
            while ( !(c == '(' || c == ')' || c == ',' || c == ':' || c == ';') )
                {
                tempStr += c;
                j++;
                c = ns[j];
                }
            i = j-1;
            tks.push_back(tempStr);
            }
        }
#   if 0
    std::cout << "The Newick string, broken into its parts:" << std::endl;
    for (int i=0; i<tks.size(); i++)
        std::cout << "   tks[" << i << "] = \"" << tks[i] << "\"" << std::endl;
#   endif
    
    return tks;
}

void ParameterTree::passDown(Node* p, Node* anc) {

    if (p != NULL)
        {
        std::set<Node*> neighbors = p->getNeighbors();
        for (std::set<Node*>::iterator it = neighbors.begin(); it != neighbors.end(); it++)
            {
            if ( (*it) != anc )
                passDown( (*it), p );
            }
        p->setMyBranch( findBranch(p, anc) );
        p->setAncestor(anc);
        downPassSequence.push_back(p);
        }
}

void ParameterTree::print(void) {

    listNodes(root, root, 3);
}

void ParameterTree::print(std::string header) {

    std::cout << header << std::endl;
    print();
}

void ParameterTree::printNodes(void) {

    for (std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); it++)
        (*it)->print();
    std::cout << "Down Pass Sequence: ";
    for (size_t i=0; i<downPassSequence.size(); i++)
        std::cout << downPassSequence[i]->getIndex() << " ";
    std::cout << std::endl;
}

Branch* ParameterTree::randomBranch(void) {

    size_t whichBranch = (int)(rv->uniformRv()*branches.size()), k = 0;
    for (Branch* b : branches)
        {
        if (k == whichBranch)
            return b;
        k++;
        }
    return NULL;
}

Branch* ParameterTree::randomInteriorBranch(void) {

    for (;;)
        {
        size_t whichBranch = (int)(rv->uniformRv()*branches.size()), k = 0;
        for (Branch* b : branches)
            {
            if (k == whichBranch && b->getDescendantNode()->getIsLeaf() == false)
                return b;
            k++;
            }
        }
    return NULL;
}

Node* ParameterTree::randomInternalNode(void) {

    for (;;)
        {
        Node* n = nodes[(int)(nodes.size()*rv->uniformRv())];
        if (n->getAncestor() != NULL && n->getIsLeaf() == false)
            return n;
        }
}

std::string ParameterTree::readNewickTree(std::string fn) {

    std::ifstream treeStream(fn.c_str());
    if (!treeStream)
        Msg::error("Cannot open tree file \"" + fn + "\"");

    std::string ns = "";
    
    if (getline(treeStream, ns).good() == false)
        {
        std::cout << ns << std::endl;
        Msg::error("Failed to read tree file \"" + fn + "\"");
        }
    
    treeStream.close();
    
    return ns;
}

void ParameterTree::removeBranch(Branch* b) {

    branches.erase(b);
    BranchFactory& bf = BranchFactory::branchFactoryInstance();
    bf.returnBranchToPool(b);
}

void ParameterTree::removeBranches(void) {

    BranchFactory& bf = BranchFactory::branchFactoryInstance();
    for (Branch* b : branches)
        bf.returnBranchToPool(b);
    branches.clear();
}

double ParameterTree::treeLength(void) {

    double sum = 0.0;
    for (Branch* b : branches)
        sum += b->getLength();
    return sum;
}

void ParameterTree::updateAllConditionalLikelihoods(bool tf)  {

    for (Node* p : nodes)
        p->setNeedsUpdate(tf);
}

void ParameterTree::updateAllTransitionProbabilities(bool tf) {

    for (Branch* b : branches)
        b->setNeedsUpdate(tf);
}

void ParameterTree::writeTree(Node* p, std::stringstream& ss) {

    if (p != NULL)
        {
        Branch* b = findBranch(p, p->getAncestor());
        if (p->getIsLeaf() == true)
            {
            ss << p->getIndex()+1;
            ss << ":" << b->getLength();
            }
        else
            {
            ss << "(";
            }
        std::vector<Node*> myDescendants = p->getDescendants();
        for (int i=0; i<(int)myDescendants.size(); i++)
            {
            writeTree(myDescendants[i], ss);
            if ( (i + 1) != (int)myDescendants.size() )
                ss << ",";
            }
        if (p->getIsLeaf() == false)
            {
            ss << ")";
            if (b != NULL)
                ss << ":" << b->getLength();
            }
        }
}

std::string ParameterTree::getParmString(int n) {

    return "";
}

std::string ParameterTree::getParmHeader(int n) {

    return "";
}


