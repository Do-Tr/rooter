#include <iostream>
#include "Node.h"



Node::Node(void) {

    index    = 0;
    isLeaf   = false;
    name     = "";
    myBranch = NULL;
    ancestor = NULL;
    activeCl = 0;
    updateCl = false;
    offset   = 0;
}

void Node::clean(void) {

    index    = 0;
    isLeaf   = false;
    name     = "";
    myBranch = NULL;
    ancestor = NULL;
    activeCl = 0;
    updateCl = false;
    offset   = 0;
    neighbors.clear();
}

void Node::flipActiveCl(void) {

    if (activeCl == 0)
        activeCl = 1;
    else
        activeCl = 0;
}

std::vector<Node*> Node::getDescendants(void) {

    std::vector<Node*> des;
    for (Node* n : neighbors)
        {
        if (n != ancestor)
            des.push_back(n);
        }
    return des;
}

std::vector<Node*> Node::getNeighborsAsVector(void) {

    std::vector<Node*> nb;
    for (Node* n : neighbors)
        nb.push_back( n );
    return nb;
}

bool Node::isDescendant(Node* p) {

    std::set<Node*>::iterator it = neighbors.find(p);
    if (it != neighbors.end() && p != ancestor)
        return true;
    return false;
}

void Node::print(void) {

    std::cout << "Node (" << this << ")" << std::endl;
    std::cout << "    Neighbors: ";
    for(Node* n : neighbors)
        {
        if (n == ancestor)
            std::cout << "a";
        std::cout << n->getIndex() << " ";
        }
    std::cout << std::endl;
    std::cout << "        index: " << index  << std::endl;
    std::cout << "       isLeaf: " << isLeaf << std::endl;
    std::cout << "         name: " << name   << std::endl;
}

