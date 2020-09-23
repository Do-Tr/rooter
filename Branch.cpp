#include <iostream>
#include "Branch.h"
#include "Node.h"

#define MIN_BRANCH_LENGTH 0.00000001



Branch::Branch(void) {

    end1 = NULL;
    end2 = NULL;
    length = 0.0;
    activeTi = 0;
    updateTi = false;
}

Branch::Branch(Node* e1, Node* e2) {

    setEnds(e1, e2);
    length = 0.0;
    activeTi = 0;
    updateTi = false;
}

Branch::Branch(Node* e1, Node* e2, double len) {

    setEnds(e1, e2);
    length = len;
    activeTi = 0;
    updateTi = false;
}

Branch& Branch::operator=(const Branch& a) {

    if (this != &a)
        {
        end1 = a.end1;
        end2 = a.end2;
        length = a.length;
        activeTi = a.activeTi;
        updateTi = a.updateTi;
        }
    return *this;
}

bool Branch::operator==(const Branch& a) const {

    if ( end1 == a.end1 && end2 == a.end2 )
        return true;
    return false;
}

bool Branch::operator<(const Branch& a) const {

    if ( a.end1 < end1 )
        return true;
    else if ( a.end1 == end1 )
        {
        if ( a.end2 < end2 )
            return true;
        }
    return false;
}

void Branch::clean(void) {

    end1 = NULL;
    end2 = NULL;
    length = 0.0;
    activeTi = 0;
    updateTi = false;
}

void Branch::flipActiveTi(void) {

    if (activeTi == 0)
        activeTi = 1;
    else
        activeTi = 0;
}

Node* Branch::getAncestralNode(void) {

    if (end1->getAncestor() == end2)
        return end2;
    else if (end2->getAncestor() == end1)
        return end1;
    return NULL;
}

Node* Branch::getDescendantNode(void) {

    if (end1->isDescendant(end2) == true)
        return end2;
    else if (end2->isDescendant(end1) == true)
        return end1;
    return NULL;
}

bool Branch::isTip(void) {

    if (end1->getIsLeaf() == true || end2->getIsLeaf() == true)
        return true;
    return false;
}

void Branch::print(std::string header) {

    std::cout  << header << std::endl;
    print();
}

void Branch::print(void) {

    std::cout << "Branch: " << this << std::endl;
    std::cout << "   * length = " << length << std::endl;
    std::cout << "   * end1   = " << end1->getIndex() << std::endl;
    std::cout << "   * end2   = " << end2->getIndex() << std::endl;
}

void Branch::setEnds(Node* e1, Node* e2) {

    if (e1 < e2)
        {
        end1 = e1;
        end2 = e2;
        }
    else if (e2 < e1)
        {
        end1 = e2;
        end2 = e1;
        }
    else
        {
        std::cout << "Error: A branch is defined by two different nodes (" << e1 << ", " << e2 << ")" << std::endl;
        exit(1);
        }
}

void Branch::setLength(double x) {

    if (x < MIN_BRANCH_LENGTH)
        x = MIN_BRANCH_LENGTH;
    length = x;
}

bool CompBranch::operator()(Branch* b1, Branch* b2) const {

    if ( b1->getEnd1() < b2->getEnd1() )
        return true;
    else if ( b1->getEnd1() == b2->getEnd1() )
        {
        if ( b1->getEnd2() < b2->getEnd2() )
            return true;
        }
    return false;
}

