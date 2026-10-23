#pragma once

#include "Type.h"

#include <vector>
#include <map>
#include <string>

struct Node
{
    int nodeID;
    double x, y, z;

    //constructor
    Node(int nodeID_, double x_, double y_, double z_) : nodeID(nodeID_), x(x_), y(y_), z(z_) {}
};

struct Element
{
    int elementID;
    ElementType elementType;
    int physicalTag;

    std::vector<int> nodeIDs; //what the elements are made of

    //constructor
    Element(int eID_, ElementType eType_, int pTag_, std::vector<int> nID_) : 
                    elementID(eID_), elementType(eType_), physicalTag(pTag_), nodeIDs(nID_) {}
};

class Grid{
private:
    std::vector<Node> nodes_;
    std::vector<Element> elements_;
    std::map<int, std::string> boundaries_;
    std::map<int, std::string> domains_;

public:

    //construct from raw data read by mesh reader
    Grid
    (
        std::vector<Node> n,
        std::vector<Element> e,
        std::map<int, std::string> b,
        std::map<int, std::string> d
    )
    :
        nodes_(std::move(n)),
        elements_(std::move(e)),
        boundaries_(std::move(b)),
        domains_(std::move(d))
    {}

    //getters for other mesh classes
    const std::vector<Node>& nodes() const {return nodes_; }
    const std::vector<Element>& elements() const {return elements_; }
    const std::map<int, std::string>& boundaries() const {return boundaries_; }
    const std::map<int, std::string>& domians() const {return domains_; }
};