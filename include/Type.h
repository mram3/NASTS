#pragma once

#include <string>

enum class ElementType
{
    NONE = -1,
    Line = 1,
    Triangle = 2,
    Quad = 3
};

enum class BoundaryConditionType
{
    Dirichlet,
    Neumann,
    None
};

template<typename T>
struct BoundaryCondition
{
    BoundaryConditionType bc = BoundaryConditionType::None;
    T value =  {};

};

struct TwoVector
{
    double x = 0.0, y = 0.0;

};

struct SolverPerf
{
    double initRes;
    double finalRes;
    int iter;
};