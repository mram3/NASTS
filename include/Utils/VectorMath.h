#pragma once

#include <vector>

namespace VectorMath{

    double dot
    (
        const std::vector<double>& a, 
        const std::vector<double>& b
    );

    void add
    (
        const std::vector<double>& a, 
        const std::vector<double>& b,
        std::vector<double>& r
    );

    void sub
    (
        const std::vector<double>& a, 
        const std::vector<double>& b,
        std::vector<double>& r
    );

    void daxpy
    (
        double a,
        const std::vector<double>& x,
        std::vector<double>& y
    );

    double L2Norm
    (
        const std::vector<double>& a, 
        const std::vector<double>& b
    );

}