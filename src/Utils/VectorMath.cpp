#include "Utils/VectorMath.h"

#include <cmath>

namespace VectorMath
{

double dot 
(
    const std::vector<double>& a,
    const std::vector<double>& b
)
{
    double sum = 0.0;

    for(size_t i = 0; i < a.size(); ++i){
        sum += a[i] * b[i];
    }

    return sum;
}

void add
(
    const std::vector<double>& a,
    const std::vector<double>& b,
    std::vector<double>& r
)
{
    for(size_t i = 0; i < a.size(); ++i){
        r[i] = a[i] + b[i];
    }
}

void sub
(
    const std::vector<double>& a,
    const std::vector<double>& b,
    std::vector<double>& r
)
{
    for(size_t i = 0; i < a.size(); ++i){
        r[i] = a[i] - b[i];
    }
}

void daxpy
(
    double a,
    const std::vector<double>& x,
    std::vector<double>& y
)
{
    for(size_t i = 0; i < x.size(); ++i){
        y[i] += a * x[i];
    }
}

double L2Norm
(
    const std::vector<double>& a, 
    const std::vector<double>& b
)
{
    double norm = 0.0;
    for(int i = 0; i < a.size(); ++i)
    {
        norm += (a[i] - b[i]) * (a[i] - b[i]);
    }

    return std::sqrt(norm);
}

double L2Norm
(
    const std::vector<double>& r
)
{
    double norm = 0.0;
    for(int i = 0; i < r.size(); ++i)
    {
        norm += r[i]*r[i];
    }

    return std::sqrt(norm);
}
}