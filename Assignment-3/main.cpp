#include <iostream>
#include <Eigen/Dense>
#include "math/math.h"

//Assignment 3
//Task 1: Functions and algorithms
//T.1 a) Create a function that converts between vector types.
Eigen::VectorXd std_vector_to_eigen(const std::vector<double> &v) {
    Eigen::VectorXd r(v.size());
    for (int i = 0; i < v.size(); i++)
        r[i]=v[i];
    return r;
}

//T.1 b) Compare a number to the average value of the last items in a collection.
bool is_average_below_eps(const std::vector<double> &values, double eps, uint8_t n_values) {
    if (values.size() < n_values)
        return false;
    else {
        double average = 0.0;
        for (int i = (values.size()-n_values); i < values.size(); i++) {average=average+values[i];}
        average = average/n_values;
        if (average > eps)
            return false;
        else
            return true;
    }
}


int main() {
    double eps = 0.1;
    return 0;
}