#include <iostream>
#include "tools.h"
#include <math.h>
using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    VectorXd rmse(estimations[0].size());
    rmse.fill(0.0);

    if (estimations.size()==0){
        std::cout << "estimation size 0" << std::endl;
        return rmse;
    }
    if (estimations.size()!=ground_truth.size()){
        std::cout << "estimation size small" << std::endl;
        return rmse;
    }
    //accumulate squared residuals
    for(int i=0; i < estimations.size(); ++i){
        VectorXd residual = estimations[i] - ground_truth[i];
        residual = residual.array()*residual.array();
        rmse += residual;
    }
    //calculate the mean
    rmse /= estimations.size();
    //calculate the squared root
    rmse=rmse.array().sqrt();

    //return the result
    return rmse;
}

double Tools::normalize_angle(double angle) {
    if (angle > M_PI){
        angle = fmod((angle - M_PI),(2*M_PI)) - M_PI;
    }
    if (angle < -M_PI){
        angle = fmod((angle + M_PI),(2*M_PI)) + M_PI;
    }
    return angle;
}
