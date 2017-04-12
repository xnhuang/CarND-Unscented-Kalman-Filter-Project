#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    VectorXd rmse(4);
    rmse << 0,0,0,0;

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

VectorXd Tools::Polar2Cartesian(const VectorXd& polar) {
    double rho = polar(0);
    double phi = polar(1);
    double rho_dot = polar(2);
    VectorXd cartesian = VectorXd(4);
    cartesian << rho*cos(phi), rho*sin(phi), 0, 0;
    return cartesian;
}

VectorXd Tools::Cartesian2Polar(const VectorXd& cartesian) {
    double rho2 = pow(cartesian(0),2)+pow(cartesian(1),2);
    double phi = atan2(cartesian(1),cartesian(0));
    double rho_dot = cos(phi)*cartesian(2)+sin(phi)*cartesian(3);
    VectorXd polar = VectorXd(3);
    polar << sqrt(rho2), phi, rho_dot;
    return polar;
}