#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

class Tools {
public:
    /**
    * Constructor.
    */
    Tools();

    /**
    * Destructor.
    */
    virtual ~Tools();

    /**
    * A helper method to calculate RMSE.
    */
    Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, const std::vector<Eigen::VectorXd> &ground_truth);
    /**
     * helper method to convert coordinate of radar
     */
    Eigen::VectorXd Polar2Cartesian(const VectorXd& polar);
    Eigen::VectorXd Cartesian2Polar(const VectorXd& cartesian);
};

#endif /* TOOLS_H_ */
