#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include "measurement_package.h"
#include <iostream>
#include <math.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    // initial state vector
    x_ = VectorXd(5);

    // initial covariance matrix
    P_ = MatrixXd(5, 5);

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 0.8;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 0.6;

    // Laser measurement noise standard deviation position1 in m
    std_laspx_ = 0.15;

    // Laser measurement noise standard deviation position2 in m
    std_laspy_ = 0.15;

    // Radar measurement noise standard deviation radius in m
    std_radr_ = 0.3;

    // Radar measurement noise standard deviation angle in rad
    std_radphi_ = 0.03;

    // Radar measurement noise standard deviation radius change in m/s
    std_radrd_ = 0.3;

    // State dimension
    n_x_ = 5;

    // Augmented state dimension
    n_aug_ = n_x_+2;

    // predicted sigma points matrix
    Xsig_aug_ = MatrixXd(n_aug_, 2*n_aug_ + 1);

    // predicted sigma points matrix
    Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_ + 1);

    // Sigma point spreading parameter
    lambda_ = 3 - n_aug_;

    // Weights of sigma points
    weights_ = VectorXd(2*n_aug_ + 1);
    weights_.fill(0.5/(lambda_+n_aug_));
    weights_(0)=lambda_/(lambda_+n_aug_);

    // the current NIS for radar
    NIS_radar_ = 0.0;

    // the current NIS for laser
    NIS_laser_ = 0.0;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    if(meas_package.raw_measurements_.isZero(0)) return;
    if (!is_initialized_) {
        /**
          * Initialize the state ekf_.x_ with the first measurement.
          * Create the covariance matrix.
          * Remember: you'll need to convert radar from polar to cartesian coordinates.
        */
        // first measurement
        cout << "weight= " << weights_ << endl;
        cout << "UKF: " << endl;
        VectorXd x_initial_ = VectorXd(n_x_);
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            x_initial_[0] = meas_package.raw_measurements_(0)*cos(meas_package.raw_measurements_(1));
            x_initial_[1] = meas_package.raw_measurements_(0)*sin(meas_package.raw_measurements_(1));
            x_initial_[2]=0;
            x_initial_[3]=0;
            x_initial_[4]=0;
        }
        else{
            x_initial_<< meas_package.raw_measurements_(0),meas_package.raw_measurements_(1),0,0,0;;
        }
        //state covariance matrix P
        P_ << 1, 0, 0, 0, 0,
                0, 1, 0, 0, 0,
                0, 0, 1, 0, 0,
                0, 0, 0, 1, 0,
                0, 0, 0, 0, 1;
        x_ = x_initial_;
        cout << "Kalman Filter Initialization " << endl;

        time_us_ = meas_package.timestamp_;

        // done initializing, no need to predict or update
        is_initialized_ = true;
        cout << "x_initial_ = " << x_ << endl;
        cout << "P_initial_ = " << P_ << endl;
        return;
    }

    /*****************************************************************************
     *  Prediction
     ****************************************************************************/
    //compute the time elapsed between the current and previous measurements
    if ((meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
        || (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)) {
        double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;    //dt - expressed in seconds
        cout << "dt = " << dt << endl;
        time_us_ = meas_package.timestamp_;

        Prediction(dt);
    }
    /*****************************************************************************
    *  Update
    ****************************************************************************/

    /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
    */

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
        // Radar updates
        UpdateRadar(meas_package);
        cout << "x_update_radar_ = " << x_ << endl;
        cout << "P_update_radar__ = " << P_ << endl;
    } else {
        // Laser updates
        if (use_laser_){
            UpdateLidar(meas_package);
            cout << "x_update_laser_ = " << x_ << endl;
            cout << "P_update_laser__ = " << P_ << endl;
        }
    }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
    //create augmented mean vector
    VectorXd x_aug = VectorXd(n_aug_);

    //create augmented state covariance
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

    //create augmented mean state
    x_aug.fill(0.0);
    x_aug.head(n_x_) = x_;

    //create augmented covariance matrix
    P_aug.fill(0.0);
    P_aug.topLeftCorner(n_x_, n_x_)=P_;
    P_aug(n_x_,n_x_)=std_a_*std_a_;
    P_aug(n_x_+1,n_x_+1)=std_yawdd_*std_yawdd_;

    MatrixXd A = P_aug.llt().matrixL();
    Xsig_aug_.col(0) = x_aug;
    for (int i = 0; i < n_aug_; i++)
    {
        Xsig_aug_.col(i+1)        = x_aug + sqrt(lambda_+n_aug_) * A.col(i);
        Xsig_aug_.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * A.col(i);
    }
    P_.fill(0.0);
    //predict state mean
    x_.fill(0.0);
    for(int aug_id=0;aug_id<(2*n_aug_+1);aug_id++){

        double p_x = Xsig_aug_(0,aug_id);
        double p_y = Xsig_aug_(1,aug_id);
        double v = Xsig_aug_(2,aug_id);
        double phi = Xsig_aug_(3,aug_id);
        double phi_dot = Xsig_aug_(4,aug_id);
        double nu_spd = Xsig_aug_(5,aug_id);
        double nu_angle = Xsig_aug_(6,aug_id);

        double p_x_update,p_y_update;
        if(fabs(phi_dot)> 0.001){
            p_x_update = p_x + v/phi_dot * (sin(phi + phi_dot*delta_t) - sin(phi));
            p_y_update = p_y + v/phi_dot * (cos(phi) - cos(phi + phi_dot*delta_t));
        }else{
            p_x_update = p_x + v*cos(phi)*delta_t;
            p_y_update = p_y + v*sin(phi)*delta_t;
        }
        double v_update = v;
        double phi_update = phi + phi_dot*delta_t;
        double phi_dot_update = phi_dot;

        p_x_update += 0.5*cos(phi)*delta_t*delta_t*nu_spd;
        p_y_update += 0.5*sin(phi)*delta_t*delta_t*nu_spd;
        v_update += delta_t*nu_spd;
        phi_update += 0.5*delta_t*delta_t*nu_angle;
        phi_dot_update += delta_t*nu_angle;

        Xsig_pred_(0,aug_id) = p_x_update;
        Xsig_pred_(1,aug_id) = p_y_update;
        Xsig_pred_(2,aug_id) = v_update;
        Xsig_pred_(3,aug_id) = phi_update;
        Xsig_pred_(4,aug_id) = phi_dot_update;

        x_ += weights_(aug_id) * Xsig_pred_.col(aug_id);
    }
    //predict state covariance matrix
    x_(3) = tools_.normalize_angle(x_(3));
    for(int col_id=0;col_id<(2 * n_aug_ + 1);col_id++){
        VectorXd state_center = Xsig_pred_.col(col_id) - x_;
        state_center(3) = tools_.normalize_angle(state_center(3));
        P_ += weights_(col_id)*state_center*(state_center.transpose());
    }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

    int n_z = 2;
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    VectorXd measure_prediction = VectorXd(n_z);
    MatrixXd S = MatrixXd(n_z,n_z);
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

    Tc.fill(0.0);
    S.fill(0.0);
    S.diagonal() << pow(std_laspx_,2),pow(std_laspy_,2);
    measure_prediction.fill(0.0);

    for(int col_id=0;col_id<(2 * n_aug_ + 1);col_id++){
        double px = Xsig_pred_(0,col_id);
        double py = Xsig_pred_(1,col_id);

        VectorXd lidar_point(n_z);
        lidar_point << px,py;

        measure_prediction += weights_(col_id)*lidar_point;
        Zsig.col(col_id) = lidar_point;
    }

    for(int col_id=0;col_id<(2 * n_aug_ + 1);col_id++){
        VectorXd z_center = Zsig.col(col_id) - measure_prediction;
        VectorXd x_center = Xsig_pred_.col(col_id) - x_;
        x_center(3) = tools_.normalize_angle(x_center(3));
        S += weights_(col_id)*z_center*(z_center.transpose());
        Tc += weights_(col_id)*x_center*(z_center.transpose());
    }
//    cout << "determinate of S " << S.determinant() << endl;
//    cout << "measure_prediction " << measure_prediction << endl;
//    cout << "measure " << meas_package.raw_measurements_ << endl;
    //inverse S
    MatrixXd S_inverse = S.llt().solve(MatrixXd::Identity(n_z,n_z));
    MatrixXd K = Tc*S_inverse;
    //update state mean and covariance matrix
    x_ += K*(meas_package.raw_measurements_-measure_prediction);
    x_(3) = tools_.normalize_angle(x_(3));
    P_ -= K*S*K.transpose();

    NIS_laser_ = (Zsig.col(0)-measure_prediction).transpose()*S_inverse*(Zsig.col(0)-measure_prediction);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

    int n_z = 3;
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    VectorXd measure_prediction = VectorXd(n_z);
    MatrixXd S = MatrixXd(n_z,n_z);
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

    Tc.fill(0.0);
    S.fill(0.0);
    S.diagonal()<<pow(std_radr_,2),pow(std_radphi_,2),pow(std_radrd_,2);

    measure_prediction.fill(0.0);

    for(int col_id=0;col_id<(2 * n_aug_ + 1);col_id++){
        double px = Xsig_pred_(0,col_id);
        double py = Xsig_pred_(1,col_id);
        double v = Xsig_pred_(2,col_id);
        double yaw = Xsig_pred_(3,col_id);

        double rho = sqrt(pow(px,2)+pow(py,2));
        double radar_angle = atan2(py,px);
        double rho_dot = v*cos(yaw)*cos(radar_angle) + v*sin(yaw)*sin(radar_angle);

        VectorXd radar_point(n_z);
        radar_point << rho,radar_angle,rho_dot;

        measure_prediction += weights_(col_id)*radar_point;
        Zsig.col(col_id) = radar_point;
    }
    measure_prediction(1) = tools_.normalize_angle(measure_prediction(1));
    for(int col_id=0;col_id<(2 * n_aug_ + 1);col_id++){
        VectorXd z_center = Zsig.col(col_id) - measure_prediction;
        VectorXd x_center = Xsig_pred_.col(col_id) - x_;
        z_center(1) = tools_.normalize_angle(z_center(1));
        x_center(3) = tools_.normalize_angle(x_center(3));
        S += weights_(col_id)*z_center*(z_center.transpose());
        Tc += weights_(col_id)*x_center*(z_center.transpose());
    }
//    cout << "determinate of S " << S.determinant() << endl;
//    cout << "measure_prediction " << measure_prediction << endl;
//    cout << "measure " << meas_package.raw_measurements_ << endl;
    //calculate Kalman gain K;
    MatrixXd S_inverse = S.llt().solve(MatrixXd::Identity(n_z,n_z));
    MatrixXd K = Tc*S_inverse;
    //update state mean and covariance matrix
    x_ += K*(meas_package.raw_measurements_-measure_prediction);
    x_(3) = tools_.normalize_angle(x_(3));
    P_ -= K*S*K.transpose();
    NIS_radar_ = (Zsig.col(0)-measure_prediction).transpose()*S_inverse*(Zsig.col(0)-measure_prediction);
}
