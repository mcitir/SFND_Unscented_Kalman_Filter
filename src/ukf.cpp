#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

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
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */

  is_initialized_ = false;

  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_.fill(0.5 / (lambda_ + n_aug_));
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  R_radar_ = MatrixXd(3, 3);
  R_radar_.fill(0.0);
  R_radar_(0, 0) = std_radr_ * std_radr_;
  R_radar_(1, 1) = std_radphi_ * std_radphi_;
  R_radar_(2, 2) = std_radrd_ * std_radrd_;

  R_lidar_ = MatrixXd(2, 2);
  R_lidar_.fill(0.0);
  R_lidar_(0, 0) = std_laspx_ * std_laspx_;
  R_lidar_(1, 1) = std_laspy_ * std_laspy_;

  debug_ = false;

  if(debug_) {
    cout << "UKF::UKF()" << std::endl;
    cout << "n_x_ = " << n_x_ << std::endl;
    cout << "n_aug_ = " << n_aug_ << std::endl;
    cout << "lambda_ = " << lambda_ << std::endl;
    cout << "weights_ = " << weights_ << std::endl;
    cout << "R_radar_ = " << R_radar_ << std::endl;
    cout << "R_lidar_ = " << R_lidar_ << std::endl;
  }

}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */

  if(!is_initialized_){
    if(debug_) cout << "UKF::ProcessMeasurement() - Initializing" << std::endl;

    if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
      if(debug_) cout << "UKF::ProcessMeasurement() - Initializing with RADAR" << std::endl;
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double rho_dot = meas_package.raw_measurements_[2];

      double px = rho * cos(phi);
      double py = rho * sin(phi);
      double vx = rho_dot * cos(phi);
      double vy = rho_dot * sin(phi);
      double v = sqrt(vx * vx + vy * vy);

      x_ << px, py, v, 0, 0;
    } else if(meas_package.sensor_type_ == MeasurementPackage::LASER){
      if(debug_) cout << "UKF::ProcessMeasurement() - Initializing with LASER" << std::endl;
      double px = meas_package.raw_measurements_[0];
      double py = meas_package.raw_measurements_[1];

      x_ << px, py, 0, 0, 0;
    }

    P_ << 1, 0, 0, 0, 0,
          0, 1, 0 ,0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0 ,1, 0,
          0, 0, 0, 0, 1;

    cout << "UKF::ProcessMeasurement() - Initialized x_ = " << x_ << std::endl;
    cout << "UKF::ProcessMeasurement() - Initialized P_ = " << P_ << std::endl;
    

    time_us_ = meas_package.timestamp_;

    is_initialized_ = true;
    return;
  }

  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;
  cout << "UKF::ProcessMeasurement() - delta_t = " << delta_t << std::endl;

  Prediction(delta_t);
  cout << "UKF::ProcessMeasurement() - Prediction done" << std::endl;

  if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
    cout << "UKF::ProcessMeasurement() - Update Radar" << std::endl;
    UpdateRadar(meas_package);
  } else if(meas_package.sensor_type_ == MeasurementPackage::LASER){
    cout << "UKF::ProcessMeasurement() - Update Lidar" << std::endl;
    UpdateLidar(meas_package);
  }



}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;

  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.fill(0.0);

  MatrixXd L = P_aug.llt().matrixL();

  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; i++){
    Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  for (int i = 0; i < 2 * n_aug_ + 1; i++){
    double px = Xsig_aug(0, i);
    double py = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    double px_p, py_p, v_p, yaw_p, yawd_p;

    if(fabs(yawd) > 0.001){
      px_p = px + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
      py_p = py + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    } else {
      px_p = px + v * delta_t * cos(yaw);
      py_p = py + v * delta_t * sin(yaw);
    }

    v_p = v;
    yaw_p = yaw + yawd * delta_t;
    yawd_p = yawd;

    // add noise
    px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;

    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p = yawd_p + nu_yawdd * delta_t;

    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;

  }

  // predict mean and covariance
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++){
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++){
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while(x_diff(3) > M_PI) x_diff(3) -= 2.0 * M_PI;
    while(x_diff(3) < -M_PI) x_diff(3) += 2.0 * M_PI;
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }

  cout << "UKF::Prediction() - x_ = " << x_ << std::endl;
  cout << "UKF::Prediction() - P_ = " << P_ << std::endl;

}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */

  int n_z = 2;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  Zsig.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++){
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);

    Zsig(0, i) = px;
    Zsig(1, i) = py;
  }

  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++){
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++){
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  S = S + R_lidar_;

  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++){
    VectorXd z_diff = Zsig.col(i) - z_pred;
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  MatrixXd K = Tc * S.inverse();

  int n_z_ = 2;
  VectorXd z = VectorXd(n_z_);
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];
  VectorXd z_diff = z - z_pred;

  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  cout << "UKF::UpdateLidar() - x_ = " << x_ << std::endl;
  cout << "UKF::UpdateLidar() - P_ = " << P_ << std::endl;

  nis_lidar_ = z_diff.transpose() * S.inverse() * z_diff;

}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

  int n_z = 3;
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], meas_package.raw_measurements_[2];

  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  Zsig.fill(0.0);

  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);

  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++){
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    double vx = cos(yaw) * v;
    double vy = sin(yaw) * v;

    Zsig(0, i) = sqrt(px * px + py * py);
    Zsig(1, i) = atan2(py, px);
    Zsig(2, i) = (px * vx + py * vy) / sqrt(px * px + py * py);

    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  for (int i = 0; i < 2 * n_aug_ + 1; i++){
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    while(z_diff(1) > M_PI) z_diff(1) -= 2.0 * M_PI;
    while(z_diff(1) < -M_PI) z_diff(1) += 2.0 * M_PI;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  S = S + R_radar_;

  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++){
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    while(z_diff(1) > M_PI) z_diff(1) -= 2.0 * M_PI;
    while(z_diff(1) < -M_PI) z_diff(1) += 2.0 * M_PI;

    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while(x_diff(3) > M_PI) x_diff(3) -= 2.0 * M_PI;
    while(x_diff(3) < -M_PI) x_diff(3) += 2.0 * M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  MatrixXd K = Tc * S.inverse();

  VectorXd z_diff = z - z_pred;

  // angle normalization
  while(z_diff(1) > M_PI) z_diff(1) -= 2.0 * M_PI;
  while(z_diff(1) < -M_PI) z_diff(1) += 2.0 * M_PI;

  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  cout << "UKF::UpdateRadar() - x_ = " << x_ << std::endl;
  cout << "UKF::UpdateRadar() - P_ = " << P_ << std::endl;

  nis_radar_ = z_diff.transpose() * S.inverse() * z_diff;
}