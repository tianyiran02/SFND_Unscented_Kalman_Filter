#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

#define N_RADAR_MEAS    3
#define N_LIDAR_MEAS    2
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
  std_a_ = 5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.5;
  
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
  // Init state parameters
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;

  // Init weights
  weights = VectorXd(2 * n_aug_ + 1);
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  double weight = 0.5 / (lambda_ + n_aug_);
  weights(0) = weight_0;
  for (int i=1; i< 2 * n_aug_ + 1; ++i) {  
    weights(i) = weight;
  }

  // Init matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.0);

  P_aug_ = MatrixXd(n_aug_, n_aug_);
  P_aug_.fill(0.0);
  x_.fill(0.0);

  overcnt_ = 0;
  allcnt_ = 0;

  std::cout << "==========CP1" << std::endl;

  is_initialized_ = false;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  std::cout << "============= ROUND ===================" << std::endl;
  std::cout << "previous x_:\n" << x_ << std::endl;
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  std::cout << "Measurements: \n" << meas_package.raw_measurements_ << std::endl;
  if (meas_package.sensor_type_ == meas_package.LASER)
  {

    if (is_initialized_ == false)
    {
      std::cout << "==========CP2.1.0" << std::endl;

      x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 0, 0, 0;

      P_ << std_laspx_ * std_laspx_, 0, 0, 0, 0,
            0, std_laspy_ * std_laspy_, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;

      // update time
      time_us_ = meas_package.timestamp_;

      is_initialized_ = true;
    }
    else
    {
      // here we need to update hte timestamp. The timestamp origin unit is us. We need to map to s.
      double dt = double(meas_package.timestamp_ - time_us_) / 1e6;
      time_us_ = meas_package.timestamp_;
      Prediction(dt);

      std::cout << "==========CP2.2" << std::endl;
      // do radar process
      // create matrix for sigma points in measurement space
      MatrixXd Zsig = MatrixXd(N_LIDAR_MEAS, 2 * n_aug_ + 1);

      // mean predicted measurement
      VectorXd z_pred = VectorXd(N_LIDAR_MEAS);
      
      // measurement covariance matrix S
      MatrixXd S = MatrixXd(N_LIDAR_MEAS, N_LIDAR_MEAS);

      // transform sigma points into measurement space
      for (int i = 0; i < 2 * n_aug_ + 1; i ++)
      {
          double px = Xsig_pred_(0, i);
          double py = Xsig_pred_(1, i);

          Zsig(0, i) = px;
          Zsig(1, i) = py;
      }
      std::cout << "==========CP3" << std::endl;
      std::cout << "Zsig is:\n" << Zsig << std::endl;
      // calculate mean predicted measurement
      z_pred = Zsig * weights;

      // calculate innovation covariance matrix S
      MatrixXd R = MatrixXd(N_LIDAR_MEAS, N_LIDAR_MEAS);
      R << pow(std_laspx_, 2), 0,
          0, pow(std_laspy_, 2);

      // predicted state covariance matrix
      S.fill(0.0);
      for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
        // state difference
        VectorXd z_diff = Zsig.col(i) - z_pred;
        // angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

        S = S + weights(i) * z_diff * z_diff.transpose() ;
      }
      
      S = S + R;
    std::cout << "==========CP4" << std::endl;
      // write result
      z_pred_lidar_ = z_pred;
      S_lidar_ = S;
      Zsig_lidar_ = Zsig;

      std::cout << "z_pred_lidar_\n" << z_pred_lidar_ << std::endl;
      std::cout << "S_lidar_\n" << S_lidar_ << std::endl;
      std::cout << "Zsig_lidar_\n" << Zsig_lidar_ << std::endl;

      UpdateLidar(meas_package);
    }
  }
  else
  {
    std::cout << "==========CP2.0" << std::endl;

    if (is_initialized_ == false)
    {
      std::cout << "==========CP2.1" << std::endl;
      double rho = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
      double rho_d = meas_package.raw_measurements_(2);
      double p_x = rho * cos(phi);
      double p_y = rho * sin(phi);
      // To note: this is inaccurate value, aim to initialize velocity which's magnitude/order is close to real value
      double vx = rho_d * cos(phi);
      double vy = rho_d * sin(phi);
      double v = sqrt(vx * vx + vy * vy);
      x_ << p_x, p_y, v,0, 0;

      P_ << 0.5, 0, 0, 0, 0,
            0, 0.5, 0, 0, 0,
            0, 0, 0.75, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;

      // update time
      time_us_ = meas_package.timestamp_;

      is_initialized_ = true;
    }
    else
    {
      // here we need to update hte timestamp. The timestamp origin unit is us. We need to map to s.
      double dt = double(meas_package.timestamp_ - time_us_) / 1e6;
      time_us_ = meas_package.timestamp_;
      Prediction(dt);

      std::cout << "==========CP2.2" << std::endl;
      // do radar process
      // create matrix for sigma points in measurement space
      MatrixXd Zsig = MatrixXd(N_RADAR_MEAS, 2 * n_aug_ + 1);

      // mean predicted measurement
      VectorXd z_pred = VectorXd(N_RADAR_MEAS);
      
      // measurement covariance matrix S
      MatrixXd S = MatrixXd(N_RADAR_MEAS, N_RADAR_MEAS);

      // transform sigma points into measurement space
      for (int i = 0; i < 2 * n_aug_ + 1; i ++)
      {
          double px = Xsig_pred_(0, i);
          double py = Xsig_pred_(1, i);
          double v = Xsig_pred_(2, i);
          double yaw = Xsig_pred_(3, i);

          Zsig(0, i) = sqrt(pow(px, 2) + pow(py, 2));
          Zsig(1, i) = atan2(py, px);
          Zsig(2, i) = (px * cos(yaw) * v + py * sin(yaw) * v) / Zsig(0, i);
      }
      std::cout << "==========CP3" << std::endl;
      std::cout << "Zsig is:\n" << Zsig << std::endl;
      // calculate mean predicted measurement
      z_pred = Zsig * weights;

      // calculate innovation covariance matrix S
      MatrixXd R = MatrixXd(N_RADAR_MEAS, N_RADAR_MEAS);
      R << pow(std_radr_, 2), 0, 0,
          0, pow(std_radphi_, 2), 0,
          0, 0, pow(std_radrd_, 2);

      // predicted state covariance matrix
      S.fill(0.0);
      for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
        // state difference
        VectorXd z_diff = Zsig.col(i) - z_pred;
        // angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

        S = S + weights(i) * z_diff * z_diff.transpose() ;
      }
      
      S = S + R;
    std::cout << "==========CP4" << std::endl;
      // write result
      z_pred_radar_ = z_pred;
      S_radar_ = S;
      Zsig_radar_ = Zsig;
      std::cout << "z_pred_radar\n" << z_pred_radar_ << std::endl;
      std::cout << "S_radar_\n" << S_radar_ << std::endl;
      std::cout << "Zsig_radar_\n" << Zsig_radar_ << std::endl;

      UpdateRadar(meas_package);
    }
  }
}

void UKF::Prediction(double delta_t) {

  std::cout << "==========CP5" << std::endl;
  std::cout << "delta_t is:\n" << delta_t << std::endl;
  // Step 1, predict sigma points with augumentation
  // create augmented mean state
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.fill(0.0);
  Xsig_aug.col(0) << x_,
                     0,
                     0;
  std::cout << "==========CP5.1" << std::endl;
  P_aug_.row(0) << P_.row(0), 0, 0;
  P_aug_.row(1) << P_.row(1), 0, 0;
  P_aug_.row(2) << P_.row(2), 0, 0;
  P_aug_.row(3) << P_.row(3), 0, 0;
  P_aug_.row(4) << P_.row(4), 0, 0;
  P_aug_.row(5) << 0, 0, 0, 0, 0, std_a_ * std_a_, 0;
  P_aug_.row(6) << 0, 0, 0, 0, 0, 0, std_yawdd_ * std_yawdd_;

  std::cout << "==========CP5.2" << std::endl;
  std::cout << "P_aug_\n" << P_aug_ << std::endl;

  // create augmented covariance matrix
  MatrixXd A = P_aug_.llt().matrixL();
  std::cout << "==========CP6" << std::endl;
  // create augmented sigma points
  for (int i = 0; i < n_aug_; i ++)
  {
      Xsig_aug.col(i + 1) << Xsig_aug.col(0) + 1.73205 * A.col(i);
      Xsig_aug.col(i + 1 + n_aug_) << Xsig_aug.col(0) - 1.73205 * A.col(i);
  }

  // Step 2. Sigma points prediction
  // predict sigma points
  for (int i = 0; i< 2 * n_aug_ + 1; ++i) {
    // extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    // predicted state values
    double px_p, py_p;

    // avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    } else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    // add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    // write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
  std::cout << "==========CP7" << std::endl;
  std::cout << "Xsig_pred_\n" << Xsig_pred_ << std::endl;

  // Step 3. Update mean and covariance
  // create vector for predicted state
  VectorXd x = VectorXd(n_x_);
  MatrixXd temp = Xsig_pred_ * weights;
  x = temp.col(0);
  std::cout << "==========CP7.0" << std::endl;
  std::cout << "temp\n" << temp << std::endl;
  std::cout << "x\n" << x << std::endl;

  // predict state covariance matrix
  MatrixXd tempCov = Xsig_pred_;
  for (int j = 0; j < n_x_; j ++)
  {
      for (int i = 0; i < 2 * n_aug_ + 1 ; i ++)
      {
          tempCov(j, i) = (tempCov(j, i) - x(j));
      }
  }

  // predicted state covariance matrix
  MatrixXd P = MatrixXd(n_x_, n_x_);
  P.fill(0.0);
  std::cout << "==========CP7.1" << std::endl;
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P = P + weights(i) * x_diff * x_diff.transpose() ;
  }
  std::cout << "==========CP8" << std::endl;
  // write result
  x_ = x;
  P_ = P;

  std::cout << "x_\n" << x_ << std::endl;
  std::cout << "P_\n" << P_ << std::endl;
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, N_LIDAR_MEAS);
  std::cout << "==========CP9.0.0" << std::endl;

  // calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    // state difference
    VectorXd z_diff = Zsig_lidar_.col(i) - z_pred_lidar_;
    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
  }

  // calculate Kalman gain K;
  MatrixXd Kgain = MatrixXd(n_x_, N_LIDAR_MEAS);
  Kgain = Tc * S_lidar_.inverse();
  std::cout << "==========CP10.0.0" << std::endl;

  // residual
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred_lidar_;

  // angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  // update state mean and covariance matrix
  P_ = P_ - Kgain * S_lidar_ * Kgain.transpose();
  x_ = x_ + Kgain * z_diff;

  double NIS = z_diff.transpose() * S_lidar_.inverse() * z_diff;

  if (NIS > 7.815)
  {
    overcnt_ ++;
  }

  allcnt_ ++;

  std::cout << "NIS over threashold percentage:\n" << (double)overcnt_ / allcnt_ << std::endl;
  std::cout << "updated z_diff:\n" << z_diff << std::endl;
  std::cout << "updated x_:\n" << x_ << std::endl;

  //sleep(2);
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, N_RADAR_MEAS);
  std::cout << "==========CP9" << std::endl;
  // calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    // state difference
    VectorXd z_diff = Zsig_radar_.col(i) - z_pred_radar_;
    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
  }
  // calculate Kalman gain K;
  MatrixXd Kgain = MatrixXd(n_x_, N_RADAR_MEAS);
  Kgain = Tc * S_radar_.inverse();
  std::cout << "==========CP10" << std::endl;

  // residual
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred_radar_;

  // angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  // update state mean and covariance matrix
  P_ = P_ - Kgain * S_radar_ * Kgain.transpose();
  x_ = x_ + Kgain * z_diff;

  double NIS = z_diff.transpose() * S_radar_.inverse() * z_diff;

  if (NIS > 7.815)
  {
    overcnt_ ++;
  }

  allcnt_ ++;

  std::cout << "NIS over threashold percentage:\n" << (double)overcnt_ / allcnt_ << std::endl;
  std::cout << "updated z_diff:\n" << z_diff << std::endl;
  std::cout << "updated x_:\n" << x_ << std::endl;

  //sleep(2);
}

void UKF::Update(MeasurementPackage meas_package)
{
  // if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER)
  // {
  //   // update for Lidar
  //   UpdateLidar(meas_package);
  // }
  // else if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR)
  // {
  //   // update for Radar
  //   UpdateRadar(meas_package);
  // }
}