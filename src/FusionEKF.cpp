#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  // measurement covariance matrix - laser
  H_laser_ <<
      1, 0, 0, 0,
      0, 1, 0, 0;

  // measurement covariance matrix - radar - must do later because it changes with each measurement

  // transition matrix F_
  ekf_.F_ = MatrixXd(4,4);
  ekf_.F_ << 1,0,1,0,
             0,1,0,1,
             0,0,1,0,
             0,0,0,1;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    cout << "Initializing" << endl;
    ekf_.x_ = VectorXd(4);

    // Initialize state covariance with low confidence in position and very low confidence in velocity.
    // Ideally vx and vy should have a lot of covariance in radar case, but not dealing with that complication yet.
    ekf_.P_ = MatrixXd(4,4);
    ekf_.P_ << 1,0,0,0,
               0,1,0,0,
               0,0,1000,0,
               0,0,0,1000;

    previous_timestamp_ = measurement_pack.timestamp_;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      double ro = measurement_pack.raw_measurements_[0];
      double theta = measurement_pack.raw_measurements_[1];
      double ro_dot = measurement_pack.raw_measurements_[2];

      // Convert radar from polar to cartesian coordinates
      double px = ro * sin(theta);
      double py = ro * cos(theta);

      // Low estimates of velocity, assuming no motion perpendicular to radar beam.
      double vx = ro_dot * sin(theta);
      double vy = ro_dot * cos(theta);

      ekf_.x_ << px, py, vx, vy;
      cout << "Initializing with radar: " << ekf_.x_ << endl;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      double px = measurement_pack.raw_measurements_[0];
      double py = measurement_pack.raw_measurements_[1];

      ekf_.x_ << px, py, 0, 0;
      cout << "Initializing with lidar: " << ekf_.x_ << endl;
    }

    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  // Convert time elapsed from microseconds to seconds
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;

  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;

  double noise_ax = 9;
  double noise_ay = 9;

  ekf_.Q_ <<
      0.25 * noise_ax * pow(dt,4), 0, 0.5 * noise_ax * pow(dt,3), 0,
      0, 0.25 * noise_ay * pow(dt,4), 0, 0.5 * noise_ay * pow(dt,3),
      0.5 * noise_ax * pow(dt,3), 0, noise_ax * pow(dt,2), 0,
      0, 0.5 * noise_ay * pow(dt,3), 0, noise_ay * pow(dt,2);

  ekf_.Predict();

  cout << "After prediction with t=" << dt << endl;
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
  cout << endl;

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    cout << "Updating with radar" << measurement_pack.raw_measurements_ << endl;
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    cout << "Updating with lidar" << measurement_pack.raw_measurements_ << endl;
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  previous_timestamp_ = measurement_pack.timestamp_;

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
  cout << endl;
}
