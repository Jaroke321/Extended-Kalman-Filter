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

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  H_laser_ << 1, 0, 0, 0,
  			 0, 1, 0, 0;
  
  //state transition matrix
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
  			 0, 1, 0, 1,
  			 0, 0, 1, 0,
  			 0, 0, 0, 1;
  
  noise_ay = 9.0;
  noise_ax = 9.0;

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
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    //cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    //state covariance matrix
  	ekf_.P_ = MatrixXd(4, 4);
 	ekf_.P_ << 1, 0, 0, 0, //revert back to 1 when stuck
  		   	   0, 1, 0, 0, //revert back to 1 when stuck
  		   	   0, 0, 1000, 0, //revert back to 1000 when stuck
  		   	   0, 0, 0, 1000; //revert back to 1000 when stuck
    previous_timestamp_ = measurement_pack.timestamp_;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float r = measurement_pack.raw_measurements_[0];
      float theta = measurement_pack.raw_measurements_[1];
      
      float x = r * cos(theta);
      float y = r * sin(theta);
      
      ekf_.x_ << x, y, 1, 1;
      
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 1, 1;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;
  
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt_4/4.0*noise_ax, 0, dt_3/2.0*noise_ax, 0,
			   0, dt_4/4.0*noise_ay, 0, dt_3/2.0*noise_ay,
			   dt_3/2.0*noise_ax, 0, dt_2*noise_ax, 0,
			   0, dt_3/2.0*noise_ay, 0, dt_2*noise_ay;
    
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    ekf_.H_ = MatrixXd(3, 4);
    ekf_.H_ << tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = MatrixXd(3, 3);
    ekf_.R_ << R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_); // Radar updates....Extended Kalman filter equations...Uses jacobian matrix
  } else {
    // Laser updates....Standard Kalman filter equations...doesnt use jacobian
    ekf_.H_ = MatrixXd(2, 4);
    ekf_.H_ << H_laser_;
    ekf_.R_ = MatrixXd(2, 2);
    ekf_.R_ << R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}