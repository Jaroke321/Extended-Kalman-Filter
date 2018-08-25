#include "kalman_filter.h"

#define _USE_MATH_DEFINES
#include <math.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  //get individual values from predicted state
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  //shorten repeated terms in h(x)
  float px_2 = px * px;
  float py_2 = py * py;
  float p_sqrt = sqrt(px_2 + py_2);
  float phi = atan2(py, px);
  
  //create h(x) function and apply it to predicted values to convert cartesian values to polar coordinates
  VectorXd Hx;
  Hx = VectorXd(3);
  Hx << p_sqrt,
  		phi,
  		((px*vx) + (py*vy)) / p_sqrt;
  
  VectorXd y = z - Hx;
  //make sure phi in the y Vector is between -pi and pi
  while(y[1] < -M_PI || y[1] > M_PI){
    if(y[1] < -M_PI){
      y[1] += 2.0*M_PI;
    }
    else if(y[1] > M_PI){
      y[1] -= 2.0*M_PI;
    }
  }
  
  //rest of Update EKF algorithm here
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
