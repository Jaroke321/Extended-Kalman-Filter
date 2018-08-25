#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;
  
  if(estimations.size() != ground_truth.size() || estimations.size() == 0)
  {
    std::cout<<"invalid estimations vector or ground truth";
    return rmse;
  }
  
  for(unsigned int i; i < estimations.size(); ++i){
    
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse += residual;
  }
  
  //calculate mean
  rmse = rmse.array() / estimations.size();
  
  //calculate the squared root
  rmse = rmse.array().sqrt();
  
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3, 4);
  
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);
  
  float c1 = px*px + py*py;
  float c2 = sqrt(c1);
  //float c3 = (c1*c2);
  
  if(fabs(c1) < 0.0001){
    std::cout<<"Jacobian matrix error, division by zero";
    return Hj;
  }
  
  Hj << (px/c2), (py/c2), 0, 0,
  		-(py/c1), (px/c1), 0, 0,
  		py*((vx*py)-(vy*px))/pow(c1, 1.5), px*((vy*px)-(vx*py))/pow(c1, 1.5), px/c2, py/c2;
  
  return Hj;
  
}
