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
  rmse << 0, 0, 0, 0;
  if (estimations.size() == 0) {
    return rmse;
  }
  if (estimations.size() != ground_truth.size()) {
    throw "estimations and ground_truth are different length";
  }
  for(int i = 0; i < estimations.size(); i++) {
    VectorXd residual = estimations[i] - ground_truth[i];
    VectorXd residual_squared = residual.array() * residual.array();
    rmse += residual_squared;
  }
  rmse /= estimations.size();
  rmse = rmse.array().sqrt();
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj(3,4);

  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);

  double pxy = pow(px,2) + pow(py,2);
  if(pxy < 0.00001) {
    return Hj;
  }
  double pxy_sqrt = pow(pxy,0.5);
  double pxy_three_halves = pow(pxy,1.5);
  double z = (vx * py - vy * px) / pxy_three_halves;

  Hj << (px / pxy_sqrt), (py / pxy_sqrt), 0, 0,
        -(py / pxy), (px / pxy), 0, 0,
        (py * z), -(px * z), (px / pxy_sqrt), (py / pxy_sqrt);

  return Hj;
}
