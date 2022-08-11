//
// Created by ChenYL on 2022/8/1.
// Reference:https://www.guyuehome.com/15356
//
#include <iostream>
#include "KF.h"
using namespace std;



KalmanFilter::KalmanFilter(int stateSize_, int measSize_, int uSize_) :stateSize(stateSize_), measSize(measSize_), uSize(uSize_)
{
    if (stateSize == 0 || measSize == 0)
    {
        std::cerr << "Error, State size and measurement size must bigger than 0\n";
    }

    x.resize(stateSize);
    x.setZero();

    A.resize(stateSize, stateSize);
    A.setIdentity();

    u.resize(uSize);
    u.transpose();
    u.setZero();

    B.resize(stateSize, uSize);
    B.setZero();

    P.resize(stateSize, stateSize);
    P.setIdentity();

    H.resize(measSize, stateSize);
    H.setZero();

    z.resize(measSize);
    z.setZero();

    Q.resize(stateSize, stateSize);
    Q.setZero();

    R.resize(measSize, measSize);
    R.setZero();
}

void KalmanFilter::init(Eigen::VectorXd &x_, Eigen::MatrixXd& P_, Eigen::MatrixXd& R_, Eigen::MatrixXd& Q_)
{
    x = x_;
    P = P_;
    R = R_;
    Q = Q_;
}

// 已知控制矩阵B和控制向量u
Eigen::VectorXd KalmanFilter::predict(Eigen::MatrixXd& A_, Eigen::MatrixXd &B_, Eigen::VectorXd &u_)
{
    A = A_;
    B = B_;
    u = u_;
    x = A * x + B * u;
    Eigen::MatrixXd F_T = A.transpose();
    P = A * P * F_T + Q;
    return x;
}

// 未知控制矩阵B和控制向量u
Eigen::VectorXd KalmanFilter::predict(Eigen::MatrixXd& A_)
{
    A = A_;
    x = A * x;
    Eigen::MatrixXd A_T = A.transpose();
    P = A * P * A_T + Q;
//    cout << "P-=" << P<< endl;
    return x;
}

void KalmanFilter::KFUpdate(Eigen::MatrixXd& H_, Eigen::VectorXd z_meas)
{
    H = H_;
    Eigen::MatrixXd S,Ht;
    Ht = H.transpose();
    S = H * P * Ht + R;
    Eigen::MatrixXd K = P*Ht*S.inverse(); //卡尔曼增益
    z = H * x;
    x = x + K * (z_meas - z);
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(stateSize, stateSize);
    P = (I - K * H) * P;
//    cout << "P=" << P << endl;
}

void KalmanFilter::EKFUpdate(Eigen::MatrixXd& H_,  Eigen::VectorXd z)
{
    H = H_;
    // 根据实际情况修改
    double rho = sqrt(x(0)*x(0) + x(1)*x(1));
    double theta = atan2(x(1), x(0));
    double rho_dot = (x(0)*x(2) + x(1)*x(3)) / rho;
    Eigen::VectorXd h = Eigen::VectorXd(3);
    h << rho, theta, rho_dot;
    Eigen::VectorXd y = z - h;

    CalculateJacobianMatrix();

    Eigen::MatrixXd Ht = H.transpose();
    Eigen::MatrixXd S = H * P * Ht + R;
    Eigen::MatrixXd Si = S.inverse();
    Eigen::MatrixXd K =  P * Ht * Si;
    x = x + (K * y);
    int x_size = x.size();
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(x_size, x_size);
    P = (I - K * H) * P;
}

void KalmanFilter::CalculateJacobianMatrix()
{
    //measSize=3, stateSize=4的情况
    Eigen::MatrixXd Hj(3, 4);

    // get state parameters

    float px = x(0);
    float py = x(1);
    float vx = x(2);
    float vy = x(3);

    // pre-compute a set of terms to avoid repeated calculation
    float c1 = px * px + py * py;
    float c2 = sqrt(c1);
    float c3 = (c1 * c2);

    // Check division by zero
    if(fabs(c1) < 0.0001){
        H = Hj;
        return;
    }

    Hj << (px/c2), (py/c2), 0, 0,
            -(py/c1), (px/c1), 0, 0,
            py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;
    H = Hj;
}