//
// Created by ChenYL on 2022/8/1.
// Reference:https://www.guyuehome.com/15356
//

#ifndef KF_H
#define KF_H
#include "Eigen/Dense"


class KalmanFilter
{
private:
    int stateSize; //state variable's dimenssion
    int measSize; //measurement variable's dimession
    int uSize; //control variables's dimenssion
    Eigen::VectorXd x; //状态
    Eigen::VectorXd z; //测量值

    Eigen::MatrixXd A; //状态转移矩阵
    Eigen::MatrixXd B; //输入控制矩阵

    Eigen::VectorXd u; //控制向量，外界对系统的作用

    Eigen::MatrixXd P; //coveriance
    Eigen::MatrixXd H; //观测矩阵
    Eigen::MatrixXd R; //measurement noise covariance
    Eigen::MatrixXd Q; //process noise covariance

    void CalculateJacobianMatrix();

public:
    KalmanFilter(int stateSize_=0, int measSize_=0,int uSize_=0);//构造函数
    ~KalmanFilter(){} //析构函数
    void init(Eigen::VectorXd &x_, Eigen::MatrixXd& P_,Eigen::MatrixXd& R_, Eigen::MatrixXd& Q_);
    Eigen::VectorXd predict(Eigen::MatrixXd& A_);
    Eigen::VectorXd predict(Eigen::MatrixXd& A_, Eigen::MatrixXd &B_, Eigen::VectorXd &u_);
    void KFUpdate(Eigen::MatrixXd& H_, Eigen::VectorXd z_meas);
    void EKFUpdate(Eigen::MatrixXd& H_, Eigen::VectorXd z_meas);
};

#endif //KF_H
