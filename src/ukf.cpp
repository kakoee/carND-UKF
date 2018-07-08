#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
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
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_ ;
  //set weights
  weights_ = VectorXd(2*n_aug_+1);
  weights_.fill(1/(2*(lambda_ + n_aug_)));
  weights_[0] = lambda_ / (lambda_+ n_aug_);
  
  // first measurement
  x_ << 1, 1, 1, 1, 0.1;

  // init covariance matrix
  P_ << 0.1,    0, 0, 0, 0,
               0, 0.1, 0, 0, 0,
               0,    0, 1, 0, 0,
               0,    0, 0, 1, 0,
               0,    0, 0, 0, 1;
          
 Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
 time_us_ = 0.0;

 // new variables

 //add measurement noise covariance matrix
  R_laser = MatrixXd(2,2);
  R_laser <<    std_laspx_*std_laspx_, 0,
          0, std_laspy_*std_laspy_;

  R_radar = MatrixXd(3,3);
  R_radar <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;


}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
/*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
     // first measurement

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      //need to convert polar coordinates to cartesian coordinates
      float x0 = meas_package.raw_measurements_[0] * cos(meas_package.raw_measurements_[1]);
      float x1 = meas_package.raw_measurements_[0] * sin(meas_package.raw_measurements_[1]);
      x_(0) = x0;
      x_(1) = x1;
      x_(2) = 1;
      x_(3) = 1;
      x_(4) = 0.1;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
   		x_(0) = meas_package.raw_measurements_[0];
      x_(1) = meas_package.raw_measurements_[1];
      x_(2) = 1;
      x_(3) = 1;
      x_(4) = 0.1;
     	
        
    }
    time_us_ = meas_package.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

	//compute the time elapsed between the current and previous measurements
	double dt = (meas_package.timestamp_ - time_us_)/1000000.0;	//dt in second
    
	time_us_ = meas_package.timestamp_;

  

  Prediction(dt);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // Radar updates
    
      UpdateRadar(meas_package);
     
    
    } else {
      UpdateLidar(meas_package);
    
    }

}

// newly added funciton to calculate sigmaprediciton

void UKF::CalculatePrediction(std::vector<double>& sp,double dt,std::vector<double>& spred){

    if(sp.size()!=7)
        return   ;    

    double px= sp[0];
    double py= sp[1];
    double v=sp[2];
    double sy = sp[3];
    double sydot=sp[4];
    double v_ak = sp[5];
    double v_sik = sp[6];
    
    double px2,py2,v2,sy2,sydot2;
    if(fabs(sydot)>0.001){
        px2= px+ (v/sydot)*(sin(sy+dt*sydot)-sin(sy))   ;
        py2= py + (v/sydot)*(-cos(sy+dt*sydot)+cos(sy));
        v2=v+0;
        sy2 = sy + dt*sydot;
        sydot2=sydot+0;
    }
    else{
        px2= px+ (v*cos(sy)*dt)   ;
        py2= py + (v*sin(v)*dt);
        v2=v+0;
        sy2 = sy ;
        sydot2=sydot+0;
    }
    px2= px2 + 0.5*(dt*dt) * cos(sy)* v_ak;
    py2 = py2+ 0.5 *(dt*dt) * sin(sy)*v_ak;
    v2 = v2 + dt *v_ak;
    sy2 = sy2  + 0.5 * (dt*dt) * v_sik;
    sydot2 = sydot2 + dt * v_sik;
    spred[0]=(px2);
    spred[1]=(py2);
    spred[2]=(v2);
    spred[3]=(sy2);
    spred[4]=(sydot2);

    return ;
    
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

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
 
  //create augmented mean state
  x_aug.head(n_x_)  = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;
  //create augmented covariance matrix;
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug.bottomRows(2) = MatrixXd::Zero(2, n_aug_);
  P_aug.rightCols(2) = MatrixXd::Zero(n_aug_,2);
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;
  //create square root matrix

  MatrixXd A = P_aug.llt().matrixL();
  //create augmented sigma points
  Xsig_aug.col(0) =x_aug;
  for(int i=1;i<n_aug_+1;i++)
  {
      Xsig_aug.col(i) = (x_aug + A.col(i-1) * sqrt(lambda_+n_aug_));
  }
  for(int i=n_aug_+1;i<2*n_aug_+1;i++)
  {
      Xsig_aug.col(i) = (x_aug - A.col(i-n_aug_-1) * sqrt(lambda_+n_aug_));
  }

 /*******************************************************************************
 * predict sigma points
 ******************************************************************************/
  std::vector<double> spoint(n_aug_);
  std::vector<double> spoint_pred(n_x_);
  for(int i = 0;i<2*n_aug_+1;i++){
      for(int j=0;j<n_aug_;j++){
         spoint[j]=Xsig_aug(j,i);
      }
      CalculatePrediction(spoint,delta_t,spoint_pred);
      for(int j=0;j<n_x_;j++){
          Xsig_pred_(j,i) = spoint_pred[j];
      }  
  } 

 /*******************************************************************************
 * Predict Mean and Covariance 
 ******************************************************************************/
  

  //predict state mean
  double t=0;
  for(int i=0;i<n_x_;i++)
  {
      t=0;
      for(int j=0;j<2*n_aug_+1;j++){
          t = t + (weights_[j] * Xsig_pred_(i,j));
      }
      x_(i) = t;
      //std::cout<<"test**************X"<< i << " = "<<x_[i]<<std::endl;
  }
  //predict state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI){
      x_diff(3)-=2.*M_PI;
      //std::cout<<"test**************"<<x_diff(3)<<std::endl;
    } 
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }

  
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  /*******************************************************************************
    * Predict Lidar Measurement
   *******************************************************************************/
  int n_z=2;
  VectorXd z;
  z = meas_package.raw_measurements_;
  //create matrix for sigma points in radar measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  
    //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
   
    Zsig(0,i) = p_x;                        //r
    Zsig(1,i) = p_y;                                 //phi
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

     S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  
  S = S + R_laser;
    
  /*******************************************************************************
    * Update Radar Measurement
   *******************************************************************************/

  
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }


  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;


  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  //calculate NIS
  laser_nis_ = z_diff.transpose() * S.inverse() * z_diff;
  // NIS will be printed in main.cpp in msgjson
  //std::cout<< "***Laser NIS=" << laser_nis_  <<std::endl; 
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  
  /*******************************************************************************
    * Predict Radar Measurement
   *******************************************************************************/

  //create matrix for sigma points in radar measurement space
  VectorXd z;
  z = meas_package.raw_measurements_;

  int n_z = 3;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  
    //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    //while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    //while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    z_diff(1) = atan2(sin(z_diff(1)),cos(z_diff(1)));

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  
  S = S + R_radar;
    
 
  /*******************************************************************************
    * Update Radar Measurement
   *******************************************************************************/
  
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    //while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    //while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    z_diff(1) = atan2(sin(z_diff(1)),cos(z_diff(1)));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    //while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    //while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
    x_diff(3) = atan2(sin(x_diff(3)),cos(x_diff(3)));


    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  //while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  //while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  z_diff(1) = atan2(sin(z_diff(1)),cos(z_diff(1)));

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  radar_nis_ = z_diff.transpose() * S.inverse() * z_diff;

  // NIS wil lbe printed in main.cpp in msgjson
  //std::cout<< "radar NIS=" << radar_nis_  <<std::endl; 

    
}
