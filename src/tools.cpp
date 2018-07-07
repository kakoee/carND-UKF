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
	rmse<<0,0,0,0;
    

	if(estimations.size()==0 || (estimations.size()!=ground_truth.size())){
		cout << "Invalid estimation or ground_truth data" << endl;
		return rmse;
    }


	for(unsigned int i=0; i < estimations.size(); ++i){
        VectorXd item_ext = estimations[i];
        VectorXd item_real = ground_truth[i];
        for(int j=0;j<item_ext.size();j++){
            float sub = item_ext(j) - item_real(j);            
            rmse(j) = rmse(j) + sub*sub;
        }
	}


	rmse = rmse/estimations.size();

	rmse = rmse.array().sqrt(); 

	//return the result  
    return rmse;
}