# pragma once

#include <vector>
#include <map>
#include <memory>

#include <math.h>
#include <Eigen/QR>
#include <Eigen/Dense>

using std::vector;

namespace GRANSAC
{

class Regression2D
{
public:

	static void calculateLinearSquaresWithQR(const std::vector<double> &xv, const std::vector<double> &yv, std::vector<double> &coeff, int order)
    {
		Eigen::MatrixXd A(xv.size(), order+1);
		Eigen::VectorXd yv_mapped = Eigen::VectorXd::Map(&yv.front(), yv.size());
		Eigen::VectorXd result;

		assert(xv.size() == yv.size());
		assert(xv.size() >= order+1);

		// create matrix
		for (size_t i = 0; i < xv.size(); i++)
		for (size_t j = 0; j < order+1; j++)
			A(i, j) = pow(xv.at(i), j);

		// solve for linear least squares fit
		result = A.householderQr().solve(yv_mapped);

		coeff.resize(order+1);
		for (size_t i = 0; i < order+1; i++)
			coeff[i] = result[i];
	}

    static void calculateLinearSquaresNormal(const std::vector<double> &xv, const std::vector<double> &yv, std::vector<double> &coeff, int order)
    {
        // faster but with larger condition number compared with QR method
		Eigen::MatrixXd A(xv.size(), order+1);
		Eigen::VectorXd yv_mapped = Eigen::VectorXd::Map(&yv.front(), yv.size());
		Eigen::VectorXd result;

		assert(xv.size() == yv.size());
		assert(xv.size() >= order+1);

		// create matrix
		for (size_t i = 0; i < xv.size(); i++)
		for (size_t j = 0; j < order+1; j++)
			A(i, j) = pow(xv.at(i), j);

		// solve for linear least squares fit
		result = (A.transpose() * A).ldlt().solve(A.transpose() * yv_mapped);

		coeff.resize(order+1);
		for (size_t i = 0; i < order+1; i++)
			coeff[i] = result[i];
	}


};


}
