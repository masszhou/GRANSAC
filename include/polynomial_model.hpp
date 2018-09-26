#ifndef POLYNOMIAL_MODEL_HPP
#define POLYNOMIAL_MODEL_HPP

#include "AbstractModel.hpp"
#include "basic_types.hpp"
#include <math.h>
#include <Eigen/QR>

#include <iostream>

namespace GRANSAC
{

// model paramter number = 3
template<int t_param_num>
class PolynomialModel: public AbstractModel<t_param_num>
{
protected:    
	//[1  x0  x0^2] [m_a0] = [y0]
    //[1  x1  x1^2] [m_a1]   [y1]
    //[1  x2  x2^2] [m_a2]   [y2]
    //...                    ...
    //[1  xn  xn^2]          [yn]
    // std::vector<VPFloat> m_a;

	// build a lookup table to calculate point to curve distance
	// e.g. target fitting is points in a 400x400 image
	// then grid of lookup table is 40x40 cell, each cell is 10x10 pixel
	int m_grid_size; // cell size of grid, e.g. m_grid_size = 10 -> each cell is 10x10 pixel
	std::vector<std::vector<float> > m_occupied_list;

	virtual VPFloat ComputeDistanceMeasure(std::shared_ptr<AbstractParameter> Param) override
	{
		auto ExtPoint2D = std::dynamic_pointer_cast<Point2D>(Param);
		if (ExtPoint2D == nullptr)
			throw std::runtime_error("PolynomialModel::ComputeDistanceMeasure() - Passed parameter are not of type Point2D.");

		// distance of point (x0, y0) to quadratic curve y(x)=a0+a1*x+a2*x*x
		// http://mathworld.wolfram.com/Point-QuadraticDistance.html
		// solve a cubic equation to find (x',y') on y(x)
		// 2*a2^2*x^3 + 3*a1a2*x^2 + (a1^2 + 2*a0a2 - 2*a2y0 + 1 )x + (a0a1 - a1y0 - x0)=0
		// correct but not efficient, I don't really need a high accuracy distance.
		// build a lookup table for distance calculation
		float min_dist = std::numeric_limits<float>::max();
        for (auto each : m_occupied_list){
            float dist = fabs(ExtPoint2D->m_Point2D[0]/10-each[0]) + fabs(ExtPoint2D->m_Point2D[1]/10-each[1]); // p-1 distance, 10 is grid size
            if (min_dist > dist)
                min_dist = dist;
        }
		return min_dist; // distance in grid, 10 times smaller than actual distance 
	};

public:
	PolynomialModel(const std::vector<std::shared_ptr<AbstractParameter>> &input_params, 
	               const std::vector<float>& additional_params)
	{
		Initialize(input_params, additional_params);
	};

	virtual void Initialize(const std::vector<std::shared_ptr<AbstractParameter>> &InputParams, 
	                        const std::vector<float>& additional_params) override
	{
		int img_width;
		int img_height;
		int grid_num_x;
		int grid_num_y;

		if (additional_params.size() > 0){
			img_width = int(additional_params[0]);
			img_height = int(additional_params[1]);
			grid_num_x = int(additional_params[2]);
			grid_num_y = int(additional_params[2]);
		}else{
			img_width = 400;
			img_height = 400;
			grid_num_x = 40;
			grid_num_y = 40;
			// each grid cell is 10x10
		}
		int grid_size_x = img_width / grid_num_x; // e.g. 10
		int grid_size_y = img_height / grid_num_y;

		// alway calculate curve with three points, since y = ax^2+bx+c
		if (InputParams.size() != t_param_num)
			throw std::runtime_error("PolynomialModel - Number of input parameters does not match minimum number required for this model.");

		// Check for AbstractParamter types
		// auto Point1 = std::dynamic_pointer_cast<Point2D>(InputParams[0]);
		// auto Point2 = std::dynamic_pointer_cast<Point2D>(InputParams[1]);
		// auto Point3 = std::dynamic_pointer_cast<Point2D>(InputParams[2]);
		// if (Point1 == nullptr || Point2 == nullptr || Point3 == nullptr)
		// 	throw std::runtime_error("PolynomialModel - InputParams type mismatch. It is not a Point2D.");

		std::copy(InputParams.begin(), InputParams.end(), AbstractModel<t_param_num>::m_MinModelParams.begin());

		// compute deterministic curve parameters with 3 points
		std::vector<double> x_values, y_values, coeff;
		for (int i=0; i< InputParams.size(); i++){
			auto point = std::dynamic_pointer_cast<Point2D>(InputParams[i]);
			x_values.push_back(point->m_Point2D[0]);
			y_values.push_back(point->m_Point2D[1]);
		}
		polyfit(x_values, y_values, coeff, t_param_num-1);

		//e.g. img_size=400x400, grid_size=40x40, cell_size=10x10
		for (int i=0; i<grid_num_x; i++){ // e.g. i=0; i<40
			
			float x_0 = float(i*grid_size_x); // e.g. i*10
			float x_1 = float((i+1)*grid_size_x);

            float y_0 = 0;
			float y_1 = 0;
            for (int j = 0; j < t_param_num; j++){
                y_0 += coeff[j]*std::pow(x_0, float(j));
                y_1 += coeff[j]*std::pow(x_1, float(j));
            }        

			bool condi_1 = y_0 <= img_height && y_0 >=0; // e.g. condi_1 = y_0<=400 && y_0>=0
			bool condi_2 = y_1 <= img_height && y_1 >=0;
			if (condi_1 || condi_2){
				y_0 = (y_0 > img_height) ? img_height : y_0;
				y_0 = (y_0 < 0) ? 0 : y_0;
				y_1 = (y_1 > img_height) ? img_height : y_1;
				y_1 = (y_1 < 0) ? 0 : y_1;

				int y_0_idx = floor(y_0 / grid_size_y); // e.g. y_0_idx = floor(y_0 / 10)
				int y_1_idx = floor(y_1 / grid_size_y);

				if (y_0_idx < y_1_idx){
					for (int j=y_0_idx; j<y_1_idx; j++){
						std::vector<float> xy_pos{float(i), float(j)};
						m_occupied_list.push_back(xy_pos);
					}
				} else if (y_0_idx > y_1_idx){
					for (int j=y_1_idx; j<y_0_idx; j++){
						std::vector<float> xy_pos{float(i), float(j)};
						m_occupied_list.push_back(xy_pos);
					}
				} else{
					std::vector<float> xy_pos{float(i), float(y_0_idx)};
					m_occupied_list.push_back(xy_pos);
				}
			}
		}
		// //if curve go through cell, fill cell = 1, otherwise = 0
		
	};

	virtual std::pair<VPFloat, std::vector<std::shared_ptr<AbstractParameter>>> Evaluate(const std::vector<std::shared_ptr<AbstractParameter>>& EvaluateParams, VPFloat Threshold)
	{
		std::vector<std::shared_ptr<AbstractParameter>> Inliers;
		int nTotalParams = EvaluateParams.size();
		int nInliers = 0;

		for (auto& Param : EvaluateParams)
		{
			if (ComputeDistanceMeasure(Param) < Threshold)
			{
				Inliers.push_back(Param);
				nInliers++;
			}
		}

		VPFloat InlierFraction = VPFloat(nInliers) / VPFloat(nTotalParams); // This is the inlier fraction

		return std::make_pair(InlierFraction, Inliers);
	};


	static void polyfit(const std::vector<double> &xv, const std::vector<double> &yv, std::vector<double> &coeff, int order)
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

};


}


#endif /* POLYNOMIAL_MODEL_HPP */
