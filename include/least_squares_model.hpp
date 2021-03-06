#ifndef LINEAR_LEAST_SQUARE_MODEL_HPP
#define LINEAR_LEAST_SQUARE_MODEL_HPP

#include <iostream>

#include <math.h>
#include <Eigen/QR>
#include <Eigen/Dense>

#include "basic_types.hpp"
#include "regression2d.hpp"


namespace GRANSAC
{

// model paramter number = 3
template<int t_param_num>
class LinearLeastSquaresModel: public AbstractModel<t_param_num>
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
	int m_grid_size_x; // cell size of grid, e.g. m_grid_size = 10 -> each cell is 10x10 pixel
	int m_grid_size_y;
	std::vector<std::vector<float> > m_occupied_list;

	virtual float computeDistanceMeasure(std::shared_ptr<AbstractParameter> input_data) override
	{
		auto ext_point2D = std::dynamic_pointer_cast<Point2D>(input_data);
		if (ext_point2D == nullptr)
			throw std::runtime_error("PolynomialModel::ComputeDistanceMeasure() - Passed parameter are not of type Point2D.");

		// build a lookup table for distance calculation
		float min_dist = std::numeric_limits<float>::max();
        for (auto each : m_occupied_list){
            float dist = fabs(ext_point2D->m_point2D[0]/m_grid_size_x-each[0]) + fabs(ext_point2D->m_point2D[1]/m_grid_size_y-each[1]); // p-1 distance, 10 is grid size
            if (min_dist > dist)
                min_dist = dist;
        }
		return min_dist; // distance in grid, 10 times smaller than actual distance 
	};

public:
	LinearLeastSquaresModel(const std::vector<std::shared_ptr<AbstractParameter>> &input_data, 
	               const std::map<std::string, float>& additional_params)
	{
		initialize(input_data, additional_params);
	};

	virtual void initialize(const std::vector<std::shared_ptr<AbstractParameter>> &input_data, 
	                        const std::map<std::string, float>& additional_params) override
	{
		int img_width = (additional_params.count("img_width") == 1) ? additional_params.at("img_width") : 400;
		int img_height = (additional_params.count("img_height") == 1) ? additional_params.at("img_height") : 400;
		int grid_num_x = (additional_params.count("grid_num_x") == 1) ? additional_params.at("grid_num_x") : 40;
		int grid_num_y = (additional_params.count("grid_num_y") == 1) ? additional_params.at("grid_num_y") : 40;

		m_grid_size_x = img_width / grid_num_x; // e.g. 10
		m_grid_size_y = img_height / grid_num_y;

		// alway calculate curve with three points, since y = ax^2+bx+c
		if (input_data.size() < t_param_num)
			throw std::runtime_error("PolynomialModel - Number of input parameters does not match minimum number required for this model.");

		AbstractModel<t_param_num>::m_model_def_parameters = input_data;

		// compute deterministic curve parameters with 3 points
		std::vector<float> x_values, y_values, coeff;
		for (int i=0; i< input_data.size(); i++){
			auto point = std::dynamic_pointer_cast<Point2D>(input_data[i]);
			if (point == nullptr)
				throw std::runtime_error("QuadraticModel - InputParams type mismatch. It is not a Point2D.");

			x_values.push_back(point->m_point2D[0]);
			y_values.push_back(point->m_point2D[1]);
		}
		Regression2D::calculateLLS(x_values, y_values, coeff, t_param_num-1);

		std::vector<float> coeff_f(coeff.begin(), coeff.end());
		AbstractModel<t_param_num>::m_model_coeffs = coeff_f;
		
		m_occupied_list.reserve(grid_num_x*grid_num_y/2);
		//e.g. img_size=400x400, grid_size=40x40, cell_size=10x10
		for (int i=0; i<grid_num_x; i++){ // e.g. i=0; i<40
			
			float x_0 = float(i*m_grid_size_x); // e.g. i*10
			float x_1 = float((i+1)*m_grid_size_x);

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

				int y_0_idx = floor(y_0 / m_grid_size_y); // e.g. y_0_idx = floor(y_0 / 10)
				int y_1_idx = floor(y_1 / m_grid_size_y);

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

	virtual std::pair<float, std::vector<std::shared_ptr<AbstractParameter> > > evaluate(const std::vector<std::shared_ptr<AbstractParameter>>& evaluate_data, float threshold)
	{
		std::vector<std::shared_ptr<AbstractParameter>> inliers;
		int n_total_data = evaluate_data.size();
		int n_inliers = 0;

		for (auto& each : evaluate_data)
		{
			if (computeDistanceMeasure(each) < threshold)
			{
				inliers.push_back(each);
				n_inliers++;
			}
		}
		float inlier_fraction = float(n_inliers) / float(n_total_data); // This is the inlier fraction
		return std::make_pair(inlier_fraction, inliers);
	};

};


}

#endif /* LINEAR_LEAST_SQUARE_MODEL_HPP */
