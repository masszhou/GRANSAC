#pragma once

#include "AbstractModel.hpp"
#include "basic_types.hpp"
#include <map>

namespace GRANSAC
{

class Line2DModel: public AbstractModel<2>
{
protected:
	// Parametric form
	float m_a, m_b, m_c; // ax + by + c = 0
	float m_dist_denominator; // = sqrt(a^2 + b^2). Stored for efficiency reasons

	// Another parametrization y = mx + d
	float m_m; // Slope
	float m_d; // Intercept

	virtual float  computeDistanceMeasure(std::shared_ptr<AbstractParameter> input_data) override
	{
		auto ext_point2D = std::dynamic_pointer_cast<Point2D>(input_data);
		if (ext_point2D == nullptr)
			throw std::runtime_error("Line2DModel::ComputeDistanceMeasure() - Passed parameter are not of type Point2D.");

		// Return distance between passed "point" and this line
		// http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
		float numer = fabs(m_a * ext_point2D->m_point2D[0] + m_b * ext_point2D->m_point2D[1] + m_c);
		float dist = numer / m_dist_denominator;

		return dist;
	};

public:
	Line2DModel(const std::vector<std::shared_ptr<AbstractParameter>> &input_data, 
	            const std::map<std::string, float>& additional_params)
	{
		initialize(input_data, additional_params);
	};

	virtual void initialize(const std::vector<std::shared_ptr<AbstractParameter>> &input_data, 
	                        const std::map<std::string, float>& additional_params) override
	{
		if (input_data.size() != 2)
			throw std::runtime_error("Line2DModel - Number of input parameters does not match minimum number required for this model.");

		// Check for AbstractParamter types
		auto point1 = std::dynamic_pointer_cast<Point2D>(input_data[0]);
		auto point2 = std::dynamic_pointer_cast<Point2D>(input_data[1]);
		if (point1 == nullptr || point2 == nullptr)
			throw std::runtime_error("Line2DModel - InputParams type mismatch. It is not a Point2D.");

		m_model_def_parameters = input_data;

		// Compute the line parameters
		m_m = (point2->m_point2D[1] - point1->m_point2D[1]) / (point2->m_point2D[0] - point1->m_point2D[0]); // Slope
		m_d = point1->m_point2D[1] - m_m * point1->m_point2D[0]; // Intercept

		// mx - y + d = 0
		m_a = m_m;
		m_b = -1.0;
		m_c = m_d;

		m_dist_denominator = sqrt(m_a * m_a + m_b * m_b); // Cache square root for efficiency
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