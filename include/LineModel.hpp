#pragma once

#include "AbstractModel.hpp"
#include "basic_types.hpp"


namespace GRANSAC
{

class Line2DModel: public AbstractModel<2>
{
protected:
	// Parametric form
	VPFloat m_a, m_b, m_c; // ax + by + c = 0
	VPFloat m_DistDenominator; // = sqrt(a^2 + b^2). Stored for efficiency reasons

	// Another parametrization y = mx + d
	VPFloat m_m; // Slope
	VPFloat m_d; // Intercept

	virtual VPFloat ComputeDistanceMeasure(std::shared_ptr<AbstractParameter> Param) override
	{
		auto ExtPoint2D = std::dynamic_pointer_cast<Point2D>(Param);
		if (ExtPoint2D == nullptr)
			throw std::runtime_error("Line2DModel::ComputeDistanceMeasure() - Passed parameter are not of type Point2D.");

		// Return distance between passed "point" and this line
		// http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
		VPFloat Numer = fabs(m_a * ExtPoint2D->m_Point2D[0] + m_b * ExtPoint2D->m_Point2D[1] + m_c);
		VPFloat Dist = Numer / m_DistDenominator;

		//// Debug
		//std::cout << "Point: " << ExtPoint2D->m_Point2D[0] << ", " << ExtPoint2D->m_Point2D[1] << std::endl;
		//std::cout << "Line: " << m_a << " x + " << m_b << " y + "  << m_c << std::endl;
		//std::cout << "Distance: " << Dist << std::endl << std::endl;

		return Dist;
	};

public:
	Line2DModel(const std::vector<std::shared_ptr<AbstractParameter>> &InputParams, 
	            const std::vector<float>& additional_params)
	{
		Initialize(InputParams, additional_params);
	};

	virtual void Initialize(const std::vector<std::shared_ptr<AbstractParameter>> &InputParams, 
	                        const std::vector<float>& additional_params) override
	{
		if (InputParams.size() != 2)
			throw std::runtime_error("Line2DModel - Number of input parameters does not match minimum number required for this model.");

		// Check for AbstractParamter types
		auto Point1 = std::dynamic_pointer_cast<Point2D>(InputParams[0]);
		auto Point2 = std::dynamic_pointer_cast<Point2D>(InputParams[1]);
		if (Point1 == nullptr || Point2 == nullptr)
			throw std::runtime_error("Line2DModel - InputParams type mismatch. It is not a Point2D.");

		std::copy(InputParams.begin(), InputParams.end(), m_MinModelParams.begin());

		// Compute the line parameters
		m_m = (Point2->m_Point2D[1] - Point1->m_Point2D[1]) / (Point2->m_Point2D[0] - Point1->m_Point2D[0]); // Slope
		m_d = Point1->m_Point2D[1] - m_m * Point1->m_Point2D[0]; // Intercept
		// m_d = Point2->m_Point2D[1] - m_m * Point2->m_Point2D[0]; // Intercept - alternative should be the same as above

		// mx - y + d = 0
		m_a = m_m;
		m_b = -1.0;
		m_c = m_d;

		m_DistDenominator = sqrt(m_a * m_a + m_b * m_b); // Cache square root for efficiency
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
};

}