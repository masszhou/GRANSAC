#ifndef LLS_MODEL_HPP
#define LLS_MODEL_HPP

#include "AbstractModel.hpp"
#include "basic_types.hpp"
#include <math.h>
#include <Eigen/QR>


namespace GRANSAC{

class QuadraticModel: public AbstractModel<2>
{
protected:
	// Parametric form
	VPFloat m_a, m_b, m_c; // ax + by + c = 0
	VPFloat m_DistDenominator; // = sqrt(a^2 + b^2). Stored for efficiency reasons
    
	//[1  x0  x0^2] [m_a0] = [y0]
    //[1  x1  x1^2] [m_a1]   [y1]
    //[1  x2  x2^2] [m_a2]   [y2]
    //...                    ...
    //[1  xn  xn^2]          [yn]
    VPFloat m_a0, m_a1, m_a2;

	// build a lookup table to calculate point to curve distance
	// e.g. target fitting is points in a 400x400 image
	// then grid of lookup table is 40x40 cell, each cell is 10x10 pixel
	int m_ngrid_x; // resolution, grid size oof lookup table
	int m_ngrid_y;
	std::vector<std::vector<int> > m_dist_lookup;

	// Another parametrization y = mx + d
	GRANSAC::VPFloat m_m; // Slope
	GRANSAC::VPFloat m_d; // Intercept

	virtual GRANSAC::VPFloat ComputeDistanceMeasure(std::shared_ptr<GRANSAC::AbstractParameter> Param) override
	{
		auto ExtPoint2D = std::dynamic_pointer_cast<Point2D>(Param);
		if (ExtPoint2D == nullptr)
			throw std::runtime_error("Line2DModel::ComputeDistanceMeasure() - Passed parameter are not of type Point2D.");

		// Return distance between passed "point" and this line
		// http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
		// GRANSAC::VPFloat Numer = fabs(m_a * ExtPoint2D->m_Point2D[0] + m_b * ExtPoint2D->m_Point2D[1] + m_c);
		// GRANSAC::VPFloat Dist = Numer / m_DistDenominator;

		// distance of point (x0, y0) to quadratic curve y(x)=a0+a1*x+a2*x*x
		// http://mathworld.wolfram.com/Point-QuadraticDistance.html
		// solve a cubic equation to find (x',y') on y(x)
		// 2*a2^2*x^3 + 3*a1a2*x^2 + (a1^2 + 2*a0a2 - 2*a2y0 + 1 )x + (a0a1 - a1y0 - x0)=0
		// correct but not efficient, I don't really need a high accuracy distance.
		// build a lookup table for distance calculation


		//// Debug
		//std::cout << "Point: " << ExtPoint2D->m_Point2D[0] << ", " << ExtPoint2D->m_Point2D[1] << std::endl;
		//std::cout << "Line: " << m_a << " x + " << m_b << " y + "  << m_c << std::endl;
		//std::cout << "Distance: " << Dist << std::endl << std::endl;

		// return Dist;
		return 0;
	};

public:
	QuadraticModel(const std::vector<std::shared_ptr<GRANSAC::AbstractParameter>> &InputParams)
	{
		Initialize(InputParams);
	};

	virtual void Initialize(const std::vector<std::shared_ptr<GRANSAC::AbstractParameter>> &InputParams) override
	{
		if (InputParams.size() != 2)
			throw std::runtime_error("QuadraticModel - Number of input parameters does not match minimum number required for this model.");

		// Check for AbstractParamter types
		auto Point1 = std::dynamic_pointer_cast<Point2D>(InputParams[0]);
		auto Point2 = std::dynamic_pointer_cast<Point2D>(InputParams[1]);
		if (Point1 == nullptr || Point2 == nullptr)
			throw std::runtime_error("QuadraticModel - InputParams type mismatch. It is not a Point2D.");


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

		// compute curve parameters
		std::vector<double> x_values, y_values, coeff;
		for (int i=0; i< InputParams.size(); i++){
			auto point = std::dynamic_pointer_cast<Point2D>(InputParams[i]);
			x_values.push_back(point->m_Point2D[0]);
			y_values.push_back(point->m_Point2D[1]);
		}
		polyfit(x_values, y_values, coeff, 2);
		m_a0 = coeff[0];
		m_a1 = coeff[1];
		m_a2 = coeff[2];

		//size=400x400, grid=40x40,
		//x=[5, 15, 20, 25, 30, 45, ..., 355, 365, 375, 385, 395]
		//y=[5, 15, 20, 25, 30, 45, ..., 355, 365, 375, 385, 395]
		//if curve go through cell, fill cell = 1, otherwise = 0

	};

	virtual std::pair<GRANSAC::VPFloat, std::vector<std::shared_ptr<GRANSAC::AbstractParameter>>> Evaluate(const std::vector<std::shared_ptr<GRANSAC::AbstractParameter>>& EvaluateParams, GRANSAC::VPFloat Threshold)
	{
		std::vector<std::shared_ptr<GRANSAC::AbstractParameter>> Inliers;
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

		GRANSAC::VPFloat InlierFraction = GRANSAC::VPFloat(nInliers) / GRANSAC::VPFloat(nTotalParams); // This is the inlier fraction

		return std::make_pair(InlierFraction, Inliers);
	};

private:
	void polyfit(const std::vector<double> &xv, const std::vector<double> &yv, std::vector<double> &coeff, int order)
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

#endif /* LLS_MODEL_HPP */
