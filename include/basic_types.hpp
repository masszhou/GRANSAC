#ifndef BASIC_TYPES_HPP
#define BASIC_TYPES_HPP

#include <iostream>
#include <stdexcept>
#include <vector>
#include <memory>

namespace GRANSAC
{

// Each abstract model is made of abstract parameters
// Could be anything from a point (that make a 2D line or 3D plane or image correspondences) to a line
class AbstractParameter
{
public:
virtual ~AbstractParameter(void) {}; // To make this polymorphic we add dummy destructor
};

typedef std::array<float, 2> Vector2f;

class Point2D: public AbstractParameter
{
public:
	Point2D(float x, float y)
	{
		m_point2D[0] = x;
		m_point2D[1] = y;
	};

	GRANSAC::Vector2f m_point2D;
};

// Abstract model type for generic RANSAC model fitting
template <int t_num_params> /* Minimum number of parameters required to define this model*/
class AbstractModel
{
protected:

    std::vector<float> m_model_params; // for image processing, points are always 2D
    std::vector<std::shared_ptr<AbstractParameter> > m_model_def_parameters; // e.g. 2 points for line model, 3 points for quadratic line model. point could be line or plane two, since abstract type

    virtual float computeDistanceMeasure(std::shared_ptr<AbstractParameter> input_data) = 0;

public:
    virtual void initialize(const std::vector<std::shared_ptr<AbstractParameter>> &input_data, 
                            const std::map<std::string, float>& additional_params) = 0;
                            
    virtual std::pair<float, std::vector<std::shared_ptr<AbstractParameter>>> evaluate(const std::vector<std::shared_ptr<AbstractParameter>> &evaluate_params, float threshold) = 0;

    virtual std::vector<std::shared_ptr<AbstractParameter> > getModelDefParams(void) { return m_model_def_parameters; };
    virtual std::vector<float> getModelCoefficients(void) { return m_model_params; };
};


}
#endif /* BASIC_TYPES_HPP */
