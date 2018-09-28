#ifndef BASIC_TYPES_HPP
#define BASIC_TYPES_HPP

#include "AbstractModel.hpp"

namespace GRANSAC
{

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

}
#endif /* BASIC_TYPES_HPP */
