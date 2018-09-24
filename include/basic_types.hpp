#ifndef BASIC_TYPES_HPP
#define BASIC_TYPES_HPP

#include "AbstractModel.hpp"

namespace GRANSAC
{

typedef std::array<VPFloat, 2> Vector2VP;

class Point2D: public AbstractParameter
{
public:
	Point2D(VPFloat x, VPFloat y)
	{
		m_Point2D[0] = x;
		m_Point2D[1] = y;
	};

	GRANSAC::Vector2VP m_Point2D;
};

}
#endif /* BASIC_TYPES_HPP */
