#ifndef _SIMPLE_H_
#define _SIMPLE_H_

#pragma once

#include <limits>

#include <boost/geometry/core/cs.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/ring.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/multi_point.hpp>
#include <boost/geometry/geometries/multi_linestring.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>

#include <boost/geometry/algorithms/is_valid.hpp>
#include <boost/geometry/algorithms/equals.hpp>
#include <boost/geometry/algorithms/azimuth.hpp>

#include <boost/json.hpp>

#include <QPointF>

namespace geo {

using coord_type							= double;
using real_type								= double;

using boost_point2d							= typename boost::geometry::model::point<coord_type, 2, boost::geometry::cs::cartesian>;

class point2d : public boost_point2d
{
public:
	friend inline bool operator==(point2d const& lhs_, point2d const& rhs_)
	{
		if (rhs_.is_null() && lhs_.is_null())
			return true;

		return boost::geometry::equals<boost_point2d, boost_point2d>(lhs_, rhs_);
	}

	friend inline bool operator!=(point2d const& lhs_, point2d const& rhs_)
	{
		return !(rhs_ == lhs_);
	}

	friend inline boost::json::value const& operator>>(boost::json::value const& v_, point2d& p_)
	{
		assert(v_.kind() == boost::json::kind::array);
		if (boost::json::array const* arr = v_.if_array(); nullptr != arr && 2 == arr->size())
		{
			p_.set_x(static_cast<coord_type>((*arr)[0].as_double()));
			p_.set_y(static_cast<coord_type>((*arr)[1].as_double()));
		}

		return v_;
	}

	friend inline boost::json::value& operator<<(boost::json::value& v_, point2d const& p_)
	{
		boost::json::array& arr = v_.is_array() ? v_.as_array() : v_.emplace_array();

		if(2 != arr.size())
		{
			arr.resize(2);
		}

		arr[0] = p_.x();
		arr[1] = p_.y();

		return v_;
	}

	constexpr point2d() noexcept : boost_point2d(
				std::numeric_limits<coord_type>::quiet_NaN(),
				std::numeric_limits<coord_type>::quiet_NaN()) {}
	template<typename ctype>
	constexpr point2d(ctype const& x_, ctype const& y_) noexcept : boost_point2d(x_, y_) {};

	constexpr point2d(point2d const&) = default;
	constexpr point2d& operator=(point2d const&) = default;

	constexpr point2d(point2d&&) noexcept = default;
	constexpr point2d& operator=(point2d&&) noexcept = default;

	inline coord_type x() const noexcept { return get<0>(); }
	inline coord_type y() const noexcept { return get<1>(); }

	inline void set_x(coord_type const& x_) { boost_point2d::set<0>(x_); }
	inline void set_y(coord_type const& y_) { boost_point2d::set<1>(y_); }

	inline void set_null() noexcept {
		set_x(std::numeric_limits<coord_type>::quiet_NaN());
		set_y(std::numeric_limits<coord_type>::quiet_NaN()); }

	inline bool is_null() const noexcept { return std::isnan(x()) || std::isnan(y()); }

	inline real_type angle(point2d const& p_) const	{
		return !is_null() && !p_.is_null() ? boost::geometry::azimuth<boost_point2d, boost_point2d>(*this, p_) : .0; }

	constexpr point2d(QT_PREPEND_NAMESPACE(QPointF) const& qpt_) noexcept :
		boost_point2d(static_cast<coord_type>(qpt_.x()), static_cast<coord_type>(qpt_.y())) {};

	constexpr inline operator QT_PREPEND_NAMESPACE(QPointF)() const
	{ return QT_PREPEND_NAMESPACE(QPointF)(static_cast<qreal>(x()), static_cast<qreal>(y())); }

};

};//geo

BOOST_GEOMETRY_REGISTER_POINT_2D_GET_SET(geo::point2d, geo::coord_type, boost::geometry::cs::cartesian, x, y, set_x, set_y)

namespace geo {

using boost_box2d							= typename boost::geometry::model::box<point2d>;

class box2d : public boost_box2d
{
public:
	friend inline boost::json::value const& operator>>(boost::json::value const& v_, box2d& box_)
	{
		assert(v_.kind() == boost::json::kind::array);
		if(boost::json::array const* arr = v_.if_array(); nullptr != arr && 2 == arr->size())
		{
			(*arr)[0] >> box_.min_corner();
			(*arr)[1] >> box_.max_corner();
		}

		return v_;
	}

	friend inline boost::json::value& operator<<(boost::json::value& v_, box2d const& box_)
	{
		boost::json::array& arr = v_.is_array() ? v_.as_array() : v_.emplace_array();
		arr.resize(2);

		arr[0] << box_.min_corner();
		arr[1] << box_.max_corner();

		return v_;
	}

	constexpr box2d() noexcept : boost_box2d(point2d(), point2d()) {}
	template<typename ctype>
	constexpr explicit box2d(ctype const& x1_, ctype const& y1_, ctype const& x2_, ctype const& y2_) noexcept :
		boost_box2d(point2d(std::min<ctype>(x1_, x2_), std::min<ctype>(y1_, y2_)),
			point2d(std::max<ctype>(x1_, x2_), std::max<ctype>(y1_, y2_))) { }

	box2d(box2d const&) = default;
	box2d& operator=(box2d const&) = default;

	box2d(box2d&&) noexcept = default;
	box2d& operator=(box2d&&) noexcept = default;

	inline bool is_null() const noexcept { return min_corner().is_null() || max_corner().is_null(); }
	inline bool empty() const noexcept { return is_null() || 0. == width() || 0. == height(); }
	inline bool is_valid() const { return !is_null() && boost::geometry::is_valid<boost_box2d>(*this); }

	inline coord_type width() const noexcept { return abs(max_x() - min_x()); }
	inline coord_type height() const noexcept { return abs(max_y() - min_y()); }

	inline coord_type min_x() const noexcept { return min_corner().x(); }
	inline coord_type min_y() const noexcept { return min_corner().y(); }
	inline coord_type max_x() const noexcept { return max_corner().x(); }
	inline coord_type max_y() const noexcept { return max_corner().y(); }

	void set_null() noexcept { min_corner().set_null(); max_corner().set_null();	}

};

using line2d								= typename boost::geometry::model::linestring<point2d>;
using ring2d								= typename boost::geometry::model::ring<point2d>;
using polygon2d								= typename boost::geometry::model::polygon<point2d>;
using multipoint2d							= typename boost::geometry::model::multi_point<point2d>;
using multiline2d							= typename boost::geometry::model::multi_linestring<line2d>;
using multipolygon2d						= typename boost::geometry::model::multi_polygon<polygon2d>;

template<typename container>
boost::json::value& collection2json(boost::json::value& v_, container const& c_)
{
	if(!c_.empty())
	{
		boost::json::array& arr = v_.is_array() ? v_.as_array() : v_.emplace_array();
		arr.resize(c_.size());

		boost::json::array::size_type loop = 0;
		for(auto const& item : c_)
		{
			arr[loop] << item;
			++loop;
		}
	}

	return v_;
}

inline boost::json::value& operator<<(boost::json::value& v_, line2d const& l_)
{
	return collection2json<line2d>(v_, l_);
}

inline boost::json::value& operator<<(boost::json::value& v_, ring2d const& r_)
{
	return collection2json<ring2d>(v_, r_);
}

inline boost::json::value& operator<<(boost::json::value& v_, polygon2d const& pp_)
{
	boost::json::array& arr = v_.is_array() ? v_.as_array() : v_.emplace_array();

	arr.resize(pp_.inners().size() + 1);
	collection2json<polygon2d::ring_type>(arr[0], pp_.outer());
	if (!pp_.inners().empty())
	{
		std::size_t loop = 1;
		for(auto const& ring : pp_.inners() )
		{
			collection2json(arr[loop], ring);
			++loop;
		}
	}

	return v_;
}

};

#endif //_SIMPLE_H_
