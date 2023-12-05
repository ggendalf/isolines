#ifndef _ISOLINES_H_
#define _ISOLINES_H_

#pragma once

#include <vector>
#include <memory>
#include <string>
#include <map>
#include <shared_mutex>

#include "matrix.hpp"
#include "simple.h"

namespace geo::algorithm {

class isolines
{
public:
	using value_type						= geo::real_type;
	using value_array						= typename std::vector<value_type>;
	using value_matrix						= typename math::matrix2d<value_type>;
	using edge_matrix						= typename math::matrix2d<signed char>;
	using coord_array						= typename std::vector<coord_type>;

	using lines_type						= typename geo::line2d;
	using lines_array						= typename std::vector<std::unique_ptr<lines_type>>;
	using lines_result						= typename std::pair<lines_array, std::string>;

	using rings_type						= typename geo::ring2d;
	using rings_array						= typename std::vector<rings_type>;
	using rings_result						= typename std::pair<rings_array, std::string>;

	using polygones_type					= typename geo::polygon2d;
	using polygones_array					= typename std::vector<std::unique_ptr<polygones_type>>;
	using polygones_result					= typename std::pair<polygones_array, std::string>;

	using result_item						= typename std::tuple<lines_result, rings_result, polygones_result>;
	using result_map						= typename std::map<value_type, result_item>;

	isolines() noexcept;
	~isolines() = default;

	inline std::shared_mutex& get_mutex() { return m_guard; }

	bool empty() const noexcept;
	value_type get_null_value() const noexcept;
	box2d get_bound() const noexcept;
	value_matrix const& get_values() const noexcept;
	edge_matrix const& get_edges_vertical() const noexcept;
	edge_matrix const& get_edges_horizontal() const noexcept;
	coord_type get_eps_xy() const noexcept;
	value_type get_eps_value() const noexcept;
	coord_array const& get_array_xx() const noexcept;
	coord_array const& get_array_yy() const noexcept;
	value_type get_min_value() const noexcept;
	value_type get_max_value() const noexcept;

	void reset() noexcept;
	void init(value_matrix const& matrix_, box2d const& box_, value_type const& vblank_) noexcept(false);
	void set_holes(multipolygon2d const& holes_);
	value_type calc_value(point2d const& point_) const;

	void build(std::vector<value_type> const& levels_, result_map& result_, bool polygones_ = false) const noexcept(false);

protected:
	static void build_single(isolines const& source_, value_type const& level_,
							lines_array* lines_,
							rings_array* rings_ = nullptr) noexcept(false);

private:
	mutable std::shared_mutex				m_guard;

	value_type								m_null_value;
	coord_type								m_dx;
	coord_type								m_dy;
	box2d									m_bound;
	value_type								m_min_value;
	value_type								m_max_value;

	coord_type								m_eps_xy;
	value_type								m_eps_value;

	value_matrix							value_fun;   /* Массив данных [ny][nx]            */
	edge_matrix								gor_gran;    /* Массив гориз.  граней [ny][nx-1]  */
	edge_matrix								ver_gran;    /* Массив вертик. граней [ny-1][nx]  */
	coord_array								xx;          /* Массив координат сетки по X [nx]  */
	coord_array								yy;          /* Массив координат сетки по Y [ny]  */
	multipolygon2d							m_holes;

	void analyze_geometry() noexcept(false);
	void analyze_matrix() noexcept(false);
	void init_edges();
};

};

#endif//_ISOLINES_H_
