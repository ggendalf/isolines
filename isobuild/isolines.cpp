#include <memory>
#include <algorithm>
#include <iterator>
#include <array>

#include "isolines.h"

#include <boost/geometry/algorithms/correct.hpp>
#include <boost/geometry/algorithms/simplify.hpp>
#include <boost/geometry/algorithms/is_valid.hpp>
#include <boost/geometry/algorithms/covered_by.hpp>
#include <boost/geometry/algorithms/area.hpp>

namespace geo::algorithm {

using vtype												= typename isolines::value_type;

static constexpr vtype min_dxy							= 1.;//метр

constexpr static vtype default_epsilon					= 1.e-12;
constexpr static coord_type default_epsilon_xy			= 1.e-6;
constexpr static vtype default_multiplicator			= (1./default_epsilon);
constexpr static std::size_t point_buffer_size			= 128;
constexpr static std::size_t isoline_size				= 512;
constexpr static coord_type minimal_distance			= 0.00999;//for simplify
constexpr static vtype default_blank					= std::numeric_limits<vtype>::min();

//////////////////////////////////////////////////////////////////////
/// \brief detABC
/// \param p1
/// \param p2
/// \param p3
/// \return false - если направление по часовой стрелки, true -- если против.
///
inline bool detABC( point2d const& p1, point2d const& p2, point2d const& p3 )
{
	return 0. < ((p2.x() - p1.x())*(p3.y() - p1.y()) - (p2.y() - p1.y())*(p3.x() - p1.x())) ? true : false;
}

#ifdef _DEBUG
#define angle_shift  + boost::math::constants::pi<vtype>()
#else
#define angle_shift
#endif//_DEBUG

//////////////////////////////////////////////////////////////////////
/// \brief The iso_params struct
/// параметры изолинии в клетке
struct iso_params
{
	/* xx=x-sx;  yy=y-sy;  vg=xx*yy. */
	coord_type								shift_x;/* координаты центра изолинии     */
	coord_type								shift_y;
	vtype									val_giper;          /* значение гиперболы в КСК       */
	vtype									val_deter;          /* детерминант квадратичной формы */
	vtype									sign_x;/* знаки в КСК == ветвь гиперболы */
	vtype									sign_y;

	iso_params() noexcept { std::memset(this, 0, sizeof(*this)); }
};

//////////////////////////////////////////////////////////////////////
/// \brief The spline_param struct
/// параметры сплайна в клетке
struct spline_param
{
	/* xx=x-dx;   yy=y-dy;   z=a0+ax*xx+ay*yy+axy*xx*yy. */
	coord_type								dx;			/* координаты узла <ix,iy> в КСК */
	coord_type								dy;
	vtype									k_a0;        /* коэффициент сплайна - const   */
	vtype									k_ax;        /* коэффициент сплайна при "X"   */
	vtype									k_ay;        /* коэффициент сплайна при "Y"   */
	vtype									k_axy;       /* коэффициент сплайна при "X*Y" */

	spline_param() noexcept { std::memset(this, 0, sizeof(*this)); }

	vtype calc_value(point2d const& target_) const noexcept
	{
		vtype cx = static_cast<vtype>(target_.x() - dx);
		vtype cy = static_cast<vtype>(target_.y() - dy);

		return k_a0 + k_ax * cx + k_ay * cy + k_axy * cx * cy;
	}

	bool calc_iso_point(line2d& isoline_, iso_params const& giper_, vtype const& eps_xy_, std::size_t const index_)
	{
		vtype w0, w1, ww, distance, _sq;

		std::size_t i_curr = index_;
		std::size_t i_next = i_curr + 1;

		coord_type x0 = isoline_[i_curr].x();
		coord_type y0 = isoline_[i_curr].y();
		coord_type x1 = isoline_[i_next].x();
		coord_type y1 = isoline_[i_next].y();

		if (fabs(x1 - x0) >= fabs(y1 - y0))
		{
			w0 = fabs(k_ay + static_cast<vtype>(x0) * k_axy);
			w1 = fabs(k_ay + static_cast<vtype>(x1) * k_axy);
			ww = giper_.val_deter;


			_sq = boost::geometry::math::sqr<vtype>(ww) +
				boost::geometry::math::sqr<vtype>(w0) *
				boost::geometry::math::sqr<vtype>(w1);

			distance = std::fabs((k_axy * ww * static_cast<vtype>((x1 - x0) * (x1 - x0))) /
				((w0 + w1 + std::sqrt(w0 * w1)) * std::sqrt(_sq)));

			if (distance < eps_xy_)
			{
				return false;
			}

			isoline_.insert((isoline_.begin() + i_next), isoline_[i_next]);

			_sq = static_cast<vtype>((x1 - giper_.shift_x)*(x0 - giper_.shift_x));
			if (_sq < 0.)
			{
				_sq = 0.;
			}

			ww = giper_.sign_x * std::sqrt(_sq);
			isoline_[i_next].set_x(giper_.shift_x + static_cast<coord_type>(ww));
			isoline_[i_next].set_y(giper_.shift_y + static_cast<coord_type>(giper_.val_giper / ww));
		}
		else
		{
			w0 = std::fabs(k_ax + static_cast<vtype>(y0) * k_axy);
			w1 = std::fabs(k_ax + static_cast<vtype>(y1) * k_axy);
			ww = giper_.val_deter;

			distance = std::fabs((k_axy * ww * static_cast<vtype>((y1 - y0) * (y1 - y0))) /
				((w0 + w1 + std::sqrt(w0 * w1)) *
					std::sqrt(
						boost::geometry::math::sqr<vtype>(ww) +
						boost::geometry::math::sqr<vtype>(w0) *
						boost::geometry::math::sqr<vtype>(w1))));

			if (distance < eps_xy_)
			{
				return false;
			}

			isoline_.insert((isoline_.begin() + i_next), isoline_[i_next]);

			_sq = static_cast<vtype>((y1 - giper_.shift_y) * (y0 - giper_.shift_y));
			if (_sq < 0.)
			{
				_sq = 0.;
			}

			ww = giper_.sign_y * std::sqrt(_sq);
			isoline_[i_next].set_y(giper_.shift_y + static_cast<coord_type>(ww));
			isoline_[i_next].set_x(giper_.shift_x + static_cast<coord_type>(giper_.val_giper / ww));
		}

		return  true;
	}

};

//////////////////////////////////////////////////////////////////////
/// \brief The box_params class
/// параметры клетки для анализа изолиний
class box_params
{
public:
	box_params(isolines const* source_, vtype level_) :
		m_level(level_),
		m_source(source_)
	{
		m_iso_begin.set_null();
		m_iso_end.set_null();

		m_type_inp = 0;
		m_ix_gran = m_iy_gran = m_ix = m_iy = 0;
		m_hx = m_hy = m_xx = m_yy = 0.;

		m_buffer.reserve(point_buffer_size);

		assert(nullptr != m_source);
		assert(0.0 != m_level);

		m_edgesH = m_source->get_edges_horizontal();
		m_edgesV = m_source->get_edges_vertical();
		m_values = m_source->get_values();

		//Юстировка значений функции в узлах я целью избежать вырождения
		vtype eps = get_eps_value();
		for(vtype& v : m_values)
		{
			if (fabs(v - m_level) < eps)
			{
				v += ((0 <= v) ? eps : -eps);
			}
		}
	}

	void build(isolines::lines_array* lines_, isolines::rings_array* rings_ = nullptr)
	{
		if( nullptr == lines_ && nullptr == rings_)
		{
			throw std::exception("пустой запрос построения");
		}

		if(m_level < get_min_value() || m_level > get_max_value())
		{
			throw std::exception("уровень изолинии лежит вне ОДЗ");
		}

		isolines::lines_array res;
		isolines::lines_array* lines = nullptr != lines_ ? lines_ : std::addressof(res);

		int32_t ii;
		int32_t start_type, logic_gran, logic_zamk;
		point2d start_izol;

		std::size_t nx = m_values.cols();
		std::size_t ny = m_values.rows();

		isolines::lines_type buffer;
		int32_t w_closed = 0;

		buffer.reserve(isoline_size);

		for (std::size_t iy = 0u; iy < ny; ++iy)
		{
			for (std::size_t ix = 0u; ix < nx - 1u; ++ix)
			{
				if ((ii = m_edgesH.at(iy, ix)) <= 0)//?
				{
					continue;
				}

				if ((m_values.at(iy, ix) - m_level) * (m_values.at(iy, ix + 1u) - m_level) >= 0.)
				{
					m_edgesH.set(iy, ix, static_cast<isolines::edge_matrix::value_type>(-ii));
					continue;
				}

				buffer.clear();
				logic_gran = ii;
				start_type = (ii == 2) ? 3 : 1;
				w_closed = 0;
				init_first(start_type, ix, iy);

				start_izol = m_iso_begin;

				buffer.push_back(start_izol);

				if (logic_gran < 3)
				{
					do
					{
						ii = do_isoline(buffer);
						if (0 == ii)
						{
							break;
						}
						else if(0 > ii)
						{
							throw std::exception("ошибка построения изолиний");
						}
					}
					while (true);
				}
				else
				{  // logic_gran == 3
					logic_zamk = 0;
					do
					{
						ii = do_isoline(buffer);
						if (0 == ii)
						{
							break;
						}
						else if(0 > ii)
						{
							throw std::exception("ошибка построения изолиний");
						}

						if (ix == m_ix  &&  iy == m_iy + 1)
						{
							m_edgesH.set(iy, ix, static_cast<isolines::edge_matrix::value_type>(logic_gran));
						}

						if (start_type == m_type_inp && ix == m_ix && iy == m_iy)
						{
							m_edgesH.set(iy, ix, static_cast<isolines::edge_matrix::value_type>(-logic_gran));
							logic_zamk = 1;

							break;
						}
					}
					while (true);

					if (logic_zamk != 1)
					{
						start_type = 3;
						m_edgesH.set(iy, ix, static_cast<isolines::edge_matrix::value_type>(logic_gran));
						init_params(ix, iy, start_izol);

						if (-1 == w_closed)
						{
							if(finish_isoline(buffer, w_closed))
							{
								//lines->push_back(buffer);
								lines->emplace_back(std::make_unique<isolines::lines_type>(buffer));
							}

							buffer.clear();
							buffer.push_back(start_izol);
						}
						else
						{
							std::reverse(buffer.begin(), buffer.end());
						}

						do
						{
							ii = do_isoline(buffer);

							if (0 == ii)
							{
								break;
							}
							else if (0 > ii)
							{
								throw std::exception("ошибка построения изолиний");
							}
						}
						while (true);
					}
				}

				if(finish_isoline(buffer, w_closed))
				{
					//lines->push_back(buffer);
					lines->emplace_back(std::make_unique<isolines::lines_type>(buffer));
				}
			}
		}// end of ix // end of iy

		for (std::size_t iy = 0u; iy < ny - 1u; ++iy)
		{
			for (std::size_t ix = 0u; ix < nx; ++ix)
			{
				if ((ii = m_edgesV.at(iy, ix)) <= 0)
				{
					continue;
				}

				if ((m_values.at(iy, ix) - m_level) * (m_values.at(iy + 1u, ix) - m_level) >= 0.)
				{
					m_edgesV.set(iy, ix, static_cast<isolines::edge_matrix::value_type>(-ii));
					continue;
				}

				if (ii != 1)
				{
					throw std::exception("ошибка построения изолиний");
				}

				buffer.clear();
				start_type = 4;
				w_closed = 0;
				init_first(start_type, ix, iy);
				start_izol = m_iso_begin;

				buffer.push_back(start_izol);

				do
				{
					ii = do_isoline(buffer);
					if (0 == ii)
					{
						break;
					}
					else if (0 > ii)
					{
						throw std::exception("ошибка построения изолиний");
					}
				}
				while (true);

				if(finish_isoline(buffer, w_closed))
				{
					//lines->push_back(buffer);
					lines->emplace_back(std::make_unique<isolines::lines_type>(buffer));
				}
			}
		}

		if(rings_)
		{
			create_contours(*rings_, *lines, !lines_);
		}
	}

private:
	point2d									m_iso_begin; /* коор. начальной точки изолинии в клетке */
	point2d									m_iso_end; /* коор. конечной точки изолинии в клетке  */

	int32_t									m_type_inp;
	std::size_t								m_ix_gran;/* индексы текущей грани      */
	std::size_t								m_iy_gran;
	std::size_t								m_ix;/* индексы текущей клетки     */
	std::size_t								m_iy;

	coord_type								m_hx;/* размер клетки по осям "X", "Y"  */
	coord_type								m_hy;
	coord_type								m_xx;/* координаты узла  <"ix", "iy">   */
	coord_type								m_yy;

	vtype									m_level;/* значения уровня изолинии            */
	std::array<vtype, 4>					m_value_fun;// значения функции в узлах клетки
	line2d									m_buffer;

protected:
	isolines::value_matrix					m_values;
	isolines::edge_matrix					m_edgesH;
	isolines::edge_matrix					m_edgesV;

	isolines const*							m_source;

	inline coord_type get_eps_xy() const noexcept { return nullptr != m_source ? m_source->get_eps_xy() : default_epsilon_xy; }
	inline vtype get_eps_value() const noexcept { return nullptr != m_source ? m_source->get_eps_value() : default_epsilon; }
	inline box2d get_bound() const noexcept { return nullptr != m_source ? m_source->get_bound() : box2d(); }
	inline vtype get_null_value() const noexcept { return nullptr != m_source ? m_source->get_null_value() : default_blank; }
	inline vtype get_min_value() const noexcept { return nullptr != m_source ? m_source->get_min_value() : default_blank; }
	inline vtype get_max_value() const noexcept { return nullptr != m_source ? m_source->get_max_value() : default_blank; }

	inline vtype calc_value(point2d const& p_) const { return nullptr != m_source ? m_source->calc_value(p_) : get_null_value(); }
	inline coord_type get_xx(size_t idx_) const { return nullptr != m_source ? m_source->get_array_xx()[idx_] : coord_type(); }
	inline coord_type get_yy(size_t idx_) const { return nullptr != m_source ? m_source->get_array_yy()[idx_] : coord_type(); }

	void init()
	{
		m_ix = (m_type_inp == 2) ? (m_ix_gran - 1u) : m_ix_gran;
		m_iy = (m_type_inp == 3) ? (m_iy_gran - 1u) : m_iy_gran;

		m_xx = get_xx(m_ix);
		m_yy = get_yy(m_iy);
		m_hx = get_xx(m_ix + 1u) - get_xx(m_ix);
		m_hy = get_yy(m_iy + 1u) - get_yy(m_iy);

		m_value_fun[0] = m_values.at(m_iy, m_ix);
		m_value_fun[1] = m_values.at(m_iy, m_ix + 1u);
		m_value_fun[2] = m_values.at(m_iy + 1u, m_ix + 1u);
		m_value_fun[3] = m_values.at(m_iy + 1u, m_ix);
		m_iso_begin = m_iso_end;
	}

	void init_first(int32_t type_, std::size_t ix_, std::size_t iy_)
	{
		m_type_inp = type_;
		m_ix_gran = ix_;
		m_iy_gran = iy_;
		m_ix = ix_;          // type = 1 || 4 || 3
		m_iy = (3 == type_) ? (iy_ - 1u) : iy_;

		m_xx = get_xx(m_ix);
		m_yy = get_yy(m_iy);
		m_hx = get_xx(m_ix + 1u) - get_xx(m_ix);
		m_hy = get_yy(m_iy + 1u) - get_yy(m_iy);

		m_value_fun[0] = m_values.at(m_iy, m_ix);
		m_value_fun[1] = m_values.at(m_iy, m_ix + 1u);
		m_value_fun[2] = m_values.at(m_iy + 1u, m_ix + 1u);
		m_value_fun[3] = m_values.at(m_iy + 1u, m_ix);

		switch (type_)
		{
		case  1:
			if (m_value_fun[0] == m_value_fun[1])
			{
				m_iso_begin.set_x(m_xx);
			}
			else
			{
				m_iso_begin.set_x((m_value_fun[0] - m_level) /
					(m_value_fun[0] - m_value_fun[1]) * m_hx + m_xx);
			}

			m_iso_begin.set_y(m_yy);
			break;

		case  3:
			if (m_value_fun[2] == m_value_fun[3])
			{
				m_iso_begin.set_x(m_xx);
			}
			else
			{
				m_iso_begin.set_x((m_value_fun[3] - m_level) /
					(m_value_fun[3] - m_value_fun[2]) * m_hx + m_xx);
			}

			m_iso_begin.set_y(m_yy + m_hy);
			break;

		case  4:
			if (m_value_fun[0] == m_value_fun[3])
			{
				m_iso_begin.set_y(m_yy);
			}
			else
			{
				m_iso_begin.set_y((m_value_fun[0] - m_level) /
					(m_value_fun[0] - m_value_fun[3]) * m_hy + m_yy);
			}

			m_iso_begin.set_x(m_xx);
			break;
		}
	}

	void init_params(std::size_t ix_gran_, std::size_t iy_gran_, point2d const& iso_end_)
	{
		m_type_inp = 3;
		m_ix_gran = ix_gran_;
		m_iy_gran = iy_gran_;

		m_ix = m_ix_gran;
		m_iy = m_iy_gran - 1u;

		m_xx = get_xx(m_ix);
		m_yy = get_yy(m_iy);
		m_hx = get_xx(m_ix + 1u) - get_xx(m_ix);
		m_hy = get_yy(m_iy + 1u) - get_yy(m_iy);

		m_value_fun[0] = m_values.at(m_iy, m_ix);
		m_value_fun[1] = m_values.at(m_iy, m_ix + 1u);
		m_value_fun[2] = m_values.at(m_iy + 1u, m_ix + 1u);
		m_value_fun[3] = m_values.at(m_iy + 1u, m_ix);
		m_iso_begin = iso_end_;
	}

	void restore_edges()
	{
		for(isolines::edge_matrix::value_type& v : m_edgesH)
		{
			v = -v;
		};

		for(isolines::edge_matrix::value_type& v : m_edgesV)
		{
			v = -v;
		};
	}

	void init_param_giper(iso_params& giper_, spline_param& bilin_) const
	{
		bilin_.k_a0 = m_level - m_value_fun[0];
		bilin_.k_ax = (m_value_fun[1] - m_value_fun[0]) / static_cast<isolines::value_type>(m_hx);
		bilin_.k_ay = (m_value_fun[3] - m_value_fun[0]) / static_cast<isolines::value_type>(m_hy);
		bilin_.k_axy = (m_value_fun[0] - m_value_fun[1] + m_value_fun[2] - m_value_fun[3]) / (m_hx * m_hy);

		if (0.0 != bilin_.k_axy)
		{
			giper_.shift_x = static_cast<coord_type>(-bilin_.k_ay / bilin_.k_axy);
			giper_.shift_y = static_cast<coord_type>(-bilin_.k_ax / bilin_.k_axy);
			giper_.val_deter = bilin_.k_a0 * bilin_.k_axy + bilin_.k_ax * bilin_.k_ay;
			giper_.val_giper = giper_.val_deter / (bilin_.k_axy * bilin_.k_axy);
		}

		giper_.sign_x = ((m_iso_begin.x() - m_xx) > giper_.shift_x) ? 1. : -1.;
		giper_.sign_y = ((m_iso_begin.y() - m_yy) > giper_.shift_y) ? 1. : -1.;
	}

	template<class operation_T>
	void shift_points(line2d& isoline_) const
	{
		for(point2d& v : isoline_)
		{
			v.set_x(operation_T{}(v.x(), m_xx));
			v.set_y(operation_T{}(v.y(), m_yy));
		}
	}

	void create_contours(isolines::rings_array& contours_, isolines::lines_array& lines_, bool lines_int_)
	{
		class line2ring_inserter : public std::back_insert_iterator<isolines::rings_array>
		{
		public:
			explicit line2ring_inserter(container_type& c_) :
				std::back_insert_iterator<isolines::rings_array>(c_)
			{
			}

			line2ring_inserter& operator=(const typename isolines::lines_array::value_type& line_)
			{
				//container->emplace_back(std::make_unique<ring2d>(line_->begin(), line_->end()));
				container->emplace_back(container_type::value_type(line_->begin(), line_->end()));
				return (*this);
			}

			_NODISCARD line2ring_inserter& operator*()
			{
				return (*this);
			}
		};

		contours_.reserve(lines_.size());
		std::copy_if(lines_.begin(), lines_.end(),
					 line2ring_inserter(contours_),
		[](isolines::lines_array::value_type const& l)
		{
			return l->front() == l->back();
		});

		isolines::lines_array lines_buffer;
		isolines::lines_array* plines = nullptr;

		if(lines_int_)//внутренний буфер
		{
			lines_.erase(std::remove_if(lines_.begin(), lines_.end(),
			[](isolines::lines_array::value_type const& l)
			{
				return l->front() == l->back();
			}), lines_.end());

			plines = std::addressof(lines_);
		}
		else
		{
			for(isolines::lines_array::value_type const& l : lines_)
			{
				if(l->front() != l->back())
				{
					//lines_buffer.push_back(l);
					lines_buffer.emplace_back(std::make_unique<isolines::lines_type>(*l));
				}
			}

			plines = std::addressof(lines_buffer);
		}
		assert(nullptr != plines);

		//остались только линии
		if(!plines->empty())
		{
			box2d bounds = get_bound();
			point2d center(boost::geometry::return_centroid<point2d, boost_box2d>(bounds));
			ring2d area;

			using cycle_item = typename std::tuple<vtype, point2d, line2d const*>;
			using cycle_list = typename std::list<cycle_item>;
			cycle_list cycle_points;
			for(isolines::lines_array::value_type& l : (*plines))
			{
				assert(l->front() != l->back());
				if(detABC(l->front(), l->back(), center))
				{
					std::reverse(l->begin(), l->end());
				}

				cycle_points.emplace_back(std::make_tuple(center.angle(l->front())angle_shift, l->front(), l.get()));
				cycle_points.emplace_back(std::make_tuple(center.angle(l->back())angle_shift, l->back(), nullptr));
			}

			cycle_points.sort([](cycle_item const& lh, cycle_item const& rh)
			{
				return std::get<0>(lh) < std::get<0>(rh);
			});

			//line2d border_points(get_bound());
			line2d border_points;
			border_points.push_back(bounds.min_corner());
			border_points.push_back(point2d(bounds.min_corner().x(), bounds.max_corner().y()));
			border_points.push_back(bounds.max_corner());
			border_points.push_back(point2d(bounds.max_corner().x(), bounds.min_corner().y()));

			vtype angle_start = 0.;
			vtype angle_end = 0.;
			bool border_insert = false;
			for(cycle_item const& ip : cycle_points)
			{
				auto [angle, point, line] = ip;
				if(border_insert = nullptr == line; !border_insert)
				{
					if(!area.empty())
					{
						angle_start = center.angle(area.back())angle_shift;
					}
					angle_end = angle;

					std::remove_copy_if(border_points.begin(), border_points.end(),
										std::back_inserter(area), [&center, &angle_start, &angle_end](point2d const& p)
					{
						vtype a = center.angle(p)angle_shift;
						return !(a > angle_start && a < angle_end);
					});

					area.insert(area.end(), line->begin(), line->end());
				}
				else
				{
					angle_start = angle;
				}
			}

			if(border_insert)
			{
				std::remove_copy_if(border_points.begin(), border_points.end(),
								std::back_inserter(area), [&center, &angle_start](point2d const& p)
				{
					vtype a = center.angle(p)angle_shift;
					return !(a > angle_start && a < boost::math::constants::two_pi<vtype>());
				});
			}

			boost::geometry::correct<ring2d>(area);
			assert(boost::geometry::is_valid<ring2d>(area));

			restore_edges();

			//проверка на инверсию
			static constexpr coord_type step_xy = 0.1;
			bool need_inverse = false;
			ring2d::const_iterator ip = area.begin() + area.size() / 2;
			bool need_check = true;
			do
			{
				point2d test[4];
				size_t idx = 0;
				for(point2d& p : test)
				{
					p.set_x((*ip).x() + ( 0 == idx % 2 ? step_xy : -step_xy));
					p.set_y((*ip).y() + ( 0 == idx % 3 ? step_xy : -step_xy));

					//covered_by
					if( boost::geometry::within<point2d, ring2d>(p, area))
					{
						vtype cv = calc_value(p);
						if(get_null_value() != cv)
						{
							need_check = false;
							need_inverse = cv < m_level;
							break;
						}
					}

					++idx;
				}

				++ip;
			}while(need_check && area.end() != ip);

			if(need_inverse)
			{
				ring2d bound_area;
				boost::geometry::convert<boost_box2d, ring2d>(bounds, bound_area);

				std::vector<ring2d> res;
				boost::geometry::difference<ring2d, ring2d, std::vector<ring2d>>(bound_area, area, res);

				contours_.insert(contours_.end(), res.begin(), res.end());
			}
			else
			{
				contours_.push_back(area);
			}
		}
	}

	//works_izol
	int32_t do_isoline(line2d& isoline_)
	{
		int32_t res = 0;
		if(!analyze())
			return -2;  //err -> т.к. НЕ замкнуто!

		res = do_isopoints();

		isoline_.insert(isoline_.end(), m_buffer.begin() + 1u, m_buffer.end());

		if (!can_isoline())
			return 0;

		init();

		return  res;
	}

	std::size_t do_isopoints()
	{
		iso_params giper;
		spline_param bilin;

		m_buffer.clear();
		m_buffer.push_back(m_iso_begin);
		m_buffer.push_back(m_iso_end);

		vtype d_fun = m_value_fun[0] - m_value_fun[1] + m_value_fun[2] - m_value_fun[3];
		if (fabs(d_fun) < get_eps_value())
			return m_buffer.size();

		shift_points<std::minus<coord_type>>(m_buffer);

		init_param_giper(giper, bilin);

		std::size_t index_works_pnt = 0u;
		do
		{
			if (!bilin.calc_iso_point(m_buffer, giper, get_eps_xy(), index_works_pnt))
			{
				index_works_pnt++;
			}
		}
		while (index_works_pnt < (m_buffer.size() - 1u) && m_buffer.size() < m_buffer.capacity());

		shift_points<std::plus<coord_type>>(m_buffer);

		return m_buffer.size();
	}

	bool finish_isoline(line2d& buffer_, int32_t& closed_) const
	{
		if (!buffer_.empty())
		{
			point2d beg(buffer_.front());
			point2d& end(buffer_.back());

			if (-1 != closed_)
			{
				closed_ = (fabs(beg.x() - end.x()) < get_eps_xy()  &&  fabs(beg.y() - end.y()) < get_eps_xy());
				if(closed_)
				{
					end = beg;
					//buffer_.front() = buffer_.back();
				}
			}

			line2d result;
			result.reserve(buffer_.size());

			boost::geometry::simplify<line2d, line2d>(buffer_, result, minimal_distance);
			if (1 == closed_ && result.front() != result.back())
			{
				throw std::exception("ошибка построения");
			}

			buffer_ = result;

			return true;
		}

		return false;
	}

	inline point2d get_min_xy() const noexcept
	{
		return point2d(
			std::min<coord_type>(m_iso_begin.x(), m_iso_end.x()),
			std::min<coord_type>(m_iso_begin.y(), m_iso_end.y())
		);
	}

	inline point2d get_max_xy() const noexcept
	{
		return point2d(
			std::max<coord_type>(m_iso_begin.x(), m_iso_end.x()),
			std::max<coord_type>(m_iso_begin.y(), m_iso_end.y())
		);
	}

	bool coordinate_control() const noexcept
	{
		coord_type dx = fabs(m_hx) + static_cast<coord_type>(get_eps_xy());
		coord_type dy = fabs(m_hy) + static_cast<coord_type>(get_eps_xy());

		return !(fabs(m_iso_end.x() - m_xx) > dx
			|| fabs(m_iso_end.x() - m_xx - m_hx) > dx
			|| fabs(m_iso_end.y() - m_yy) > dy
			|| fabs(m_iso_end.y() - m_yy - m_hy) > dy
			|| fabs(m_iso_begin.x() - m_xx) > dx
			|| fabs(m_iso_begin.x() - m_xx - m_hx) > dx
			|| fabs(m_iso_begin.y() - m_yy) > dy
			|| fabs(m_iso_begin.y() - m_yy - m_hy) > dy);
	}

	bool can_isoline()
	{
		//int32_t ii;
		isolines::edge_matrix::value_type ii;
		switch (m_type_inp)
		{
		case  1:
		case  3:
			ii = m_edgesH.at(m_iy_gran, m_ix_gran);
			switch (ii)
			{
			case  3:
				break;

			case  1:
			case  2:
				m_edgesH.set(m_iy_gran, m_ix_gran, -ii);

			default:
				return false;
			}
			break;

		case  2:
		case  4:
			ii = m_edgesV.at(m_iy_gran, m_ix_gran);
			switch (ii)
			{
			case  3:
				break;

			case  1:
			case  2:
				m_edgesV.set(m_iy_gran, m_ix_gran, -ii);

			default:
				return false;
			}
			break;
		}

		return true;
	}

	bool analyze()
	{

		coord_type xk = 0.0, yk = 0.0;
		coord_type x1, x2, x3;
		coord_type y1, y2, y3;

		int32_t type_val = m_type_inp;
		std::size_t ix_g = m_ix_gran;
		std::size_t iy_g = m_iy_gran;

		bool log1, log2, log3;
		std::size_t ii;
		switch (type_val)
		{
		case  1:
		case  3:
			m_edgesH.set(iy_g, ix_g, -abs(m_edgesH.at(iy_g, ix_g)));

			log1 = (m_edgesV.at(m_iy, m_ix) > 0) &&
				((m_value_fun[0] - m_level) * (m_value_fun[3] - m_level) < 0.);

			log2 = (m_edgesV.at(m_iy, m_ix + 1u) > 0) &&
				((m_value_fun[1] - m_level) * (m_value_fun[2] - m_level) < 0.);

			ii = (type_val != 3) ? 1u : 0u;

			log3 = (m_edgesH.at(m_iy + ii, m_ix) > 0) &&
				((m_value_fun[0u + 2u * ii] - m_level) * (m_value_fun[1u + 2u * ii] - m_level) < 0.);

			ii = static_cast<std::size_t>(log1) +
				static_cast<std::size_t>(log2) +
				static_cast<std::size_t>(log3);

			if (ii == 0u || ii == 2u)
			{
				//	fprintf( ErrFile, "\nError_gor:  ixk=%d  iyk=%d  nn=%d\n", ix_k, iy_k, ii );
				return false;
			}

			if (log1)
			{
				y1 = static_cast<coord_type>((m_value_fun[0] - m_level) / (m_value_fun[0] - m_value_fun[3]) * m_hy + m_yy);
			}

			if (log2)
			{
				y2 = static_cast<coord_type>((m_value_fun[1] - m_level) / (m_value_fun[1] - m_value_fun[2]) * m_hy + m_yy);
			}

			if (log1 == log2)
			{
				x3 = static_cast<coord_type>((
					(m_type_inp == 3) ?
					(m_value_fun[0] - m_level) / (m_value_fun[0] - m_value_fun[1]) * m_hx + m_xx :
					(m_value_fun[3] - m_level) / (m_value_fun[3] - m_value_fun[2]) * m_hx + m_xx));
			}

			if (log1)
			{
				if (log2)
				{
					if (fabs(y1 - y2) > get_eps_xy() || fabs(x3 - m_iso_begin.x()) > get_eps_xy())
					{
						if (m_type_inp == 1)
						{
							if ((y1 - y2 <= 0.) == (m_hy > 0.))
							{ // y1<=y2
								xk = m_xx;
								yk = y1;

								m_type_inp = 2;
							}
							else
							{
								xk = m_xx + m_hx;
								yk = y2;

								m_type_inp = 4;
								m_ix_gran++;
							}
						}
						else
						{
							if ((y1 - y2 > 0.) == (m_hy > 0.))
							{ // y1>y2
								xk = m_xx;
								yk = y1;

								m_type_inp = 2;
								m_iy_gran--;
							}
							else
							{
								xk = m_xx + m_hx;
								yk = y2;

								m_type_inp = 4;
								m_ix_gran++;
								m_iy_gran--;
							}
						}
					}
					else
					{
						xk = x3;
						if (m_type_inp == 1)
						{
							m_iy_gran++;

							yk = m_yy + m_hy;
						}
						else
						{
							m_iy_gran--;
							yk = m_yy;
						}
					}
				}   // The end log2 == 1;

				else
				{
					ii = static_cast<std::size_t>(m_type_inp == 1);

					m_edgesH.set(m_iy + ii, m_ix, -abs(m_edgesH.at(m_iy + ii, m_ix)));
					m_edgesV.set(m_iy, m_ix + 1u, -abs(m_edgesV.at(m_iy, m_ix + 1u)));

					xk = m_xx;
					yk = y1;

					if (m_type_inp == 3)
					{
						m_iy_gran--;
					}

					m_type_inp = 2;
				}
			}   // The end log1 == 1;
			else
			{
				if (log2)
				{
					ii = static_cast<std::size_t>((m_type_inp == 1));

					m_edgesH.set(m_iy + ii, m_ix, -abs(m_edgesH.at(m_iy + ii, m_ix)));
					m_edgesV.set(m_iy, m_ix, -abs(m_edgesV.at(m_iy, m_ix)));

					xk = m_xx + m_hx;
					yk = y2;
					if (m_type_inp == 3)
					{
						m_iy_gran--;
					}

					m_type_inp = 4;
					m_ix_gran++;
				}
				else
				{
					m_edgesV.set(m_iy, m_ix, -abs(m_edgesV.at(m_iy, m_ix)));
					m_edgesV.set(m_iy, m_ix + 1u, -abs(m_edgesV.at(m_iy, m_ix + 1u)));

					xk = x3;
					if (m_type_inp == 1)
					{
						yk = m_yy + m_hy;
						m_iy_gran++;
					}
					else
					{
						yk = m_yy;
						m_iy_gran--;
					}
				}
			}   // The end  log1 & log2;
			break;

		case  2:
		case  4:
			m_edgesV.set(iy_g, ix_g, -abs(m_edgesV.at(iy_g, ix_g)));

			log1 = (m_edgesH.at(m_iy, m_ix) > 0) &&
				((m_value_fun[0] - m_level) * (m_value_fun[1] - m_level) < 0.);

			log2 = (m_edgesH.at(m_iy + 1u, m_ix) > 0) &&
				((m_value_fun[2] - m_level) * (m_value_fun[3] - m_level) < 0.);

			ii = (type_val != 2) ? 1u : 0u;

			log3 = (m_edgesV.at(m_iy, m_ix + ii) > 0) &&
				((m_value_fun[0 + ii] - m_level) * (m_value_fun[3 - ii] - m_level) < 0.);

			ii = static_cast<std::size_t>(log1) +
				static_cast<std::size_t>(log2) +
				static_cast<std::size_t>(log3);

			if (ii == 0u || ii == 2u)
			{
				//fprintf( ErrFile, "\nError_ver:  ixk=%d  iyk=%d  nn=%d\n", ix_k, iy_k, ii );
				return false;
			}

			if (log1)
			{
				x1 = static_cast<coord_type>((m_value_fun[0] - m_level) / (m_value_fun[0] - m_value_fun[1]) * m_hx + m_xx);
			}

			if (log2)
			{
				x2 = static_cast<coord_type>((m_value_fun[3] - m_level) / (m_value_fun[3] - m_value_fun[2]) * m_hx + m_xx);
			}

			if (log1 == log2)
			{
				y3 = static_cast<coord_type>((
					(m_type_inp == 2) ?
					(m_value_fun[0] - m_level) / (m_value_fun[0] - m_value_fun[3]) * m_hy + m_yy :
					(m_value_fun[1] - m_level) / (m_value_fun[1] - m_value_fun[2]) * m_hy + m_yy));
			}

			if (log1)
			{
				if (log2)
				{
					if (fabs(x1 - x2) > get_eps_xy() || fabs(y3 - m_iso_begin.y()) > get_eps_xy())
					{
						if (m_type_inp == 2)
						{
							if ((x1 - x2 > 0.) == (m_hx > 0.))
							{ // x1>x2
								xk = x1;
								yk = m_yy;

								m_type_inp = 3;
								m_ix_gran--;
							}
							else
							{
								xk = x2;
								yk = m_yy + m_hy;

								m_type_inp = 1;
								m_ix_gran--;
								m_iy_gran++;
							}
						}
						else
						{
							if ((x1 - x2 <= 0.) == (m_hx > 0.))
							{ // x1<=x2
								xk = x1;
								yk = m_yy;

								m_type_inp = 3;
							}
							else
							{
								xk = x2;
								yk = m_yy + m_hy;

								m_type_inp = 1;
								m_iy_gran++;
							}
						}
					}
					else
					{
						yk = y3;
						if (m_type_inp == 2)
						{
							xk = m_xx;

							m_ix_gran--;
						}
						else
						{
							xk = m_xx + m_hx;

							m_ix_gran++;
						}
					}
				}   // Є®­Ґж б«гз п log1=log2=1

				else
				{
					ii = static_cast<std::size_t>((m_type_inp == 4));

					m_edgesV.set(m_iy, m_ix + ii, -abs(m_edgesV.at(m_iy, m_ix + ii)));
					m_edgesH.set(m_iy + 1u, m_ix, -abs(m_edgesH.at(m_iy + 1u, m_ix)));

					xk = x1;
					yk = m_yy;
					if (m_type_inp == 2)
					{
						m_ix_gran--;
					}

					m_type_inp = 3;
				}   // Є®­Ґж б«гз п log1=1;  log2=0;
			}

			else
			{  // б«гз © log1=0;
				if (log2)
				{
					ii = static_cast<std::size_t>((m_type_inp == 4));

					m_edgesV.set(m_iy, m_ix + ii, -abs(m_edgesV.at(m_iy, m_ix + ii)));
					m_edgesH.set(m_iy, m_ix, -abs(m_edgesH.at(m_iy, m_ix)));

					xk = x2;
					yk = m_yy + m_hy;
					if (m_type_inp == 2)
					{
						m_ix_gran--;
					}

					m_type_inp = 1;
					m_iy_gran++;
				}   // Є®­Ґж б«гз п log1=0;  log2=1;
				else
				{
					m_edgesH.set(m_iy, m_ix, -abs(m_edgesH.at(m_iy, m_ix)));
					m_edgesH.set(m_iy + 1u, m_ix, -abs(m_edgesH.at(m_iy + 1u, m_ix)));

					yk = y3;
					if (m_type_inp == 2)
					{
						xk = m_xx;

						m_ix_gran--;
					}
					else
					{
						xk = m_xx + m_hx;

						m_ix_gran++;
					}
				}   // Є®­Ґж б«гз п log1=0;  log2=0;
			}
		}  // The end of switch.

		m_iso_end.set_x(xk);
		m_iso_end.set_y(yk);

		return coordinate_control();
	}
};

//isolines::value_type Angle(vector2d const& v1_, vector2d const& v2_) // Результат 0, если один из векторов равен 0. По умолчанию направление положительное.
//{
//	isolines::value_type a = (v1_.x() * v2_.y() - v1_.y() * v2_.x());
//
//	return 0. != a ? std::atan2(a, v1_.dot_product(v2_)) : a;
//}

inline std::size_t integer_divide(coord_type const& a_, coord_type const& b_)
{
	return 0 != b_ ? static_cast<std::size_t>(a_ / b_) : 0;
}

//////////////////////////////////////////////////////////////////////
///
void isolines::build_single(isolines const& source_, vtype const& level_,
							 lines_array* lines_,
							 rings_array* rings_) noexcept(false)
{
	_TRY_BEGIN

	box_params box(&source_, level_);
	box.build(lines_, rings_);

	_CATCH_ALL
	//write log
	_CATCH_END
}


//////////////////////////////////////////////////////////////////////
///
isolines::isolines() noexcept
{
	reset();
}

bool isolines::empty() const noexcept
{
	std::shared_lock<std::shared_mutex> lock(m_guard);
	return value_fun.empty();
}

vtype isolines::get_null_value() const noexcept
{
	std::shared_lock<std::shared_mutex> lock(m_guard);
	return m_null_value;
}

box2d isolines::get_bound() const noexcept
{
	std::shared_lock<std::shared_mutex> lock(m_guard);
	return m_bound;
}

isolines::value_matrix const& isolines::get_values() const noexcept
{
	std::shared_lock<std::shared_mutex> lock(m_guard);
	return value_fun;
}

isolines::edge_matrix const& isolines::get_edges_vertical() const noexcept
{
	std::shared_lock<std::shared_mutex> lock(m_guard);
	return ver_gran;
}

isolines::edge_matrix const& isolines::get_edges_horizontal() const noexcept
{
	std::shared_lock<std::shared_mutex> lock(m_guard);
	return gor_gran;
}

coord_type isolines::get_eps_xy() const noexcept
{
	std::shared_lock<std::shared_mutex> lock(m_guard);
	return m_eps_xy;
}

vtype isolines::get_eps_value() const noexcept
{
	std::shared_lock<std::shared_mutex> lock(m_guard);
	return m_eps_value;
}

isolines::coord_array const& isolines::get_array_xx() const noexcept
{
	std::shared_lock<std::shared_mutex> lock(m_guard);
	return xx;
}

isolines::coord_array const& isolines::get_array_yy() const noexcept
{
	std::shared_lock<std::shared_mutex> lock(m_guard);
	return yy;
}

vtype isolines::get_min_value() const noexcept
{
	std::shared_lock<std::shared_mutex> lock(m_guard);
	return m_min_value;
}

vtype isolines::get_max_value() const noexcept
{
	std::shared_lock<std::shared_mutex> lock(m_guard);
	return m_max_value;
}

void isolines::reset() noexcept
{
	std::unique_lock lock(m_guard);

	value_fun.reset();
	xx.clear();
	yy.clear();
	gor_gran.reset();
	ver_gran.reset();
	m_null_value = default_blank;
	m_dx = m_dy = 0.;
	m_bound.set_null();
	m_min_value = std::numeric_limits<vtype>::max();
	m_max_value = std::numeric_limits<vtype>::min();
	m_eps_xy = default_epsilon_xy;
	m_eps_value = default_epsilon;
}

void isolines::set_holes(multipolygon2d const& holes_)
{
	std::unique_lock lock(m_guard);

	m_holes.clear();
	std::copy_if(holes_.begin(), holes_.end(), std::back_inserter(m_holes),
	[](multipolygon2d::value_type const& poly)
	{
		return boost::geometry::is_valid<polygones_type>(poly);
	});
}

void isolines::init(value_matrix const& matrix_, box2d const& box_, vtype const& vblank_) noexcept(false)
{
	assert(sizeof(vtype) == sizeof(value_matrix::value_type));
	reset();

	_TRY_BEGIN
	std::unique_lock lock(m_guard);

	std::size_t nx = matrix_.cols();
	std::size_t ny = matrix_.rows();
	if(1 >= nx || 1 >= ny)
	{
		throw std::exception("матрица слишком мала");
	}
	value_fun = matrix_;

	assert(nx == value_fun.cols());
	assert(ny == value_fun.rows());

	if(box_.empty())
	{
		throw std::exception("не задана область расчёта");
	}

	m_bound = box_;
	assert(m_bound.is_valid());

	m_dx = static_cast<vtype>(m_bound.width() / static_cast<vtype>(nx - 1u));
	m_dy = static_cast<vtype>(m_bound.height() / static_cast<vtype>(ny - 1u));
	if(min_dxy > m_dx || min_dxy > m_dy)
	{
		throw std::exception("слишком малый шаг координатной сетки");
	}
	m_null_value = vblank_;

	m_min_value = matrix_.min_value(vblank_);
	m_max_value = matrix_.max_value(vblank_);

	if(m_min_value == vblank_ || m_max_value == vblank_)
	{
		throw std::exception("в матрице только пустые значения");
	}

	xx.resize(nx);
	assert(!xx.empty());
	for (std::size_t ix = 0u; ix < nx; ++ix)
	{
		xx[ix] = m_bound.min_x() + static_cast<coord_type>(ix * m_dx);
	}

	yy.resize(ny);
	assert(!yy.empty());
	for (std::size_t iy = 0u; iy < ny; ++iy)
	{
		yy[iy] = m_bound.min_y() + static_cast<coord_type>(iy * m_dy);
	}

	gor_gran.init(ny, nx - 1u);
	ver_gran.init(ny - 1u, nx);

#ifdef _DEBUG
	for (std::size_t ix = 0u; ix < nx; ix++)
	{
		for (std::size_t iy = 0u; iy < ny; iy++)
		{
			assert(value_fun.at(iy, ix) == matrix_.at(iy, ix));
		}
	}
#endif

	analyze_geometry();
	analyze_matrix();
	init_edges();

	_CATCH_ALL
	reset();
	_RERAISE;
	_CATCH_END
}

isolines::value_type isolines::calc_value(point2d const& point_) const
{
	vtype res = m_null_value;
	std::shared_lock<std::shared_mutex> lock(m_guard);

	_TRY_BEGIN
	box2d bound(get_bound());

	if( !empty() && !bound.empty() && boost::geometry::within<point2d, boost_box2d>(point_, bound))
	{
		//точка не в дырке
		if(!boost::geometry::within<point2d, multipolygon2d>(point_, m_holes))
		//if(!boost::geometry::covered_by<point2d, multipolygon2d>(point_, m_holes))
		{
			std::size_t cx = integer_divide(point_.x() - bound.min_x(), m_dx);
			cx = value_fun.cols() <= cx ? value_fun.cols() - 1u : cx;

			std::size_t cy = integer_divide(point_.y() - bound.min_y(), m_dy);
			cy = value_fun.rows() <= cy ? value_fun.rows() - 1u : cy;

			if( edge_matrix::value_type ii = gor_gran.at(cy, cx); ii > 0)
			{
				std::array<vtype, 4> mvalue_fun;
				std::size_t mix = cx;          // type = 1 || 4 || 3
				std::size_t miy = 2 == ii ? (cy - 1u) : cy;

				coord_type mxx = xx[mix];
				coord_type myy = yy[miy];
				coord_type mhx = xx[mix + 1u] - xx[mix];
				coord_type mhy = yy[miy + 1u] - yy[miy];

				mvalue_fun[0] = value_fun.at(miy, mix);
				mvalue_fun[1] = value_fun.at(miy, mix + 1u);
				mvalue_fun[2] = value_fun.at(miy + 1u, mix + 1u);
				mvalue_fun[3] = value_fun.at(miy + 1u, mix);

				spline_param bilin;
				bilin.k_a0 = 0.0 - mvalue_fun[0];
				bilin.k_ax = (mvalue_fun[1] - mvalue_fun[0]) / static_cast<isolines::value_type>(mhx);
				bilin.k_ay = (mvalue_fun[3] - mvalue_fun[0]) / static_cast<isolines::value_type>(mhy);
				bilin.k_axy = (mvalue_fun[0] - mvalue_fun[1] + mvalue_fun[2] - mvalue_fun[3]) / (mhx * mhy);
				bilin.dx = mxx;
				bilin.dy = myy;
				bilin.k_a0 *= -1.;

				res = std::max<vtype>(m_min_value, std::min<vtype>(m_max_value, bilin.calc_value(point_)));
			}
		}
	}
	_CATCH_ALL
	_RERAISE;
	_CATCH_END

	return res;
}

void isolines::build(std::vector<vtype> const& levels_, result_map& result_, bool polygones_) const noexcept(false)
{
	_TRY_BEGIN
		result_.clear();

		if(levels_.empty())
		{
			return;
		}

		std::vector<std::thread> threads;
		for( vtype const& level : levels_)
		{
			if(get_min_value() > level || get_max_value() < level)
			{
				continue;
			}

			//auto& [lines, rings, polygones] = result_[level];
			result_map::mapped_type& out = result_[level];

			lines_result& lines = std::get<0>(out);
			rings_result& rings = std::get<1>(out);

			lines.second = std::to_string(level);
			if(polygones_)
			{
				rings.second = "more or equal than " + std::to_string(level);
			}

			threads.emplace_back(std::thread(algorithm::isolines::build_single, std::ref(*this),
				level, std::addressof(lines.first), polygones_ ? std::addressof(rings.first) : nullptr));
			threads.back().join();
		}

		if(polygones_)
		{
			struct rings_tree
			{
				vtype							level;
				vtype							area;
				ring2d const*					p_ring;
				std::vector<ring2d const*>		inners;

				rings_tree(vtype const& l_, ring2d const* p_) : level(l_), p_ring(p_)
				{
					area = nullptr != p_ ? boost::geometry::area<ring2d>(*p_) : 0.;
				}
			};

			//формирование массива областей
			std::vector<rings_tree> array;
			for( result_map::value_type& item : result_)
			{
				rings_result& rings = std::get<1>(item.second);
				for(auto& r : rings.first)
				{
					array.emplace_back(rings_tree(item.first, std::addressof(r)));
				}
			}

			if(!array.empty())
			{
				const auto minmax = std::minmax_element(array.begin(), array.end(),
					[](rings_tree const& lh, rings_tree const& rh)
					{
						return std::less<vtype>{}(lh.level, rh.level);
					});
				vtype min_level = minmax.first->level;
				vtype max_level = minmax.second->level;

				//добавляем расчётную область
				ring2d full_area;
				boost::geometry::convert<boost_box2d, ring2d>(get_bound(), full_area);

				array.emplace_back(rings_tree(get_null_value(), std::addressof(full_area)));

				//сортировка по площади
				std::sort(array.begin(), array.end(), [](rings_tree const& lh, rings_tree const& rh)
				{
					return std::less<vtype>{}(lh.area, rh.area);
				});

				//поиск какая область находится внутри другой(только одно)
				std::vector<rings_tree>::iterator loop = array.begin(), end = array.end();
				while(end != loop)
				{
					std::vector<rings_tree>::iterator outer = (loop + 1);
					while(end != outer)
					{
	//					bool within = boost::geometry::within<ring2d, ring2d>(*(*loop)->p_ring, *(*outer)->p_ring);
	//					bool over = boost::geometry::within<ring2d, ring2d>(*(*loop)->p_ring, *(*outer)->p_ring);
	//					if(over || within)
						if(boost::geometry::within<ring2d, ring2d>(*(*loop).p_ring, *(*outer).p_ring))
						{
							(*outer).inners.push_back((*loop).p_ring);
							break;
						}

						++outer;
					}

					++loop;
				}

				//формируем области
				for(rings_tree const& tr_ring : array)
				{
					if(get_null_value() == tr_ring.level)//
						continue;

					polygones_result& polygones = std::get<2>(result_[tr_ring.level]);
					if(tr_ring.level >= max_level)
					{
						polygones.second = "more than " + std::to_string(tr_ring.level);
					}
					else if(min_level > tr_ring.level)
					{
						polygones.second = "less than " + std::to_string(min_level);
					}
					else //if(0.0 != tr_ring.level)
					{
						vtype next = *std::upper_bound(levels_.begin(), levels_.end(), tr_ring.level);
						polygones.second = "from " + std::to_string(tr_ring.level) + " to " + std::to_string(next);
					}

					polygon2d outer;
					outer.outer().assign(tr_ring.p_ring->begin(), tr_ring.p_ring->end());
					//boost::geometry::correct<polygon2d>(outer);
					assert(boost::geometry::is_valid<polygon2d>(outer));
					if(!tr_ring.inners.empty())
					{
						multipolygon2d clip;
						clip.resize(tr_ring.inners.size());
						multipolygon2d::iterator inner_poly = clip.begin();
						for(auto const* inner : tr_ring.inners)
						{
							(*inner_poly).outer().assign(inner->begin(), inner->end());
							++inner_poly;
						}

						//bool xx = boost::geometry::within<multipolygon2d, polygon2d>(clip, outer);
						multipolygon2d res;
						boost::geometry::difference<polygon2d, multipolygon2d, multipolygon2d>(outer, clip, res);
						for(auto const& p : res)
						{
							polygones.first.emplace_back(std::make_unique<polygones_type>(p));
							//polygones.first.push_back(p);
						}
					}
					else
					{
						polygones.first.emplace_back(std::make_unique<polygones_type>(outer));
						//polygones.first.push_back(outer);
					}
				}
			}
		}

		//вырежем дырки
		if(!m_holes.empty())
		{
			for( auto& level_item : result_ )
			{
				auto& [lines, rings, polygones] = level_item.second;

				lines_array clip_lines;
				for( auto& pline : lines.first)
				{
					if(!(bool)pline)
					{
						continue;
					}

					if( boost::geometry::intersects<line2d, multipolygon2d>((*pline), m_holes))
					{
						multiline2d res;
						boost::geometry::difference<line2d, multipolygon2d, multiline2d>((*pline), m_holes, res);

						pline.reset();
						for(auto const& line : res)
						{
							clip_lines.emplace_back(std::make_unique<lines_type>(line));
						}
					}
				}
				//удаляем обрезанные
				lines.first.erase(std::remove_if(lines.first.begin(), lines.first.end(),
				[](isolines::lines_array::value_type const& pline)
				{
					return !((bool)pline);
				}), lines.first.end());
				//добавляем новые
				for(auto& line : clip_lines)
				{
					lines.first.push_back(std::move(line));
				}

				polygones_array clip_poly;
				for( auto& ppoly : polygones.first)
				{
					if(!(bool)ppoly)
					{
						continue;
					}

					if(boost::geometry::intersects<polygon2d, multipolygon2d>(*ppoly, m_holes))
					{
						multipolygon2d res;
						boost::geometry::difference<polygon2d, multipolygon2d, multipolygon2d>(*ppoly, m_holes, res);

						ppoly.reset();
						for(auto const& poly : res)
						{
							clip_poly.emplace_back(std::make_unique<polygones_type>(poly));
						}

					}
				}
				//удаляем обрезанные
				polygones.first.erase(std::remove_if(polygones.first.begin(), polygones.first.end(),
				[](isolines::polygones_array::value_type const& ppoly)
				{
					return !((bool)ppoly);
				}), polygones.first.end());

				for(auto& poly : clip_poly)
				{
					polygones.first.push_back(std::move(poly));
				}
			}
		}


	_CATCH_ALL
	_RERAISE;
	_CATCH_END
}

void isolines::init_edges()
{
	std::size_t nx = value_fun.cols();
	std::size_t ny = value_fun.rows();

	edge_matrix::value_type edge_value, log0, log1;

	for (std::size_t iy = 0u; iy < ny; ++iy)
	{
		log0 = (iy > 0);
		log1 = (iy < ny - 1u);
		for (std::size_t ix = 0u; ix < nx - 1u; ++ix)
		{
			edge_value = m_null_value;
			if (m_null_value != value_fun.at(iy, ix)  && m_null_value != value_fun.at(iy, ix + 1u))
			{
				if (log1 && m_null_value != value_fun.at(iy + 1u, ix) && m_null_value != value_fun.at(iy + 1u, ix + 1u))
				{
					edge_value += 1;
				}
				if (log0 && m_null_value != value_fun.at(iy - 1u, ix) && m_null_value != value_fun.at(iy - 1u, ix + 1u))
				{
					edge_value += 2;
				}
			}
			gor_gran.set(iy, ix, edge_value);
		}
	}

	for (std::size_t ix = 0u; ix < nx; ++ix)
	{
		log0 = (ix > 0);
		log1 = (ix < nx - 1);
		for (std::size_t iy = 0; iy < ny - 1u; iy++)
		{
			edge_value = m_null_value;
			if (m_null_value != value_fun.at(iy, ix) && m_null_value != value_fun.at(iy + 1u, ix) )
			{
				if (log1 && m_null_value != value_fun.at(iy, ix + 1u) && m_null_value != value_fun.at(iy + 1u, ix + 1u))
				{
					edge_value += 1;
				}
				if (log0 && m_null_value != value_fun.at(iy, ix - 1u) && m_null_value != value_fun.at(iy + 1u, ix - 1u))
				{
					edge_value += 2;
				}
			}
			ver_gran.set(iy, ix, edge_value);
		}
	}
}

void isolines::analyze_matrix() noexcept(false)
{
	value_type eps_fun = (value_type)0.;
	value_type tmp;

	tmp = (fabs(m_max_value) + fabs(m_min_value)) / (value_type)2.;
	if (tmp > (m_max_value - m_min_value) * default_multiplicator * (value_type)10.)
	{
		throw(std::exception("Близкие значения сеточного массива"));
	}

	eps_fun = (m_max_value - m_min_value) / (default_multiplicator * (value_type)10.);
	m_eps_value = std::max<vtype>(m_eps_value, eps_fun);
}

void isolines::analyze_geometry() noexcept(false)
{
	if (!m_bound.empty())
	{
		coord_type eps_xy = 0.;

		coord_type tmp1 = m_bound.width();
		coord_type tmp2 = m_bound.height();

		if ((tmp1 > tmp2 * 1.e3) || (tmp1 < tmp2 * 1.e-3))
		{
			throw std::exception("значения координат матрицы слишком близки");
		}

		tmp1 = (fabs(m_bound.max_x()) + fabs(m_bound.min_x()) + fabs(m_bound.max_y()) + fabs(m_bound.min_y())) / 4.;
		tmp2 = (fabs(m_bound.max_x() - m_bound.min_x()) + fabs(m_bound.max_y() - m_bound.min_y())) / 2.;
		if (tmp1 < tmp2 * 1.e5)
		{
			eps_xy = static_cast<coord_type>(tmp2) / 1.e8;
		}
		else
		{
			throw std::exception("значения координат матрицы слишком близки");
		}

		m_eps_xy = std::min<coord_type>(m_eps_xy, eps_xy);
	}
}

};
