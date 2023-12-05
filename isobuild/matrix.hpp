#ifndef _MATRIX_HPP_
#define _MATRIX_HPP_

#include <assert.h>
#include <vector>
#include <algorithm>

#include <boost/json.hpp>

namespace math {

///////////////////////////////////////////////////////////////////////////////////////////////////
//
template<typename Type>
class matrix2d
{
public:
	using value_type						= Type;
	using array_type						= typename std::vector<value_type>;
	using size_type							= typename array_type::size_type;
	using iterator							= typename array_type::iterator;
	using const_iterator					= typename array_type::const_iterator;

	friend inline boost::json::value const& operator>>(boost::json::value const& v_, matrix2d<Type>& mx_)
	{
		if (boost::json::object const* obj = v_.if_object(); nullptr != obj)
		{
			if (boost::json::value const* v = obj->if_contains("rows"); nullptr != v)
			{
				assert(v->is_primitive());
				mx_.m_rows = v->to_number<size_type>();
			}
			if (boost::json::value const* v = obj->if_contains("cols"); nullptr != v)
			{
				assert(v->is_primitive());
				mx_.m_cols = v->to_number<size_type>();
			}
			if (boost::json::value const* v = obj->if_contains("data"); nullptr != v)
			{
				assert(v->is_array());
				if( boost::json::array const* arr = v->if_array(); arr && !arr->empty())
				{
					assert(arr->size() == mx_.cols() * mx_.rows());
					mx_.m_data.resize(arr->size());

					typename matrix2d<Type>::size_type loop = 0;
					for(auto const& item : (*arr))
					{
						mx_.m_data[loop] = static_cast<Type>(item.as_double());
						++loop;
					}
				}

			}
		}
		else
		{
			mx_.m_data.clear();
		}

		return v_;
	}

	friend inline boost::json::value& operator<<(boost::json::value& v_, matrix2d<Type> const& mx_)
	{
		if(!mx_.empty())
		{
			boost::json::object& obj = v_.is_object() ? v_.as_object() : v_.emplace_object();
			obj["cols"] = mx_.cols();
			obj["rows"] = mx_.rows();
			{
				boost::json::array& arr = obj["data"].emplace_array();
				arr.resize(mx_.m_data.size());

				boost::json::array::size_type loop = 0;
				for(auto const& item : mx_.m_data)
				{
					arr[loop] = item;
					++loop;
				}
			}
		}
		return v_;
	}

	inline matrix2d() noexcept : m_rows(0), m_cols(0) { }
	explicit matrix2d(size_type const& rows_/*ny*/, size_type const& cols_/*nx*/) noexcept(false)
	{
		init(rows_, cols_);
	}

	matrix2d(matrix2d const&) = default;
	matrix2d& operator=(matrix2d const&) = default;

	matrix2d(matrix2d&&) noexcept = default;
	matrix2d& operator=(matrix2d&&) noexcept = default;

	inline bool empty() const noexcept { return m_data.empty(); }
	inline size_type size() const noexcept { return m_data.size(); }

	inline iterator begin() noexcept { return m_data.begin(); }
	inline const_iterator begin() const noexcept { return m_data.begin(); }
	inline iterator end() noexcept { return m_data.end(); }
	inline const_iterator end() const noexcept { return m_data.end(); }

	inline value_type const* data() const noexcept { return m_data.data(); }
	inline value_type* data() noexcept { return m_data.data(); }

	//ny
	inline size_type rows() const noexcept { return m_rows; }
	//nx
	inline size_type cols() const noexcept { return m_cols; }

	value_type const operator() (size_type const row_, size_type const col_) const
	{
		return at(row_, col_);
	}

	value_type const& operator[](size_type const idx_) const
	{
		return m_data[idx_];
	}

	value_type& operator[](size_type const idx_)
	{
		return m_data[idx_];
	}

	void reset() noexcept
	{
		m_data.clear();
		m_rows = m_cols = 0;
	}

	bool init(size_type const rows_/*ny*/, size_type const cols_/*nx*/) noexcept(false)
	{
		if(0 < rows_ && 0 < cols_)
		{
			_TRY_BEGIN
			m_data.resize(rows_ * cols_, (value_type)0);

			m_rows = rows_;
			m_cols = cols_;
			_CATCH_ALL
			m_rows = m_cols = 0;
			m_data.clear();
			_RERAISE;
			_CATCH_END
		}

		return !empty();
	}

	inline bool into_matrix(size_type const row_, size_type const col_) const noexcept
	{
		return (m_rows > row_ && m_cols > col_);
	}

	value_type at(size_type const row_, size_type const col_) const
	{
		size_type index = row_ * m_cols + col_;
		assert(index < m_data.size());

		if (index < m_data.size())
		{
			return m_data[index];
		}

		return std::numeric_limits<value_type>::infinity();
	}

	bool set(size_type const row_, size_type const col_, value_type const& value_)
	{
		size_type index = row_ * m_cols + col_;
		if (index < m_data.size())
		{
			m_data[index] = value_;

			return true;
		}

		return false;
	}

	value_type min_value(value_type const& vblank_) const noexcept
	{
		const_iterator min = std::min_element(m_data.begin(), m_data.end(),
		[&vblank_](value_type const& lhs_, value_type const& rhs_)
		{
			if (lhs_ == vblank_ || rhs_ == vblank_)
			{
				return false;
			}

			return lhs_ < rhs_;
		});

		return m_data.end() != min ? (*min) : vblank_;
	}

	value_type max_value(value_type const& vblank_) const noexcept
	{
		const_iterator max = std::max_element(m_data.begin(), m_data.end(),
		[&vblank_](value_type const& lhs_, value_type const& rhs_)
		{
			if (lhs_ == vblank_ || rhs_ == vblank_)
			{
				return false;
			}

			return lhs_ < rhs_;
		});

		return m_data.end() != max ? (*max) : vblank_;
	}

private:
	size_type								m_rows;
	size_type								m_cols;
	array_type								m_data;
};

};//namespace lpp

#endif // _MATRIX_HPP_
