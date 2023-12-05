#include <filesystem>
#include <fstream>

#include <QPainter>
#include <QMatrix>
#include <QWheelEvent>

#include "mainwindow.h"
#include "ui_mainwindow.h"

namespace
{
static constexpr std::string_view tag_blank("blank");
static constexpr std::string_view tag_bound("bound");
static constexpr std::string_view tag_matrix("matrix");

inline QPointF point_d2(geo::point2d const& p_)
{
	return QPointF(static_cast<qreal>(p_.x()), static_cast<qreal>(p_.y()));
}

std::unique_ptr<QPainterPath> create_geometry(geo::line2d const& line_)
{
	std::unique_ptr<QPainterPath> res;
	if(!line_.empty())
	{
		res = std::make_unique<QPainterPath>(point_d2(line_.front()));
		assert((bool)res);
		if((bool)res)
		{
			std::for_each(line_.begin() + 1, line_.end(),
				[&res](geo::point2d const& pt)
				{
					res->lineTo(point_d2(pt));
				});
			res->moveTo(point_d2(line_.front()));
			res->closeSubpath();
		}
	}

	return res;
}

std::unique_ptr<QPainterPath> create_geometry(geo::ring2d const& ring_)
{
	std::unique_ptr<QPainterPath> res;
	if(!ring_.empty())
	{
		res = std::make_unique<QPainterPath>(point_d2(ring_.front()));
		assert((bool)res);
		if((bool)res)
		{
			std::for_each(ring_.begin() + 1, ring_.end(),
				[&res](geo::point2d const& pt)
				{
					res->lineTo(point_d2(pt));
				});
			res->closeSubpath();
		}
	}

	return res;
}

std::unique_ptr<QPainterPath> create_geometry(geo::box2d const& box_)
{
	geo::ring2d r;
	boost::geometry::convert<geo::boost_box2d, geo::ring2d>(box_, r);
	return create_geometry(r);
}

std::unique_ptr<QPainterPath> create_geometry(geo::polygon2d const& poly_)
{
	std::unique_ptr<QPainterPath> res(create_geometry(poly_.outer()));
	if((bool)res)
	{
		for(auto const& ring : poly_.inners())
		{
			std::unique_ptr<QPainterPath> inner = create_geometry(ring);
			if((bool)inner)
			{
				res->addPath(*inner);
			}
		}
	}

	return res;
}

};



MainWindow::MainWindow(QWidget *parent)
	: QMainWindow(parent)
	, use_holes(true)
	, with_areas(false)
	, ui(new Ui::MainWindow)
	, mscale(125000.f)
{
	ui->setupUi(this);

	ui->actWithAreas->setChecked(with_areas);
	ui->actUseHoles->setChecked(use_holes);

	m_holes[0].outer().assign({{4423.07, 6701.44}, {7864.53, 1786.53}, {6553.89, 868.80}, {3112.43, 5783.71}, {4423.07, 6701.44}});
	m_holes[1].outer().assign({{6047.95, 8475.97}, {8510.59, 5323.93}, {7249.77, 4338.87}, {4787.13, 7490.91}, {6047.95, 8475.97}});
}

MainWindow::~MainWindow()
{
	delete ui;
}

void MainWindow::build_isolines()
{
	if(m_matrix.empty())
	{
		std::filesystem::path json(L"../test_app/matrix.json");
		std::ifstream fs;
		fs.open(json, std::ios::binary | std::ios::ate);
		if(fs.is_open())
		{
			auto size = fs.tellg();
			std::string buffer(size, '\0');
			fs.seekg(0, fs.beg);
			if (fs.read(buffer.data(), size))
			{
				boost::json::error_code ec;
				boost::json::value jv = boost::json::parse(buffer, ec);
				if(!ec)
				{
					if (boost::json::object const* obj = jv.if_object(); nullptr != obj)
					{
						if (boost::json::value const* v = obj->if_contains(::tag_blank.data()); nullptr != v)
						{
							assert(v->is_primitive());
							m_blank = v->as_double();
						}
						if (boost::json::value const* v = obj->if_contains(::tag_bound.data()); nullptr != v)
						{
							assert(v->is_array());
							(*v) >> m_bound;
						}
						if (boost::json::value const* v = obj->if_contains(::tag_matrix.data()); nullptr != v)
						{
							assert(v->is_object());
							(*v) >> m_matrix;
						}
					}
				}
			}
			fs.close();
		}
	}

	if(!m_matrix.empty() && !m_bound.empty())
	{
		geo::algorithm::isolines iso;
		iso.init(m_matrix, m_bound, m_blank);

		if(use_holes)
		{
			geo::multipolygon2d pp;

			for(auto const& h : m_holes) { pp.push_back(h); }

			iso.set_holes(pp);
		}

		//iso.build({ 0.84 }, res, with_areas);
		iso.build({ 0.36, 0.37, 0.39, 0.418, 0.5, 0.84, 1.0, 8.0, 10.0 }, res, with_areas);
		//iso.build({ 0.36, 0.363, 0.367, 0.39, 0.401, 0.418, 0.5, 0.84, 1.0, 8.0, 10.0, 15.0, 18.0, 20.0, 25.0 }, res, with_areas);

		m_linespaths.clear();
		m_ringspaths.clear();
		m_polypaths.clear();

		int loop = 0;
		for(auto const& item : res)
		{
			++loop;
			if(0 == item.first)
				continue;

			auto const& [lll, rrr, ppp] = item.second;

			if(with_areas)
			{
				//for(auto const& ring : rrr.first)
				//{
				//	m_ringspaths.push_back(std::pair(create_geometry(ring),
				//		QColor(00, std::clamp<int>(255 - 24 * loop, 32, 255), 00, 125)));
				//}

				for(auto const& poly : ppp.first)
				{
					m_polypaths.push_back(std::pair(create_geometry(*poly),
						QColor(std::clamp<int>(255 - 24 * loop, 32, 255), 00, 00, 255)));
				}
			}

			for(auto const& line : lll.first)
			{
				m_linespaths.push_back(std::pair(create_geometry(*line),
					QColor(00, 00, std::clamp<int>(255 - 24 * loop, 00, 00), 255)));
			}
		}

	}
}

void MainWindow::paintEvent(QPaintEvent *event)
{
	QPainter painter(this);
	painter.setRenderHint(QPainter::Antialiasing);

	QRectF qbound(nullptr != event ? event->rect() : QRectF(0., 0.,
		static_cast<qreal>(width()),
		static_cast<qreal>(height())));
	painter.fillRect(qbound, Qt::white);

	if(!res.empty())
	{
		if(use_holes && m_holepath.empty())
		{
			for(auto const& h : m_holes)
			{
				m_holepath.emplace_back(create_geometry(h));
			}
		}
		if(!(bool)m_boundpath) {m_boundpath = create_geometry(m_bound); }

		float scaleX = static_cast<float>(physicalDpiX()) / (0.0254f * get_scale());
		float scaleY = static_cast<float>(physicalDpiY()) / (0.0254f * get_scale());

		QPointF center(qbound.center());

		QMatrix m(scaleX, 0, 0, scaleY, center.x(), center.y());
		QPointF mc(m.map(-get_center()));

		QMatrix matrix(m.m11(), m.m12(), m.m21(), m.m22(), mc.x(), mc.y());

		float width = 1.27 / scaleX;

		painter.setWorldMatrix(matrix);
		painter.begin(painter.device());


		painter.strokePath(*m_boundpath, QPen(QColor(55, 00, 00, 125), width, Qt::DashLine));
		for(auto const& p : m_ringspaths)
		{
			if((bool)p.first)
			{
				painter.fillPath(*(p.first), QBrush(p.second, Qt::Dense4Pattern));
				painter.strokePath(*(p.first), QPen(p.second, width));
			}
		}

		for(auto const& p : m_polypaths)
		{
			if((bool)p.first)
			{
				painter.fillPath(*(p.first), QBrush(p.second));
				painter.strokePath(*(p.first), QPen(p.second, width));
			}
		}

		for(auto const& p : m_linespaths)
		{
			if((bool)p.first)
			{
				painter.strokePath(*(p.first), QPen(p.second, width));
			}
		}

		if(use_holes)
		{
			for(auto const& path : m_holepath)
			{
				if((bool)path)
				{
					painter.fillPath(*path, QBrush(QColor(00, 255, 00, 255), Qt::CrossPattern));
					painter.strokePath(*path, QPen(QColor(00, 106, 00, 200), width, Qt::DotLine));
				}
			}
		}

		painter.end();
	}
}

void MainWindow::wheelEvent(QWheelEvent *event)
{
	if( Qt::NoButton == event->buttons())
	{
		int delta = event->delta();
		if(delta != 0)
		{
			double scale = get_scale();
			mscale = delta >= 0 ? scale * 1.1 : scale / 1.1;

			repaint();
			event->accept();
		}
	}
}

void MainWindow::on_actWithAreas_triggered(bool checked)
{
	with_areas = checked;

	res.clear();
	build_isolines();

	repaint();
}


void MainWindow::on_actUseHoles_triggered(bool checked)
{
	use_holes = checked;

	res.clear();
	build_isolines();

	repaint();
}

void MainWindow::on_actionSave_triggered()
{
	boost::json::value save;
	if(!res.empty())
	{
		boost::json::object& out = save.emplace_object();
		for(auto const& item : res)
		{
			std::string tag = std::to_string(item.first);
			boost::json::object& jitem = out[tag].emplace_object();

			auto const& [lines, rings, polygones] = item.second;

			boost::json::array& ls = jitem["lines"].emplace_array();
			for(auto const& line : lines.first)
			{
				if(!line->empty())
				{
					boost::json::value ij;
					ij << (*line);

					ls.emplace_back(ij);
				}
			}

			if(with_areas)
			{
				boost::json::array& rs = jitem["rings"].emplace_array();
				for(auto const& ring : rings.first)
				{
					if(!ring.empty())
					{
						boost::json::value ij;
						ij << ring;

						rs.emplace_back(ij);
					}
				}

				boost::json::array& pps = jitem["polygones"].emplace_array();
				for(auto const& poly : polygones.first)
				{
					if(!poly->outer().empty())
					{
						boost::json::value ij;
						ij << (*poly);

						pps.emplace_back(ij);
					}
				}
			}
		}

		std::filesystem::path json(L"../test_app/result.json");
		std::ofstream outf;
		outf.open(json, std::ifstream::out | std::ifstream::trunc);
		if(outf.is_open())
		{
			std::string res_data = boost::json::serialize(save);
			outf << res_data;
			outf.close();
		}
	}
}


