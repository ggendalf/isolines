#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

#include "isolines.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
	Q_OBJECT

public:
	MainWindow(QWidget *parent = nullptr);
	~MainWindow();

	void build_isolines();

	bool									use_holes;
	bool									with_areas;

private:
	Ui::MainWindow		*ui;
	float				mscale;
protected:
	void paintEvent(QPaintEvent *event) override final;
	void wheelEvent(QWheelEvent *event) override final;


	geo::algorithm::isolines::result_map	res;

	float get_scale() const { return mscale; }
	QPointF get_center() const { return QPointF(3100.f, 5350.f);}

	geo::algorithm::isolines::value_matrix	m_matrix;
	geo::box2d								m_bound;
	double									m_blank;
	geo::polygon2d							m_holes[2];

	std::vector<std::pair<std::unique_ptr<QPainterPath>, QColor>> m_linespaths;
	std::vector<std::pair<std::unique_ptr<QPainterPath>, QColor>> m_ringspaths;
	std::vector<std::pair<std::unique_ptr<QPainterPath>, QColor>> m_polypaths;
	std::vector<std::unique_ptr<QPainterPath>> m_holepath;
	std::unique_ptr<QPainterPath>			m_boundpath;

private slots:
	void on_actWithAreas_triggered(bool checked);
	void on_actUseHoles_triggered(bool checked);
	void on_actionSave_triggered();
};
#endif // MAINWINDOW_H
