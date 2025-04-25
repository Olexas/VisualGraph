#pragma once

#include <QMainWindow>
#include <QTableWidget>
#include <QHeaderView>
#include <QVBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>

#include "Graph.h"

class TableWidget : public QTableWidget {
public:
    TableWidget(const QVector<QVector<int>>& matrix, QWidget* parent = nullptr, double k = 0.026, const QColor& highlightColor = QColor(173, 216, 230));
    void RegenarateTable(const QVector<QVector<int>>& matrix, QWidget* parent, double k = 0.026);

private:
    QColor highlightColor;
    QStringList getHeaderLabels(int size);
};

class MainWindow : public QMainWindow {
    Q_OBJECT

private:
    Graph* graph;
    TableWidget* adjacency_table;
    TableWidget* weight_table;
    TableWidget* min_path_table;
    TableWidget* max_path_table;
    TableWidget* capacity_table;
    TableWidget* cost_table;

    QLineEdit* vertex_input;
    QLineEdit* shimbell_edges_input;
    QLineEdit* warshall_first_input;
    QLineEdit* warshall_second_input;

    QVBoxLayout* graph_layout;

    QWidget* inputWidget;
    QTabWidget* table_widget;

    QPushButton* generate_button;
    QPushButton* calculate_button;
    QPushButton* warshall_calculate_button;

    QLabel* warshall_label;
    QLabel* DFS_label;
    QLabel* bellmanford_label;
    QLabel* iterations_label;  
    QLabel* flow_label;
    QLabel* min_cost_flow_label;

    const double double_tables_k = 0.045;

public:
    MainWindow();

public slots:
    void generateGraph();
    void calculateShimbell();
    void onTabChanged(int index);
    void calculateWarshall();
    void calculateFlow();
};