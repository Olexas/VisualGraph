#pragma once

#include <QMainWindow>
#include <QTableWidget>
#include <QHeaderView>
#include <QVBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QHash>

#include "Graph.h"

class TableWidget : public QTableWidget {
public:
    TableWidget(const QVector<QVector<int>>& matrix, QWidget* parent = nullptr, double k = 0.020, const QColor& highlightColor = QColor(173, 216, 230));
    void RegenerateTable(const QVector<QVector<int>>& matrix, QWidget* parent, double k = 0.020);

private:
    QColor highlightColor;
    QStringList GenerateHeaderLabels(int size);
};

class MainWindow : public QMainWindow {
    Q_OBJECT

private:
    struct TableData {
        TableWidget* widget;
        QLabel* label;
        const QVector<QVector<int>>& (Graph::* getter)() const;
    };

    Graph* graph;
    QVBoxLayout* graph_layout;

    // TODO: input_widget
    QWidget* inputWidget;
    QTabWidget* tab_widget;

    QHash<QString, TableData> tables;
    QHash<QString, QLabel*> labels;

    QLineEdit* vertex_input;
    QLineEdit* shimbell_edges_input;
    QLineEdit* warshall_first_input;
    QLineEdit* warshall_second_input;

    QPushButton* generate_button;
    QPushButton* calculate_button;
    QPushButton* warshall_calculate_button;

    const double double_tables_k = 0.027;
    int current_tab;

    void setupUI();
    void createInputWidget();
    void createMatricesTab(QWidget* tab);
    void createShimbellTab(QWidget* tab);
    void createWarshallTab(QWidget* tab);
    void createFlowTab(QWidget* tab);
    void createSpanningTreeTab(QWidget* tab);
    void createEulerianTab(QWidget* tab);

    void initializeTables();
    void initializeLabels();

    void updateDynamicLabels();
    void setupConnections();
    void setupTabWidget();

    QString updatePruferInfo();
    QString joinQVector(const QVector<int>& vec) const;
    QString formatDecodedEdges(const QVector<QPair<QPair<int, int>, int>>& edges) const;
    void updateGraphViewSettings();
     
public:
    MainWindow();

private slots:
    void generateGraph();
    void calculateShimbell();
    void onTabChanged(int index);
    void calculateWarshall();
    void calculateFlow();
    void findEulerianCycle();
    void findHamiltonianCycle();
};