#pragma once

#include <QWidget>
#include <QPainter>
#include <QMainWindow>
#include <QVector>
#include <QList>
#include <iostream>
#include <queue>

class Graph : public QMainWindow {
private:
    int vertexes_count;
    int edges_count;
    QVector<QVector<int>>* adjacency_matrix;
    QVector<QVector<int>>* weight_matrix;
    QVector<QVector<int>>* capacity_matrix;
    QVector<QVector<int>>* cost_matrix;

    QList<int>* adjLists;

    QVector<int> shortestPath;

    void GenerateGraph();
    bool CheckGraph() const;
    void dfs(int v, int num, QVector<int>& component) const;

    QVector<QVector<int>>* getWayByShimbell(int edges_count, bool is_min = true) const;
    void FillRandNumbersVector(QVector<int>& vector, int max_num = 0, double lambda = 3, double a = 0.9) const;
    double getRandNumber(double lambda, double a) const;

    bool showFlowInfo = false;

public:
    Graph(int vertexes_count, QWidget* parent = nullptr);

    QVector<int>* GenerateDegrees(int vertexes_count);
    void FillAdjacencyMatrix(QVector<int>& degrees);
    void FillMatrix(QVector<QVector<int>>& matrix, int max_value, double lambda, double a, bool can_negative = false);
    void FillWeightMatrix();
    void FillCapacityMatrix();
    void FillCostMatrix();

    QVector<QVector<int>>* getMinWayByShimbell(int edges_count) const;
    QVector<QVector<int>>* getMaxWayByShimbell(int edges_count) const;

    int FindPathCount(int first_vertex, int second_vertex) const;

    int DFS(QString& way, int vertex2, int vertex1 = -1, QVector<QVector<bool>>* visited = nullptr) const;

    int BellmanFord(int startVertex, int endVertex, QVector<int>& distances, QVector<int>& path) const;
    void clearHighlight();

    int FordFulkerson(int source, int sink);
    std::pair<int, int> FindMinCostFlow(int source, int sink, int max_flow = -1);

    void paintEvent(QPaintEvent* event) override;

    void drawGraph(QPainter& painter);
    void drawEdge(QPainter& painter, int i, int j, int centerX, int centerY,
        int ellipseWidth, int ellipseHeight, int vertexRadius, bool isHighlighted, bool isForward) const;
    void drawWeight(QPainter& painter, int i, int j, int centerX, int centerY,
        int ellipseWidth, int ellipseHeight, int vertexRadius) const;
    void setShowFlowInfo(bool show);

    const QVector<QVector<int>>& getAdjacencyMatrix() const;
    const QVector<QVector<int>>& getWeightMatrix() const;
    const int getVertexesCount() const;
    const QVector<QVector<int>>& getCapacityMatrix() const;
    const QVector<QVector<int>>& getCostMatrix() const;
};


void CoutMatrix(QVector<QVector<int>>& matrix);
void CoutMatrix(QVector<QVector<bool>>& matrix);