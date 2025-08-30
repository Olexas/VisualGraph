#pragma once

#include <QWidget>
#include <QPainter>
#include <QMainWindow>
#include <QVector>
#include <QList>
#include <iostream>
#include <queue>
#include <QQueue>
#include <QStack>
#include <QSet>
#include <QTextStream>

class Graph : public QMainWindow {
private:
    int vertexes_count;
    int edges_count;
    QVector<QVector<int>>* adjacency_matrix;
    QVector<QVector<int>>* weight_matrix;
    QVector<QVector<int>>* capacity_matrix;
    QVector<QVector<int>>* cost_matrix;
    QVector<QVector<int>>* kirchhoff_matrix;

    QVector<QVector<int>>* empty_matrix;

    QVector<QPair<int, int>>* adj_vector;

    QVector<int> shortestPath;

    QVector<QVector<int>>* all_cycles;
    QVector<int>* cycle_weights;
    bool is_hamiltonian;
    QVector<int> hamiltonianAddedEdgesWeights;
    
    int spanning_trees_count;
    int mst_weight;
    QVector<QPair<int, int>>* mst_vector;
    QVector<int>* colors;
    int unique_colors;

    void GenerateGraph();
    bool CheckGraph() const;
    void dfs(int v, int num, QVector<int>& component) const;

    QVector<QVector<int>>* getWayByShimbell(int edges_count, bool is_min = true) const;
    void FillRandNumbersVector(QVector<int>& vector, int max_num = 0, double lambda = 3, double a = 0.9) const;
    double getRandNumber(double lambda, double a) const;

    double calculateDeterminant(QVector<QVector<double>>& matrix) const;
    void Boruvka();
    QVector<QVector<int>>* GetNeorAdjacencyMatrix() const;

    QVector<QPair<int, int>> addedEdges;
    QVector<QPair<int, int>> removedEdges;

    void findAllHamiltonianCycles(QVector<QVector<int>>& allCycles, QVector<int>& cycleWeights) const;
    void hamiltonianDFS(int start, int current, QVector<bool>& visited,
        QVector<int>& path, QVector<QVector<int>>& allCycles,
        QVector<int>& cycleWeights, QVector<QVector<int>>& temp_adj, QVector<QVector<int>>& temp_weight) const;
    void removeDuplicateCycles(QVector<QVector<int>>& allCycles, QVector<int>& cycleWeights) const;
    QString cycleToHash(const QVector<int>& cycle) const;

    QVector<QPair<int, int>> hamiltonianAddedEdges;

    bool showFlowInfo = false;
    bool neor_flag = false;
    bool show_MST = false;
    bool show_eulerian = false;

    void drawGraph(QPainter& painter);
    void drawEdge(QPainter& painter, int i, int j, int centerX, int centerY,
        int ellipseWidth, int ellipseHeight, int vertexRadius, bool isHighlighted, bool isForward) const;
    void drawWeight(QPainter& painter, int i, int j, int centerX, int centerY,
        int ellipseWidth, int ellipseHeight, int vertexRadius,
        int weight, bool isAddedEdge) const;
    void drawMST(QPainter& painter, const QVector<QPair<int, int>>& mstEdges) const;

public:
    Graph(int vertexes_count, QWidget* parent = nullptr);
    ~Graph();

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

    int BellmanFord(int startVertex, int endVertex, QVector<QString>& distances_result, QVector<int>& path) const;
    void clearHighlight();

    int FordFulkerson(int source, int sink);
    std::pair<int, int> FindMinCostFlow(int source, int sink, QVector<QVector<int>>& flow, int max_flow = -1);

    void GenerateSpanningTreesCount();
    QPair<QVector<int>, QVector<int>> encodePrufer();
    QVector<QPair<QPair<int, int>, int>> decodePrufer(const QVector<int>& pruferCode, const QVector<int>& weights) const;
    void GreedyColoring();  

    bool isEulerian() const;
    void makeEulerian();
    QVector<QPair<int, int>> findEulerianCycle() const;
    void clearEulerian();

    bool isHamiltonian();
    void makeHamiltonian();
    QPair<int, QVector<int>> solveTSP(const QString& filename) const;
    void clearHamiltonian();

    void paintEvent(QPaintEvent* event) override;

    void setShowFlowInfo(bool show);
    void setNeorFlag(bool show);
    void setShowMST(bool show);

    const QVector<QVector<int>>& getEmptyMatrix() const;
    const QVector<QVector<int>>& getAdjacencyMatrix() const;
    const QVector<QVector<int>>& getWeightMatrix() const;
    const int getVertexesCount() const;
    const QVector<QVector<int>>& getCapacityMatrix() const;
    const QVector<QVector<int>>& getCostMatrix() const;
    int getMSTWeight() const;
    const QVector<QVector<int>>& getKirchhoffMatrix() const;
    int getSpanningTreesCount() const;
    int getColorsCount() const;
    const QVector<QPair<int, int>>& getAddedEdges() const;
    const QVector<QPair<int, int>>& getHamiltonianAddedEdges() const;
    const QVector<int>& getHamiltonianAddedEdgesWeights() const;
};


void CoutMatrix(QVector<QVector<int>>& matrix);
void CoutMatrix(QVector<QVector<double>>& matrix);
void CoutMatrix(QVector<QVector<bool>>& matrix);
void CoutVector(QVector<int>& vector);
void CoutVector(QVector<QPair<int, int>>& vector);
void CoutVector(QVector<QPair<QPair<int, int>, int>>& vector);