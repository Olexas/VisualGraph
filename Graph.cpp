#include <QPainter>
#include <QWidget>
#include <vector>
#include <QVector>
#include <QList>
#include <qalgorithms.h>
#include <unordered_set>
#include <cmath>

#include "Graph.h"

#include <iostream>

using std::cout;
using std::endl;

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

Graph::Graph(int vertexes_count, QWidget* parent) : QMainWindow(parent) {
    this->vertexes_count = vertexes_count;

    do GenerateGraph();
    while (!CheckGraph());
}

void Graph::GenerateGraph() {
    QVector<int>* degrees = this->GenerateDegrees(vertexes_count);

    FillAdjacencyMatrix(*degrees);
    FillWeightMatrix();
    FillCapacityMatrix();
    FillCostMatrix();

    delete degrees;
}

bool Graph::CheckGraph() const{
    QVector<int> component(this->vertexes_count);

    int num = 0;
    for (int v = 0; v < this->vertexes_count; v++) {
        if (!component[v]) dfs(v, ++num, component);

        if (num > 1) {
            return false;
        }
    }

    return true;
}

void Graph::dfs(int v, int num, QVector<int>& component) const {
    component[v] = num;
    for (int u = 0; u < this->vertexes_count; u++) {
        if (((*this->adjacency_matrix)[v][u] != 0 || (*this->adjacency_matrix)[u][v] != 0) && !component[u]) dfs(u, num, component);
    }
}

QVector<int>* Graph::GenerateDegrees(int vertexes_count) {
    srand(time(0));

    QVector<int>* degrees = new QVector<int>(vertexes_count);
    FillRandNumbersVector(*degrees);

    return degrees;
}

void Graph::FillRandNumbersVector(QVector<int>& vector, int max_num, double lambda, double a) const {
    int n = vector.size();
    QVector<double> temp_vector(n);

    for (int i = 0; i < n; i++) {
        temp_vector[i] = getRandNumber(lambda, a);
    }

    auto max_element = std::max_element(temp_vector.begin(), temp_vector.end());

    for (int i = 0; i < n; i++) {
        int k = !max_num ? n - i - 1 : max_num;
        vector[i] = round(temp_vector[i] / *max_element * k);
    }
}

double Graph::getRandNumber(double lambda, double a) const {
    double b = exp(1) / (exp(1) + a);

    double v;

    while (true) {
        double r1 = double(rand() % 100) / 100;
        double r2 = double(rand() % 100) / 100;

        if (r1 < b) {
            v = pow(r1 / b, 1 / a);
            if (r2 > exp(-v)) continue;
        }

        else {
            v = 1 - log((1 - r1) / (1 - b));
            if (r2 > pow(v, a - 1)) continue;
        }

        return v / lambda;
    }
}

void Graph::FillAdjacencyMatrix(QVector<int>& degrees) {
    if (this->adjacency_matrix != nullptr) delete this->adjacency_matrix;
    if (this->adjLists != nullptr) delete [] this->adjLists;

    this->adjacency_matrix = new QVector<QVector<int>>(degrees.size());
    this->adjLists = new QList<int>[vertexes_count];

    std::unordered_set<int> zero_row_indexes;
    std::unordered_set<int> zero_column_indexes;

    for (int i = 0; i < degrees.size(); i++) {
        for (int j = 0; j < degrees.size(); j++) {
            (*this->adjacency_matrix)[i].push_back(0);
        }
        zero_column_indexes.insert(i);
    }
    
    zero_row_indexes.insert(degrees.size() - 1);

    for (int i = 0; i < degrees.size() - 1; i++) {
        int k = degrees[i]; int j;

        if (k == 0) zero_row_indexes.insert(i);

        while (k != 0) {
            j = rand() % (degrees.size() - i - 1) + i + 1;
            if ((*this->adjacency_matrix)[i][j] == 0) {
                (*this->adjacency_matrix)[i][j] = 1;
                adjLists[i].push_front(j);
                this->edges_count += 1;
                zero_column_indexes.erase(j);
                k--;
            }
        }
    }

    int j;
    for (const auto& index : zero_row_indexes) {
        if (zero_column_indexes.count(index)) {
            if (index != degrees.size() - 1) {
                j = rand() % (degrees.size() - index - 1) + index + 1;
                (*this->adjacency_matrix)[index][j] = 1;
                adjLists[index].push_front(j);
            }
            else {
                j = rand() % (index);
                (*this->adjacency_matrix)[j][index] = 1;
                adjLists[j].push_front(index);
            }
            this->edges_count += 1;

            //if ((rand() % 2 || index == 0) && index != degrees.size() - 1) {
            /*}
            else {
                j = rand() % (index);
                (*this->adjacency_matrix)[j][index] = 1;
            }*/
        }
    }

    for (int i = 0; i < this->vertexes_count; i++) {
        qSort(this->adjLists[i]);
    }
}

void Graph::FillMatrix(QVector<QVector<int>>& matrix, int max_value, double lambda, double a, bool can_negative) {
    QVector<int>* weights = new QVector<int>(this->edges_count);
    FillRandNumbersVector(*weights, max_value, lambda, a);

    int s = 0;
    for (int i = 0; i < vertexes_count; i++) {
        for (int j = 0; j < vertexes_count; j++) {
            if ((*this->adjacency_matrix)[i][j]) {
                matrix[i][j] = can_negative ? (weights->at(s) - 7 ? weights->at(s) - 7 : 1) : weights->at(s) ? weights->at(s) : 1;
                s++;
            }
            else matrix[i][j] = 0;
        }
    }

    delete weights;
}

void Graph::FillWeightMatrix() {
    if (this->weight_matrix != nullptr) delete this->weight_matrix;

    this->weight_matrix = new QVector<QVector<int>>(this->vertexes_count, QVector<int>(this->vertexes_count));
    FillMatrix(*weight_matrix, 15, 0.1, 10, true);
}

void Graph::FillCapacityMatrix() {
    if (this->capacity_matrix != nullptr) delete this->capacity_matrix;

    this->capacity_matrix = new QVector<QVector<int>>(this->vertexes_count, QVector<int>(this->vertexes_count));
    FillMatrix(*capacity_matrix, 15, 0.1, 10);
}

void Graph::FillCostMatrix() {
    if (this->cost_matrix != nullptr) delete this->cost_matrix;

    this->cost_matrix = new QVector<QVector<int>>(this->vertexes_count, QVector<int>(this->vertexes_count));
    FillMatrix(*cost_matrix, 15, 0.1, 10);
}

QVector<QVector<int>>* Graph::getMinWayByShimbell(int edges_count) const{
    return this->getWayByShimbell(edges_count);
}

QVector<QVector<int>>* Graph::getMaxWayByShimbell(int edges_count) const{
    return this->getWayByShimbell(edges_count, false);
}

QVector<QVector<int>>* Graph::getWayByShimbell(int edges_count, bool is_min) const{
    if (edges_count == 1) {
        QVector<QVector<int>>* weight_matrix_copy = new QVector<QVector<int>>(*this->weight_matrix);
        return weight_matrix_copy;
    }

    QVector<QVector<int>>* matrix = new QVector<QVector<int>>(this->vertexes_count);
    QVector<QVector<int>>* copy_matrix = new QVector<QVector<int>>(*this->weight_matrix);

    for (int i = 0; i < this->vertexes_count; i++) {
        (*matrix)[i] = QVector<int>(this->vertexes_count);
    }

    for (int edge = 0; edge < edges_count - 1; edge++) {
        for (int i = 0; i < this->vertexes_count; i++) {
            for (int j = 0; j < this->vertexes_count; j++) {
                if (j > i) {
                    for (int l = 0; l < this->vertexes_count; l++) {
                        if ((*copy_matrix)[i][l] && (*this->weight_matrix)[l][j]) {
                            if (is_min) (*matrix)[i][j] = (*matrix)[i][j] ? fmin((*copy_matrix)[i][l] + (*this->weight_matrix)[l][j], (*matrix)[i][j]) : (*copy_matrix)[i][l] + (*this->weight_matrix)[l][j];
                            else (*matrix)[i][j] = (*matrix)[i][j] ? fmax((*copy_matrix)[i][l] + (*this->weight_matrix)[l][j], (*matrix)[i][j]) : (*copy_matrix)[i][l] + (*this->weight_matrix)[l][j];
                        }
                    }
                }
            }
        }
        if (edge != edges_count - 2) {
            delete copy_matrix;
            copy_matrix = matrix;
            matrix = new QVector<QVector<int>>(this->vertexes_count);
            for (int i = 0; i < this->vertexes_count; i++) {
                (*matrix)[i] = QVector<int>(this->vertexes_count);
            }
        }
    }

    delete copy_matrix;
    return matrix;
}

int Graph::FindPathCount(int first_vertex, int second_vertex) const {
    QVector<QVector<int>>* temp_matrix = this->adjacency_matrix;
    QVector<QVector<int>>* matrix = new QVector<QVector<int>>(this->vertexes_count, QVector<int>(this->vertexes_count));

    for (int i = 0; i < this->vertexes_count; i++) {
        for (int j = 0; j < this->vertexes_count; j++) {
            for (int k = 0; k < this->vertexes_count; k++) {
                (*matrix)[j][k] = (*temp_matrix)[j][k] + (*temp_matrix)[j][i] * (*temp_matrix)[i][k];
            }
        }
        temp_matrix = matrix;
    }

    return (*matrix)[first_vertex][second_vertex];
}

int Graph::DFS(QString& way, int vertex2, int vertex1, QVector<QVector<bool>>* visited) const{
    if (visited == nullptr) visited = new QVector<QVector<bool>>(this->vertexes_count, QVector<bool>(this->vertexes_count));
    
    QList<int> adjList = this->adjLists[vertex2];

    if (vertex1 != -1) {
        (*visited)[vertex1][vertex2] = true;
        way.push_back("(" + QString::number(vertex1) + ", " + QString::number(vertex2) + "), ");
    }

    // TODO: спросить про итерации
    int iter = 0;
    QList<int>::iterator i;
    for (i = adjList.begin(); i != adjList.end(); ++i) {
        iter += 1;
        if (vertex1 == -1) iter += DFS(way, *i, vertex2, visited);
        else if (!(*visited)[vertex2][*i]) iter += DFS(way, *i, vertex2, visited);
    }

    if (vertex1 == -1) {
        delete visited;
        way.chop(2);
    }
    return iter;
}

int Graph::BellmanFord(int startVertex, int endVertex, QVector<int>& distances, QVector<int>& path) const {
    int n = this->vertexes_count;

    distances = QVector<int>(n, INT_MAX);
    QVector<int> predecessors = QVector<int>(n, -1);

    distances[startVertex] = 0;

    int s = 0;

    for (int i = 0; i < n - 1; ++i) {
        for (int u = 0; u < n; ++u) {
            for (int v = 0; v < n; ++v) {
                s += 1;
                if ((*this->weight_matrix)[u][v] != 0 && distances[u] != INT_MAX) {
                    int newDist = distances[u] + (*this->weight_matrix)[u][v];
                    if (newDist < distances[v]) {
                        distances[v] = newDist;
                        predecessors[v] = u;
                    }
                }
            }
        }
    }

    //for (int u = 0; u < n; ++u) {
    //    for (int v = 0; v < n; ++v) {
    //        if ((*this->weight_matrix)[u][v] != 0 && distances[u] != INT_MAX) {
    //            if (distances[u] + (*this->weight_matrix)[u][v] < distances[v]) {
    //                const_cast<Graph*>(this)->shortestPath = QVector<int>(this->vertexes_count);
    //                return -1;
    //            }
    //        }
    //    }
    //}

    path.clear();
    int current = endVertex;
    while (current != -1) {
        path.push_back(current);
        current = predecessors[current];
    }
    std::reverse(path.begin(), path.end());

    for (int i = 0; i < distances.size(); i++) {
        if (distances[i] == INT_MAX) distances[i] = 0;
    }
    const_cast<Graph*>(this)->shortestPath = path;

    return s;
}

int Graph::FordFulkerson(int source, int sink) {
    if (source == sink) {
        return 0;
    }

    QVector<QVector<int>> residualCapacity(vertexes_count, QVector<int>(vertexes_count));
    for (int i = 0; i < vertexes_count; ++i) {
        for (int j = 0; j < vertexes_count; ++j) {
            residualCapacity[i][j] = (*capacity_matrix)[i][j];
        }
    }

    QVector<int> parent(vertexes_count);
    int maxFlow = 0;

    while (true) {
        std::fill(parent.begin(), parent.end(), -1);
        parent[source] = -2;
        std::queue<std::pair<int, int>> q;
        q.push({ source, INT_MAX });

        bool foundPath = false;
        while (!q.empty() && !foundPath) {
            int current = q.front().first;
            int flow = q.front().second;
            q.pop();

            for (int next = 0; next < vertexes_count; ++next) {
                if (parent[next] == -1 && residualCapacity[current][next] > 0) {
                    parent[next] = current;
                    int new_flow = std::min(flow, residualCapacity[current][next]);
                    if (next == sink) {
                        maxFlow += new_flow;
                        int v = next;
                        while (v != source) {
                            int u = parent[v];
                            residualCapacity[u][v] -= new_flow;
                            residualCapacity[v][u] += new_flow;
                            v = u;
                        }
                        foundPath = true;
                        break;
                    }
                    q.push({ next, new_flow });
                }
            }
        }

        if (!foundPath) {
            break;
        }
    }

    cout << "№1" << endl;;
    CoutMatrix(residualCapacity);
    cout << endl;

    return maxFlow;
}

std::pair<int, int> Graph::FindMinCostFlow(int source, int sink, int max_flow) {
    max_flow = max_flow == -1 ? FordFulkerson(source, sink) : max_flow;
    int target_flow = (2 * max_flow) / 3;

    if (max_flow == 0) return { 0, 0 };

    QVector<QVector<int>> residualCapacity(vertexes_count, QVector<int>(vertexes_count));
    QVector<QVector<int>> flow(vertexes_count, QVector<int>(vertexes_count, 0));

    for (int i = 0; i < vertexes_count; ++i) {
        for (int j = 0; j < vertexes_count; ++j) {
            residualCapacity[i][j] = (*capacity_matrix)[i][j];
        }
    }

    int total_flow = 0;
    int total_cost = 0;

    while (total_flow < target_flow) {
        QVector<int> distances(vertexes_count, INT_MAX);
        QVector<int> predecessors(vertexes_count, -1);
        distances[source] = 0;

        for (int k = 0; k < vertexes_count - 1; ++k) {
            for (int i = 0; i < vertexes_count; ++i) {
                for (int j = 0; j < vertexes_count; ++j) {
                    if (residualCapacity[i][j] > 0 && distances[i] != INT_MAX &&
                        distances[j] > distances[i] + (*cost_matrix)[i][j]) {
                        distances[j] = distances[i] + (*cost_matrix)[i][j];
                        predecessors[j] = i;
                    }
                }
            }
        }

        for (int i = 0; i < vertexes_count; ++i) {
            for (int j = 0; j < vertexes_count; ++j) {
                if (residualCapacity[i][j] > 0 && distances[i] != INT_MAX &&
                    distances[j] > distances[i] + (*cost_matrix)[i][j]) {
                }
            }
        }

        if (predecessors[sink] == -1) {
            break;
        }

        int path_flow = INT_MAX;
        for (int v = sink; v != source; v = predecessors[v]) {
            int u = predecessors[v];
            path_flow = std::min(path_flow, residualCapacity[u][v]);
        }

        path_flow = std::min(path_flow, target_flow - total_flow);

        for (int v = sink; v != source; v = predecessors[v]) {
            int u = predecessors[v];
            residualCapacity[u][v] -= path_flow;
            residualCapacity[v][u] += path_flow;
            flow[u][v] += path_flow;
            total_cost += path_flow * (*cost_matrix)[u][v];
        }

        total_flow += path_flow;
    }

    cout << "№2" << endl;;
    CoutMatrix(flow);
    cout << endl;

    return { total_flow, total_cost };
}

void Graph::clearHighlight() {
    shortestPath.clear();
    update(); // Перерисовываем граф
}

void Graph::paintEvent(QPaintEvent* event) {
    QPainter painter(this);
    painter.fillRect(rect(), Qt::white);
    drawGraph(painter);
}

void Graph::drawGraph(QPainter& painter) {
    int n = this->vertexes_count;
    int vertexRadius = 20;

    painter.setRenderHint(QPainter::Antialiasing);

    int availableWidth = width();
    int availableHeight = height();

    int ellipseWidth = availableWidth - 2 * vertexRadius - 40;
    int ellipseHeight = availableHeight - 2 * vertexRadius - 40;

    int centerX = availableWidth / 2;
    int centerY = availableHeight / 2;

    painter.setRenderHint(QPainter::Antialiasing);

    for (int i = 0; i < n; ++i) {
        double angle = 2 * M_PI * i / n;
        int x = centerX + (ellipseWidth / 2) * cos(angle);
        int y = centerY + (ellipseHeight / 2) * sin(angle);

        painter.setBrush(Qt::lightGray);
        painter.drawEllipse(x - vertexRadius, y - vertexRadius, vertexRadius * 2, vertexRadius * 2);
        painter.setPen(Qt::black);
        painter.drawText(x - vertexRadius / 5.5, y + vertexRadius / 4, QString::number(i));
    }

    painter.setPen(QPen(Qt::black, 1));
    painter.setBrush(Qt::black);

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if ((*adjacency_matrix)[i][j] == 1) {
                drawEdge(painter, i, j, centerX, centerY, ellipseWidth, ellipseHeight, vertexRadius, false, true);
            }
        }
    }

    if (!shortestPath.isEmpty()) {
        painter.setPen(QPen(Qt::red, 3));
        painter.setBrush(Qt::red);

        for (int k = 0; k < shortestPath.size() - 1; ++k) {
            int i = shortestPath[k];
            int j = shortestPath[k + 1];

            if ((*adjacency_matrix)[i][j] == 1 || (*adjacency_matrix)[j][i] == 1) {
                bool isForward = (*adjacency_matrix)[i][j] == 1;
                drawEdge(painter, i, j, centerX, centerY, ellipseWidth, ellipseHeight, vertexRadius, true, isForward);
            }
        }
    }

    painter.setPen(Qt::black);
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if ((*adjacency_matrix)[i][j] == 1) {
                drawWeight(painter, i, j, centerX, centerY, ellipseWidth, ellipseHeight, vertexRadius);
            }
        }
    }
}

void Graph::drawEdge(QPainter& painter, int i, int j, int centerX, int centerY,
    int ellipseWidth, int ellipseHeight, int vertexRadius,
    bool isHighlighted, bool isForward) const {
    int n = this->adjacency_matrix->size();
    // Координаты вершин
    double angle1 = 2 * M_PI * i / n;
    int x1c = centerX + (ellipseWidth / 2) * cos(angle1);
    int y1c = centerY + (ellipseHeight / 2) * sin(angle1);

    double angle2 = 2 * M_PI * j / n;
    int x2c = centerX + (ellipseWidth / 2) * cos(angle2);
    int y2c = centerY + (ellipseHeight / 2) * sin(angle2);

    // Вычисляем точки на границе вершин
    double edgeAngle = atan2(y2c - y1c, x2c - x1c);
    double x1 = x1c + vertexRadius * cos(edgeAngle);
    double y1 = y1c + vertexRadius * sin(edgeAngle);
    double x2 = x2c - vertexRadius * cos(edgeAngle);
    double y2 = y2c - vertexRadius * sin(edgeAngle);

    // Рисуем линию
    painter.drawLine(x1, y1, x2, y2);

    // Рисуем стрелку (только для ориентированных рёбер)
    if (isHighlighted || isForward) {
        double arrowLength = isHighlighted ? 12 : 10;
        double arrowAngle = M_PI / 6;

        QPointF tip(isForward ? x2 : x1, isForward ? y2 : y1);
        QPointF arrowPoint(isForward ? (x2 - vertexRadius * cos(edgeAngle))
            : (x1 + vertexRadius * cos(edgeAngle)),
            isForward ? (y2 - vertexRadius * sin(edgeAngle))
            : (y1 + vertexRadius * sin(edgeAngle)));

        painter.drawLine(tip, arrowPoint);

        QPointF arrow1(arrowPoint.x() - arrowLength * cos(edgeAngle - arrowAngle),
            arrowPoint.y() - arrowLength * sin(edgeAngle - arrowAngle));
        QPointF arrow2(arrowPoint.x() - arrowLength * cos(edgeAngle + arrowAngle),
            arrowPoint.y() - arrowLength * sin(edgeAngle + arrowAngle));

        painter.drawPolygon(QPolygonF() << tip << arrow1 << arrow2);
    }
}

void Graph::drawWeight(QPainter& painter, int i, int j, int centerX, int centerY,
    int ellipseWidth, int ellipseHeight, int vertexRadius) const {
    int n = this->adjacency_matrix->size();

    double angle1 = 2 * M_PI * i / n;
    int x1c = centerX + (ellipseWidth / 2) * cos(angle1);
    int y1c = centerY + (ellipseHeight / 2) * sin(angle1);

    double angle2 = 2 * M_PI * j / n;
    int x2c = centerX + (ellipseWidth / 2) * cos(angle2);
    int y2c = centerY + (ellipseHeight / 2) * sin(angle2);

    double edgeAngle = atan2(y2c - y1c, x2c - x1c);
    double x1 = x1c + vertexRadius * cos(edgeAngle);
    double y1 = y1c + vertexRadius * sin(edgeAngle);
    double x2 = x2c - vertexRadius * cos(edgeAngle);
    double y2 = y2c - vertexRadius * sin(edgeAngle);

    int rectWidth = 30;
    int rectHeight = 20;

    if (showFlowInfo) {
        int capacity = (*capacity_matrix)[i][j];
        if (capacity != 0) {
            double posRatio = 0.4;
            int textX = x1 + (x2 - x1) * posRatio;
            int textY = y1 + (y2 - y1) * posRatio;

            // Рисуем прямоугольник для пропускной способности
            painter.save();
            painter.translate(textX, textY);
            painter.setBrush(QColor(255, 255, 200, 200)); // Цвет для capacity_matrix
            painter.setPen(Qt::NoPen);
            painter.drawRoundedRect(-rectWidth / 2, -rectHeight / 2, rectWidth, rectHeight, 5, 5);

            painter.setPen(Qt::black);
            painter.setFont(QFont("Arial", 8, QFont::Bold));
            painter.drawText(QRect(-rectWidth / 2, -rectHeight / 2, rectWidth, rectHeight),
                Qt::AlignCenter, QString::number(capacity));

            painter.restore();
        }

        // Отображаем стоимость
        int cost = (*cost_matrix)[i][j];
        if (cost != 0) {
            double posRatio = 0.4; // Увеличиваем значение для рисования ниже
            int textX = x1 + (x2 - x1) * posRatio + 25;
            int textY = y1 + (y2 - y1) * posRatio; // Смещаем вниз для визуального отделения

            // Рисуем прямоугольник для стоимости
            painter.save();
            painter.translate(textX, textY);
            painter.setBrush(QColor(144, 238, 144)); // Цвет для cost_matrix
            painter.setPen(Qt::NoPen);
            painter.drawRoundedRect(-rectWidth / 2, -rectHeight / 2, rectWidth, rectHeight, 5, 5);

            painter.setPen(Qt::black);
            painter.setFont(QFont("Arial", 8, QFont::Bold));
            painter.drawText(QRect(-rectWidth / 2, -rectHeight / 2, rectWidth, rectHeight),
                Qt::AlignCenter, QString::number(cost));

            painter.restore();
        }
    }

    else {
        int weight = (*weight_matrix)[i][j];
        if (weight != 0) {
            double posRatio = 0.4;
            int textX = x1 + (x2 - x1) * posRatio;
            int textY = y1 + (y2 - y1) * posRatio;

            // Поворачиваем систему координат для выравнивания текста вдоль ребра
            painter.save();
            painter.translate(textX, textY);
            //painter.rotate(edgeAngle * 180 / M_PI);

            painter.setBrush(QColor(173, 216, 230, 200));
            painter.setPen(Qt::NoPen);
            painter.drawRoundedRect(-rectWidth / 2, -rectHeight / 2, rectWidth, rectHeight, 5, 5);

            painter.setPen(Qt::black);
            painter.setFont(QFont("Arial", 8, QFont::Bold));
            painter.drawText(QRect(-rectWidth / 2, -rectHeight / 2, rectWidth, rectHeight),
                Qt::AlignCenter, QString::number(weight));

            painter.restore();
        }
    }
}

void Graph::setShowFlowInfo(bool show) {
    this->showFlowInfo = show;
}

const QVector<QVector<int>>& Graph::getAdjacencyMatrix() const {
    return *this->adjacency_matrix;
}

const QVector<QVector<int>>& Graph::getWeightMatrix() const {
    return *this->weight_matrix;
}

const QVector<QVector<int>>& Graph::getCapacityMatrix() const {
    return *this->capacity_matrix;
}

const QVector<QVector<int>>& Graph::getCostMatrix() const {
    return *this->cost_matrix;
}

const int Graph::getVertexesCount() const {
    return this->vertexes_count;
}

void CoutMatrix(QVector<QVector<int>>& matrix) {
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix.size(); j++) {
            cout << matrix[i][j] << "\t";
        }
        cout << endl;
    }
}

void CoutMatrix(QVector<QVector<bool>>& matrix) {
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix.size(); j++) {
            if (matrix[i][j]) cout << 1 << " ";
            else cout << 0 << " ";
        }
        cout << endl;
    }
}
