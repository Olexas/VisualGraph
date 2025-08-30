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
    this->empty_matrix = new QVector<QVector<int>>(vertexes_count, QVector<int>(vertexes_count, 0));

    do GenerateGraph();
    while (!CheckGraph());

    FillWeightMatrix();
    FillCapacityMatrix();
    FillCostMatrix();
    GenerateSpanningTreesCount();
    Boruvka();
    GreedyColoring();
}

Graph::~Graph() {
    delete this->adjacency_matrix;
    delete this->weight_matrix;
    delete this->capacity_matrix;
    delete this->cost_matrix;
    delete this->adj_vector;
    delete this->mst_vector;
    delete this->kirchhoff_matrix;
    delete this->colors;
    delete this->empty_matrix;
}

const QVector<QVector<int>>& Graph::getEmptyMatrix() const {
    return *this->empty_matrix;
}

void Graph::GenerateGraph() {
    QVector<int>* degrees = this->GenerateDegrees(vertexes_count);

    FillAdjacencyMatrix(*degrees);

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
        QVector<double> temp_vector(n + 1);

        for (int i = 0; i < n + 1; i++) {
            temp_vector[i] = getRandNumber(lambda, a);
        }

        auto max_iter = std::max_element(temp_vector.begin(), temp_vector.end());
        double max_value = *max_iter;

        temp_vector.erase(max_iter);

        for (int i = 0; i < n; i++) {
            int k = (max_num == 0) ? n - i - 1 : max_num;
            vector[i] = round(temp_vector[i] / max_value * k);
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
    if (this->adjacency_matrix != nullptr) delete this->adjacency_matrix;;
    if (this->adj_vector != nullptr) delete this->adj_vector;

    this->adjacency_matrix = new QVector<QVector<int>>(degrees.size());
    this->adj_vector = new QVector<QPair<int, int>>;

    std::unordered_set<int> zero_row_indexes;
    std::unordered_set<int> zero_column_indexes;

    for (int i = 0; i < degrees.size(); i++) {
        (*this->adjacency_matrix)[i] = QVector<int>(degrees.size(), 0);
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
                this->adj_vector->append({ i, j });
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
                this->adj_vector->append({ index, j });
            }
            else {
                j = rand() % (index);
                (*this->adjacency_matrix)[j][index] = 1;
                this->adj_vector->append({ j, index });
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

    std::sort(this->adj_vector->begin(), this->adj_vector->end(),
        [](const QPair<int, int>& a, const QPair<int, int>& b) {
            if (a.first != b.first) {
                return a.first < b.first;
            }
            else {
                return a.second < b.second;
            }
        });
}

void Graph::FillMatrix(QVector<QVector<int>>& matrix, int max_value, double lambda, double a, bool can_negative) {
    QVector<int>* weights = new QVector<int>(this->edges_count);
    FillRandNumbersVector(*weights, max_value, lambda, a);

    int s = 0;
    for (int i = 0; i < vertexes_count; i++) {
        for (int j = 0; j < vertexes_count; j++) {
            if ((*this->adjacency_matrix)[i][j]) {
                matrix[i][j] = can_negative ? (weights->at(s) - 5 ? weights->at(s) - 5 : 1) : weights->at(s) ? weights->at(s) : 1;
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

    if (vertex1 != -1) {
        (*visited)[vertex1][vertex2] = true;
        way.push_back("(" + QString::number(vertex1) + ", " + QString::number(vertex2) + "), ");
    }

    int iter = 0;
    for (const auto& edge : *this->adj_vector) {
        if (edge.first == vertex2) {
            iter += 1;
            if (vertex1 == -1) iter += DFS(way, edge.second, vertex2, visited);
            else if (!(*visited)[vertex2][edge.second]) iter += DFS(way, edge.second, vertex2, visited);
        }
    }

    if (vertex1 == -1) {
        delete visited;
        way.chop(2);
    }
    return iter;
}

int Graph::BellmanFord(int startVertex, int endVertex, QVector<QString>& distances_result, QVector<int>& path) const {
    int n = this->vertexes_count;

    QVector<int> distances(n, INT_MAX);
    QVector<int> predecessors = QVector<int>(n, -1);

    distances[startVertex] = 0;

    int s = 0;

    for (int i = 1; i < n; ++i) {
        for (const auto& edge : *this->adj_vector) {
            s += 1;
            if (distances[edge.first] != INT_MAX) {
                int new_dist = distances[edge.first] + (*this->weight_matrix)[edge.first][edge.second];
                if (distances[edge.second] > new_dist) {
                    distances[edge.second] = new_dist;
                    predecessors[edge.second] = edge.first;
                }
            }
        }
    }

    path.clear();
    int current = endVertex;
    while (current != -1) {
        path.push_back(current);
        current = predecessors[current];
    }
    std::reverse(path.begin(), path.end());

    for (int i = 0; i < distances.size(); i++) {
        distances_result[i] = distances[i] == INT_MAX ? "--" : QString::number(distances[i]);
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
    int max_flow = 0;

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
                if (parent[next] == -1 && residualCapacity[current][next]) {
                    parent[next] = current;
                    int new_flow = std::min(flow, residualCapacity[current][next]);
                    if (next == sink) {
                        max_flow += new_flow;
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

    return max_flow;
}


std::pair<int, int> Graph::FindMinCostFlow(int source, int sink, QVector<QVector<int>>& flow, int max_flow) {
    max_flow = max_flow == -1 ? FordFulkerson(source, sink) : max_flow;
    int target_flow = (2 * max_flow) / 3;

    if (max_flow == 0) return { 0, 0 };

    QVector<QVector<int>> residualCapacity(vertexes_count, QVector<int>(vertexes_count));

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

    return { total_flow, total_cost };
}

void Graph::GenerateSpanningTreesCount() {
    int n = vertexes_count;

    this->kirchhoff_matrix = new QVector<QVector<int>>(n, QVector<int>(n, 0));
    auto& kirchhoff_matrix = *this->kirchhoff_matrix;

    for (int i = 0; i < n; ++i) {
        int degree = 0;
        for (int j = 0; j < n; ++j) {
            if ((*adjacency_matrix)[i][j] == 1) {
                kirchhoff_matrix[i][i] += 1;
                kirchhoff_matrix[j][j] += 1;
                kirchhoff_matrix[i][j] = -1;
                kirchhoff_matrix[j][i] = -1;
            }
        }
        
    }

    int size = n - 1;
    QVector<QVector<double>> minor(size, QVector<double>(size));

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            minor[i][j] = kirchhoff_matrix[i][j];
        }
    }

    this->spanning_trees_count = static_cast<int>(round(calculateDeterminant(minor)));
}

double Graph::calculateDeterminant(QVector<QVector<double>>& matrix) const {
    int n = matrix.size();
    double det = 1;

    for (int i = 0; i < n; ++i) {
        int pivot = i;
        for (int j = i + 1; j < n; ++j) {
            if (fabs(matrix[j][i]) > fabs(matrix[pivot][i])) {
                pivot = j;
            }
        }

        if (pivot != i) {
            std::swap(matrix[i], matrix[pivot]);
            det *= -1;
        }

        if (fabs(matrix[i][i]) < 1e-10) {
            return 0;
        }

        for (int j = i + 1; j < n; ++j) {
            double factor = matrix[j][i] / matrix[i][i];
            for (int k = i; k < n; ++k) {
                matrix[j][k] -= factor * matrix[i][k];
            }
        }

        det *= matrix[i][i];
    }

    return det;
}

void Graph::Boruvka(){
    if (this->mst_vector != nullptr) delete this->mst_vector;

    QVector<QPair<QPair<int, int>, int>> edges;

    for (auto edge : *this->adj_vector) {
        edges.append({ edge, (*weight_matrix)[edge.first][edge.second] });
    }

    QVector<int> component(vertexes_count);
    this->mst_vector = new QVector<QPair<int, int>>;

    for (int i = 0; i < vertexes_count; ++i) {
        component[i] = i;
    }

    int numComponents = vertexes_count;

    int weight_final = 0;
    while (numComponents > 1) {
        QVector<QPair<QPair<int, int>, int>> cheapest(vertexes_count, { { -1, INT_MAX}, -1 });

        for (const auto& edge : edges) {
            int u = edge.first.first;
            int v = edge.first.second;
            int weight = edge.second;

            int compU = component[u];
            int compV = component[v];

            if (compU != compV) {
                if (weight < cheapest[compU].first.second) {
                    cheapest[compU] = { {v, weight}, u };
                }
                if (weight < cheapest[compV].first.second) {
                    cheapest[compV] = { {u, weight}, v };
                }
            }
        }

        for (int u = 0; u < vertexes_count; ++u) {
            if (cheapest[u].first.first != -1) {
                int v = cheapest[u].first.first;
                int u_in = cheapest[u].second;

                int compU = component[u];
                int compV = component[v];

                if (compU != compV) {

                    mst_vector->append({ u_in, v });
                    weight_final += cheapest[u].first.second;

                    for (int i = 0; i < vertexes_count; ++i) {
                        if (component[i] == compV) {
                            component[i] = compU;
                        }
                    }

                    numComponents--;
                }
            }
        }
    }

    this->mst_weight = weight_final;
}

QPair<QVector<int>, QVector<int>> Graph::encodePrufer()
{
    QMap<int, QSet<int>> tree;
    QVector<int> degrees(vertexes_count, 0);

    for (const auto& edge : *this->mst_vector) {
        int u = edge.first;
        int v = edge.second;
        tree[u].insert(v);
        tree[v].insert(u);
        degrees[u]++;
        degrees[v]++;
    }

    QVector<int> pruferCode;
    QVector<int> weights;
    int n = vertexes_count;

    for (int i = 0; i < n - 1; ++i) {
        int leaf = 0;
        while (leaf < n - 1 && degrees[leaf] != 1) {
            leaf++;
        }

        int parent = *tree[leaf].begin();
        pruferCode.push_back(parent);

        int weight = 0;
        for (const auto& edge : *this->mst_vector) {
            if ((edge.first == leaf && edge.second == parent) ||
                (edge.first == parent && edge.second == leaf)) {
                weight = edge.second > edge.first ? (*weight_matrix)[edge.first][edge.second] : (*weight_matrix)[edge.second][edge.first];
                break;
            }
        }
        weights.push_back(weight);

        degrees[leaf]--;
        degrees[parent]--;
        tree[leaf].remove(parent);
        tree[parent].remove(leaf);
    }

    return qMakePair(pruferCode, weights);
}

QVector<QPair<QPair<int, int>, int>> Graph::decodePrufer(const QVector<int>& pruferCode, const QVector<int>& weights) const {
    int n = pruferCode.size() + 2;
    QVector<int> degree(n, 1);

    for (int v : pruferCode) {
        degree[v]++;
    }

    QVector<QPair<QPair<int, int>, int>> edges;

    for (int i = 0; i < pruferCode.size(); ++i) {
        int code = pruferCode[i];
        int leaf = 0;
        int weight = weights[i];
        while (leaf < n && degree[leaf] != 1) {
            leaf++;
        }

        edges.push_back({ {leaf, code}, weight });

        degree[leaf]--;
        degree[code]--;
    }

    return edges;
}

QVector<QVector<int>>* Graph::GetNeorAdjacencyMatrix() const {
    QVector<QVector<int>>* neor_adj_matrix = new QVector<QVector<int>>(this->vertexes_count, QVector<int>(this->vertexes_count));

    for (int i = 0; i < vertexes_count; i++) {
        for (int j = 0; j < vertexes_count; j++) {
            if ((*this->adjacency_matrix)[i][j]) {
                (*neor_adj_matrix)[i][j] = 1;
                (*neor_adj_matrix)[j][i] = 1;
            }
        }
    }

    return neor_adj_matrix;
}

void Graph::GreedyColoring() {
    if (this->colors != nullptr) delete this->colors;
    this->colors = new QVector<int>(vertexes_count, -1);

    QVector<QVector<int>>* neor_adj_matrix = GetNeorAdjacencyMatrix();
    QVector<int> best_coloring(vertexes_count);
    int min_colors = INT_MAX;
    
    for (int start_vertex = 0; start_vertex < vertexes_count; ++start_vertex) {
        QVector<int> current_colors(vertexes_count, -1);
        QVector<bool> available(vertexes_count, false);

        current_colors[start_vertex] = 0;

        for (int u = 1; u < vertexes_count; ++u) {
            int vertex_to_color = -1;
            for (int v = 0; v < vertexes_count; ++v) {
                if (current_colors[v] == -1) {
                    vertex_to_color = v;
                    break;
                }
            }
            if (vertex_to_color == -1) break;

            for (int v = 0; v < vertexes_count; ++v) {
                if ((*neor_adj_matrix)[vertex_to_color][v] && current_colors[v] != -1) {
                    available[current_colors[v]] = true;
                }
            }

            int cr;
            for (cr = 0; cr < vertexes_count; ++cr) {
                if (!available[cr]) break;
            }

            current_colors[vertex_to_color] = cr;

            available.fill(false);
        }

        int colors_used = *std::max_element(current_colors.begin(), current_colors.end()) + 1;

        if (colors_used < min_colors) {
            min_colors = colors_used;
            best_coloring = current_colors;
        }
    }

    *this->colors = best_coloring;
    std::unordered_set<int> uniqueSet(colors->begin(), colors->end());
    this->unique_colors = uniqueSet.size();

    delete neor_adj_matrix;
}

bool Graph::isEulerian() const {
    QVector<QVector<int>> tempAdj(vertexes_count, QVector<int>(vertexes_count, 0));

    for (int i = 0; i < vertexes_count; i++) {
        for (int j = 0; j < vertexes_count; j++) {
            if (i < j) {
            tempAdj[i][j] = (*adjacency_matrix)[i][j];
            tempAdj[j][i] = (*adjacency_matrix)[i][j];
            }
        }
    }

    for (const auto& edge : addedEdges) {
        tempAdj[edge.first][edge.second] = 1;
        tempAdj[edge.second][edge.first] = 1;
    }

    for (const auto& edge : removedEdges) {
        tempAdj[edge.first][edge.second] = 0;
        tempAdj[edge.second][edge.first] = 0;
    }

    for (int v = 0; v < vertexes_count; v++) {
        int degree = 0;
        for (int u = 0; u < vertexes_count; u++) {
            degree += tempAdj[v][u];
        }
        if (degree % 2 != 0) return false;
    }
    return true;
}

void Graph::makeEulerian() {
    if (isEulerian()) return;

    addedEdges.clear();
    removedEdges.clear();

    QVector<QVector<int>> tempAdj(vertexes_count, QVector<int>(vertexes_count, 0));
    QVector<QSet<int>> adj_list(vertexes_count);
    QVector<int> degrees(this->vertexes_count);

    for (int i = 0; i < vertexes_count; i++) {
        for (int j = 0; j < vertexes_count; j++) {
            if ((*adjacency_matrix)[i][j]) {
                degrees[i]++;
                degrees[j]++;
                adj_list[i].insert(j);
                adj_list[j].insert(i);
            }

            tempAdj[i][j] = (*adjacency_matrix)[i][j];
        }
    }

    QVector<int> res_degrees(degrees);

    for (int i = 0; i < res_degrees.size(); i++) {
        bool break_flag = false;
        if (res_degrees[i] % 2 != 0) {
            for (int j = i + 1; j < res_degrees.size(); j++) {
                if (i != j && res_degrees[j] % 2 != 0) {
                    if (tempAdj[i][j]) {
                        if (res_degrees[i] == 1) {
                            for (int vertex : adj_list[j]) {
                                if (res_degrees[vertex] % 2 == 0 && vertex != i && res_degrees[vertex] != 1) {
                                    tempAdj[j][vertex] = 0;
                                    this->removedEdges.append({ j, vertex });
                                    adj_list[j].remove(vertex);

                                    if (vertex > i) {
                                        tempAdj[i][vertex] = 1;
                                        this->addedEdges.append({ i, vertex });
                                        adj_list[i].insert(vertex);
                                    }
                                    else {
                                        tempAdj[vertex][i] = 1;
                                        this->addedEdges.append({ vertex, i });
                                        adj_list[vertex].insert(i);
                                    }

                                    break_flag = true;
                                    break;
                                }
                            }
                        }

                        else if (res_degrees[j] == 1) {
                            for (int vertex : adj_list[i]) {
                                if (res_degrees[vertex] % 2 == 0 && vertex != j && res_degrees[vertex] != 1) {
                                    tempAdj[i][vertex] = 0;
                                    this->removedEdges.append({ i, vertex });
                                    adj_list[i].remove(vertex);
                                   // res_degrees[i]--;

                                    if (vertex > j) {
                                        tempAdj[j][vertex] = 1;
                                        this->addedEdges.append({ j, vertex });
                                        adj_list[j].insert(vertex);
                                    }
                                    else {
                                        tempAdj[vertex][j] = 1;
                                        this->addedEdges.append({ vertex, j });
                                        adj_list[vertex].insert(j);
                                    }

                                    break_flag = true;
                                    break;
                                }
                            }
                        }

                        else {

                            tempAdj[i][j] = 0;
                            this->removedEdges.append({ i, j });
                            adj_list[i].remove(j);
                            res_degrees[i]--;
                            res_degrees[j]--;
                            break_flag = true;
                            break;

                        }
                    }

                    else if (!tempAdj[i][j]) {
                        tempAdj[i][j] = 1;
                        this->addedEdges.append({ i, j });
                        adj_list[i].insert(j);
                        res_degrees[i]++;
                        res_degrees[j]++;
                        break_flag = true;
                        break;
                    }
                }
            }

            if (break_flag) continue;   
        } 
    }
}

QVector<QPair<int, int>> Graph::findEulerianCycle() const {
    QStack<int> stack;
    QVector<QVector<int>> tempAdj(vertexes_count, QVector<int>(vertexes_count, 0));

    for (int i = 0; i < vertexes_count; i++) {
        for (int j = 0; j < vertexes_count; j++) {
            if (i < j) {
                tempAdj[i][j] = (*adjacency_matrix)[i][j];
                tempAdj[j][i] = (*adjacency_matrix)[i][j];
            }
        }
    }

    for (const auto& edge : addedEdges) {
        tempAdj[edge.first][edge.second] = 1;
        tempAdj[edge.second][edge.first] = 1;
    }

    for (const auto& edge : removedEdges) {
        tempAdj[edge.first][edge.second] = 0;
        tempAdj[edge.second][edge.first] = 0;
    }

    QVector<QPair<int, int>> edges;
    int tmp_vertex = -1;

    stack.push(0);
    while (!stack.empty()) {
        int v = stack.top();
        bool found_edge = false;

        for (int u = 0; u < vertexes_count; ++u) {
            if (tempAdj[v][u]) {
                stack.push(u);
                tempAdj[v][u] = 0;
                tempAdj[u][v] = 0;
                found_edge = true;
                break;
            }
        }

        if (!found_edge) {
            stack.pop();
            
            if (tmp_vertex != -1) edges.append({ tmp_vertex, v });
            tmp_vertex = v;
        }
    }

    return edges;
}

void Graph::clearEulerian() {
    this->addedEdges.clear();
    this->removedEdges.clear();
}

void Graph::clearHamiltonian() {
    this->hamiltonianAddedEdges.clear();
    this->is_hamiltonian = false;
}

bool Graph::isHamiltonian() {
    if (all_cycles != nullptr) {
        delete this->all_cycles;
        delete this->cycle_weights;
    }

    this->all_cycles = new QVector<QVector<int>>;
    this->cycle_weights = new QVector<int>;

    findAllHamiltonianCycles(*this->all_cycles, *this->cycle_weights);
    
    this->is_hamiltonian = !all_cycles->empty();

    return is_hamiltonian;

}

void Graph::makeHamiltonian() {
    if (is_hamiltonian) return;

    hamiltonianAddedEdges.clear();
    hamiltonianAddedEdgesWeights.clear();

    QVector<QVector<int>> tempAdj = *adjacency_matrix;

    QVector<QSet<int>> adj_list(vertexes_count);
    QVector<int> degrees(this->vertexes_count);

    for (int i = 0; i < vertexes_count; i++) {
        for (int j = 0; j < vertexes_count; j++) {
            if ((*adjacency_matrix)[i][j]) {
                degrees[i]++;
                degrees[j]++;

                tempAdj[j][i] = 1;
            }

            else if (i != j){
                adj_list[i].insert(j);
                adj_list[j].insert(i);
            }
        }
    }

    int condition = vertexes_count / 2;
    bool do_while_flag = false;

    do {
        int potentialEdgesCount = 0;
        for (int i = 0; i < vertexes_count; ++i) {
            if (degrees[i] >= condition) continue;
            for (int j = i + 1; j < vertexes_count; ++j) {
                if (!tempAdj[i][j] && degrees[j] < condition) {
                    potentialEdgesCount++;
                }
            }
        }

        QVector<int> newWeights(potentialEdgesCount);
        for (int i = 0; i < newWeights.size(); i++) {
            if (newWeights[i] == 0) newWeights[i] == 1;
        }

        FillRandNumbersVector(newWeights, 15, 0.5, 5);
        int weightIndex = 0;

        for (int i = 0; i < vertexes_count; ++i) {
            if (degrees[i] >= condition) continue;

            bool break_flag = false;

            //while ()

            for (int j = i + 1; j < vertexes_count; ++j) {
                if (!tempAdj[i][j] && degrees[j] < condition) {
                    tempAdj[i][j] = tempAdj[j][i] = 1;
                    degrees[i] += 1;
                    degrees[j] += 1;
                    adj_list[i].remove(j);
                    adj_list[j].remove(i);
                    hamiltonianAddedEdges.append({ i, j });
                    hamiltonianAddedEdgesWeights.append(newWeights[weightIndex++]);

                    break_flag = true;

                    if (isHamiltonian()) {
                        do_while_flag = true;
                        break;
                    }

                    break;

                }
            }

            if (do_while_flag) break;
            if (break_flag) continue;

            while (condition - degrees[i]) {
                int temp_vertex = *adj_list[i].begin();
                tempAdj[i][temp_vertex] = tempAdj[temp_vertex][i] = 1;
                degrees[i]++;
                degrees[temp_vertex]++;
                adj_list[i].remove(temp_vertex);

                if (isHamiltonian()) {
                    do_while_flag = true;
                    break;
                }
            }
        }
        condition++;
    } while (!do_while_flag);
}

//void Graph::makeHamiltonian() {
//    if (is_hamiltonian) return;
//
//    hamiltonianAddedEdges.clear();
//
//    QVector<QVector<int>> tempAdj = *adjacency_matrix;
//
//    QVector<QSet<int>> adj_list(vertexes_count);
//    QVector<int> degrees(this->vertexes_count);
//
//    for (int i = 0; i < vertexes_count; i++) {
//        for (int j = 0; j < vertexes_count; j++) {
//            if ((*adjacency_matrix)[i][j]) {
//                degrees[i]++;
//                degrees[j]++;
//
//                tempAdj[j][i] = 1;
//            }
//
//            else if (i != j) {
//                adj_list[i].insert(j);
//                adj_list[j].insert(i);
//            }
//        }
//    }
//
//    int condition = vertexes_count / 2;
//    bool do_while_flag = false;
//    do {
//        cout << 111111 << endl;
//        for (int i = 0; i < vertexes_count; ++i) {
//            if (degrees[i] >= condition) continue;
//
//            bool break_flag = false;
//            for (int j = i + 1; j < vertexes_count; ++j) {
//                if (!tempAdj[i][j] && degrees[j] < condition) {
//                    tempAdj[i][j] = tempAdj[j][i] = 1;
//                    degrees[i] += 1;
//                    degrees[j] += 1;
//                    adj_list[i].remove(j);
//                    adj_list[j].remove(i);
//                    hamiltonianAddedEdges.append({ i, j });
//
//                    break_flag = true;
//
//                    if (isHamiltonian()) do_while_flag = true;
//                }
//            }
//
//            if (do_while_flag) break;
//            if (break_flag) continue;
//
//            while (condition - degrees[i]) {
//                int temp_vertex = *adj_list[i].begin();
//                tempAdj[i][temp_vertex] = tempAdj[temp_vertex][i] = 1;
//                degrees[i]++;
//                degrees[temp_vertex]++;
//                adj_list[i].remove(temp_vertex);
//
//                if (isHamiltonian()) {
//                    do_while_flag = true;
//                    break;
//                }
//            }
//        }
//        condition++;
//    } while (!do_while_flag);
//    cout << endl;
//
//}

void Graph::findAllHamiltonianCycles(QVector<QVector<int>>& allCycles, QVector<int>& cycleWeights) const {
    QVector<int> path;
    QVector<bool> visited(vertexes_count, false);

    QVector<QVector<int>> temp_adj(vertexes_count, QVector<int>(vertexes_count));
    QVector<QVector<int>> temp_weight(vertexes_count, QVector<int>(vertexes_count));
    for (int i = 0; i < vertexes_count; ++i) {
        for (int j = 0; j < vertexes_count; ++j) {
            if (i < j) {
                temp_adj[i][j] = (*this->adjacency_matrix)[i][j];
                temp_adj[j][i] = (*this->adjacency_matrix)[i][j];
                temp_weight[i][j] = (*this->weight_matrix)[i][j];
                temp_weight[j][i] = (*this->weight_matrix)[i][j];
            }
        }
    }

    for (int k = 0; k < hamiltonianAddedEdges.size(); ++k) {
        const auto& edge = hamiltonianAddedEdges[k];
        temp_adj[edge.first][edge.second] = 1;
        temp_adj[edge.second][edge.first] = 1;

        temp_weight[edge.first][edge.second] = hamiltonianAddedEdgesWeights[k];
        temp_weight[edge.second][edge.first] = hamiltonianAddedEdgesWeights[k];
    }

    for (int start = 0; start < vertexes_count; ++start) {
        path.clear();
        visited.fill(false);
        path.append(start);
        visited[start] = true;
        hamiltonianDFS(start, start, visited, path, allCycles, cycleWeights, temp_adj, temp_weight);
    }

    removeDuplicateCycles(allCycles, cycleWeights);
}

void Graph::hamiltonianDFS(int start, int current, QVector<bool>& visited,
    QVector<int>& path, QVector<QVector<int>>& allCycles,
    QVector<int>& cycleWeights, QVector<QVector<int>>& temp_adj, QVector<QVector<int>>& temp_weight) const {

    if (path.size() == vertexes_count) {
        if (temp_adj[current][start]) {
            QVector<int> cycle = path;
            cycle.append(start);
            allCycles.append(cycle);

            int weight = 0;
            for (int i = 0; i < cycle.size() - 1; ++i) {
                weight = temp_weight[cycle[i]][cycle[i + 1]] ? weight + temp_weight[cycle[i]][cycle[i + 1]] : weight + 1;
            }
            cycleWeights.append(weight);
        }
        return;
    }

    for (int next = 0; next < vertexes_count; ++next) {
        if (temp_adj[current][next] && !visited[next]) {
            visited[next] = true;
            path.append(next);

            hamiltonianDFS(start, next, visited, path, allCycles, cycleWeights, temp_adj, temp_weight);

            path.removeLast();
            visited[next] = false;
        }
    }
}

void Graph::removeDuplicateCycles(QVector<QVector<int>>& allCycles, QVector<int>& cycleWeights) const {
    if (allCycles.isEmpty()) return;

    QVector<QVector<int>> uniqueCycles;
    QVector<int> uniqueWeights;
    QSet<uint> hashes;

    uniqueCycles.reserve(allCycles.size());
    uniqueWeights.reserve(allCycles.size());

    for (int i = 0; i < allCycles.size(); ++i) {
        const auto& cycle = allCycles[i];
        uint hash = qHash(cycleToHash(cycle));

        if (!hashes.contains(hash)) {
            hashes.insert(hash);
            uniqueCycles.append(cycle);
            uniqueWeights.append(cycleWeights[i]);
        }
    }

    allCycles = std::move(uniqueCycles);
    cycleWeights = std::move(uniqueWeights);
}

QString Graph::cycleToHash(const QVector<int>& cycle) const {
    if (cycle.size() <= 1) return "";

    int minVal = cycle[0];
    int minPos = 0;
    for (int i = 1; i < cycle.size() - 1; ++i) {
        if (cycle[i] < minVal) {
            minVal = cycle[i];
            minPos = i;
        }
    }

    QString hash;
    hash.reserve(cycle.size() * 3);

    for (int i = 0; i < cycle.size() - 1; ++i) {
        int pos = (minPos + i) % (cycle.size() - 1);
        hash += QString::number(cycle[pos]) + ",";
    }

    return hash;
}

QPair<int, QVector<int>> Graph::solveTSP(const QString& filename) const {
    QFile file(filename);
    if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        QTextStream out(&file);
        out << "Find " << this->all_cycles->size() << " cycles:\n";

        for (int i = 0; i < this->all_cycles->size(); ++i) {
            out << "Cycle " << i + 1 << ": ";
            for (int v : this->all_cycles->at(i)) {
                out << v << " ";
            }
            out << "Weight: " << this->cycle_weights->at(i) << "\n";
        }

        file.close();
    }

    if (all_cycles->empty()) {
        return qMakePair(-1, QVector<int>());
    }

    int minWeight = this->cycle_weights->at(0);
    QVector<int> minCycle = this->all_cycles->at(0);

    for (int i = 1; i < this->all_cycles->size(); ++i) {
        if (this->cycle_weights->at(i) < minWeight) {
            minWeight = this->cycle_weights->at(i);
            minCycle = this->all_cycles->at(i);
        }
    }

    return qMakePair(minWeight, minCycle);
}

void Graph::clearHighlight() {
    shortestPath.clear();
    update();
}

void Graph::paintEvent(QPaintEvent* event) {
    QPainter painter(this);
    painter.fillRect(rect(), Qt::white);
    drawGraph(painter);

    if (show_MST) {
        drawMST(painter, *this->mst_vector);
    }
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

    if (this->show_MST) {
        QVector<QColor> colorPalette = {
            Qt::red, Qt::green, Qt::blue, Qt::yellow, Qt::cyan,
            Qt::magenta, Qt::gray, Qt::darkRed, Qt::darkGreen, Qt::darkBlue
        };

        auto& colors = *this->colors;

        for (int i = 0; i < n; ++i) {
            double angle = 2 * M_PI * i / n;
            int x = centerX + (ellipseWidth / 2) * cos(angle);
            int y = centerY + (ellipseHeight / 2) * sin(angle);

            painter.setBrush(colorPalette[colors[i] % colorPalette.size()]);
            painter.drawEllipse(x - vertexRadius, y - vertexRadius, vertexRadius * 2, vertexRadius * 2);
            painter.setPen(Qt::black);
            painter.drawText(x - vertexRadius / 5.5, y + vertexRadius / 4, QString::number(i));
        }
    }

    painter.setPen(QPen(Qt::black, 1));
    painter.setBrush(Qt::black);

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if ((*adjacency_matrix)[i][j] == 1) {
                drawEdge(painter, i, j, centerX, centerY, ellipseWidth, ellipseHeight, vertexRadius, false, !this->neor_flag);
            }
        }
    }

    painter.setPen(QPen(Qt::green, 2));
    for (const auto& edge : addedEdges) {
        drawEdge(painter, edge.first, edge.second, centerX, centerY, ellipseWidth, ellipseHeight, vertexRadius, false, !this->neor_flag);
    }

    painter.setPen(QPen(Qt::red, 2));
    for (const auto& edge : removedEdges) {
        drawEdge(painter, edge.first, edge.second, centerX, centerY, ellipseWidth, ellipseHeight, vertexRadius, false, !this->neor_flag);
    }

    painter.setPen(QPen(Qt::blue, 2, Qt::DashLine));
    for (const auto& edge : hamiltonianAddedEdges) {
        drawEdge(painter, edge.first, edge.second, centerX, centerY,
            ellipseWidth, ellipseHeight, vertexRadius, false, !this->neor_flag);
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
                drawWeight(painter, i, j, centerX, centerY, ellipseWidth, ellipseHeight, vertexRadius,
                    (*weight_matrix)[i][j], false);
            }
        }
    }

    for (int k = 0; k < hamiltonianAddedEdges.size(); ++k) {
        const auto& edge = hamiltonianAddedEdges[k];
        drawWeight(painter, edge.first, edge.second, centerX, centerY,
            ellipseWidth, ellipseHeight, vertexRadius,
            hamiltonianAddedEdgesWeights[k], true);
    }
}

void Graph::drawEdge(QPainter& painter, int i, int j, int centerX, int centerY,
    int ellipseWidth, int ellipseHeight, int vertexRadius,
    bool isHighlighted, bool isForward) const {
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

    painter.drawLine(x1, y1, x2, y2);

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
    int ellipseWidth, int ellipseHeight, int vertexRadius,
    int weight, bool isAddedEdge) const {

    bool isRemoved = removedEdges.contains({ i,j }) || removedEdges.contains({ j,i });
    if (isRemoved) return;

    if (weight == 0) return;

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
    double posRatio = 0.4;
    int textX = x1 + (x2 - x1) * posRatio;
    int textY = y1 + (y2 - y1) * posRatio;

    int value;

    if (showFlowInfo) {
        value = (*capacity_matrix)[i][j];
        painter.setBrush(QColor(255, 255, 200, 200));
    }

    else {
        value = weight;
        painter.setBrush(QColor(173, 216, 230, 200));
    }

    painter.save();
    int s = showFlowInfo ? 2 : 1;
    while (s) {
        painter.translate(textX, textY);

        painter.setPen(Qt::NoPen);
        painter.drawRoundedRect(-rectWidth / 2, -rectHeight / 2, rectWidth, rectHeight, 5, 5);

        painter.setPen(Qt::black);
        painter.setFont(QFont("Arial", 8, QFont::Bold));
        painter.drawText(QRect(-rectWidth / 2, -rectHeight / 2, rectWidth, rectHeight),
            Qt::AlignCenter, QString::number(value));

        painter.restore();

        s--;
        
        if (showFlowInfo && s == 1) {
            textX += 25;
            value = (*cost_matrix)[i][j];

            painter.save();
            painter.setBrush(QColor(144, 238, 144));
        }
    }
}

void Graph::drawMST(QPainter& painter, const QVector<QPair<int, int>>& mstEdges) const{
    int vertexRadius = 20;
    int availableWidth = width();
    int availableHeight = height();
    int ellipseWidth = availableWidth - 2 * vertexRadius - 40;
    int ellipseHeight = availableHeight - 2 * vertexRadius - 40;
    int centerX = availableWidth / 2;
    int centerY = availableHeight / 2;

    painter.setPen(QPen(Qt::green, 3));
    painter.setBrush(Qt::green);

    for (const auto& edge : mstEdges) {
        int i = edge.first;
        int j = edge.second;

        if (i > j) std::swap(i, j);

        drawEdge(painter, i, j, centerX, centerY, ellipseWidth, ellipseHeight, vertexRadius, false, !this->neor_flag);
        drawWeight(painter, i, j, centerX, centerY, ellipseWidth, ellipseHeight, vertexRadius,
            (*weight_matrix)[i][j], false);
    }
}

void Graph::setShowFlowInfo(bool show) {
    this->showFlowInfo = show;
}

void Graph::setNeorFlag(bool show) {
    this->neor_flag = show;
}

void Graph::setShowMST(bool show) {
    this->show_MST = show;
}

int Graph::getMSTWeight() const {
    return this->mst_weight;
}

const QVector<QVector<int>>& Graph::getKirchhoffMatrix() const {
    return *this->kirchhoff_matrix;
}

int Graph::getSpanningTreesCount() const {
    return this->spanning_trees_count;
}

int Graph::getColorsCount() const {
    return this->unique_colors;
}

const QVector<QPair<int, int>>& Graph::getAddedEdges() const {
    return addedEdges;
}

const QVector<QPair<int, int>>& Graph::getHamiltonianAddedEdges() const {
    return hamiltonianAddedEdges;
}

const QVector<int>& Graph::getHamiltonianAddedEdgesWeights() const {
    return hamiltonianAddedEdgesWeights;
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

void CoutMatrix(QVector<QVector<double>>& matrix) {
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

void CoutVector(QVector<int>& vector) {
    for (int i = 0; i < vector.size(); i++) {
        cout << vector[i] << " ";
    }
    cout << endl;
}

void CoutVector(QVector<QPair<int, int>>& vector) {
    for (int i = 0; i < vector.size(); i++) {
        cout << vector[i].first << " " << vector[i].second << endl;
    }
    cout << endl;
}

void CoutVector(QVector<QPair<QPair<int, int>, int>>& vector) {
    for (int i = 0; i < vector.size(); i++) {
        cout << vector[i].first.first << " " << vector[i].first.second << " " << vector[i].second << endl;
    }
    cout << endl;
}
