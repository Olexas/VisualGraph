#include "MainWindow.h"

MainWindow::MainWindow() {
    resize(1920, 1080);
    setupUI();
}

void MainWindow::setupUI() {
    QWidget* centralWidget = new QWidget(this);
    QHBoxLayout* mainLayout = new QHBoxLayout(centralWidget);

    createInputWidget();
    initializeTables();
    initializeLabels();

    graph_layout = new QVBoxLayout();
    graph_layout->addWidget(inputWidget);
    graph_layout->addWidget(graph);

    tab_widget = new QTabWidget(this);
    setupTabWidget();

    updateDynamicLabels();

    mainLayout->addLayout(graph_layout);
    mainLayout->addWidget(tab_widget);

    setCentralWidget(centralWidget);
    setupConnections();
}

void MainWindow::createInputWidget() {
    inputWidget = new QWidget(this);
    QHBoxLayout* inputLayout = new QHBoxLayout(inputWidget);

    auto createInputField = [&](QLineEdit*& field, const QString& placeholder, bool visible = true) {
        field = new QLineEdit(inputWidget);
        field->setPlaceholderText(placeholder);
        field->setStyleSheet("color: grey;");
        field->setVisible(visible);
        inputLayout->addWidget(field);

        QRegularExpression vertexRegex("^[0-9]*$");
        QRegularExpressionValidator* validator = new QRegularExpressionValidator(vertexRegex, this);
        field->setValidator(validator);
        };

    createInputField(vertex_input, u8"Введите количество вершин (от 2 до 100)");
    createInputField(shimbell_edges_input, "", false);
    createInputField(warshall_first_input, u8"Введите номер первой вершины", false);
    createInputField(warshall_second_input, u8"Введите номер второй вершины", false);

    auto createButton = [&](QPushButton*& button, const QString& text, bool visible = true) {
        button = new QPushButton(text, inputWidget);
        button->setVisible(visible);
        inputLayout->addWidget(button);
        };

    createButton(generate_button, u8"Сгенерировать");
    createButton(calculate_button, u8"Рассчитать", false);
    createButton(warshall_calculate_button, u8"Рассчитать", false);
}

void MainWindow::initializeTables() {
    graph = new Graph(10);
    int vertexes_count = graph->getVertexesCount();

    tables = {
        {"adjacency",
            {new TableWidget(graph->getAdjacencyMatrix(), this, double_tables_k),
             new QLabel(u8"Матрица смежности", this),
             &Graph::getAdjacencyMatrix}},

        {"weight",
            {new TableWidget(graph->getWeightMatrix(), this, double_tables_k),
             new QLabel(u8"Весовая матрица", this),
             &Graph::getWeightMatrix}},

        {"min_path",
            {new TableWidget(graph->getEmptyMatrix(), this, double_tables_k),
             new QLabel(u8"Минимальный путь", this),
             &Graph::getEmptyMatrix}},

        {"max_path",
            {new TableWidget(graph->getEmptyMatrix(), this, double_tables_k),
             new QLabel(u8"Максимальный путь", this),
             &Graph::getEmptyMatrix}},

        {"flow",
            {new TableWidget(graph->getEmptyMatrix(), this, double_tables_k, QColor(255, 255, 200)),
             new QLabel(u8"Матрица распределения потока", this),
             &Graph::getEmptyMatrix}},

        {"cost",
            {new TableWidget(graph->getCostMatrix(), this, double_tables_k, QColor(144, 238, 144)),
             new QLabel(u8"Матрица стоимостей", this),
             &Graph::getCostMatrix}},

        {"kirchhoff",
            {new TableWidget(graph->getKirchhoffMatrix(), this, double_tables_k, QColor(200, 200, 255)),
             new QLabel(u8"Матрица Кирхгофа", this),
             &Graph::getKirchhoffMatrix}}
    };

    for (auto& table : tables) {
        table.widget->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    }
}

void MainWindow::initializeLabels() {
    labels = {
        {"warshall", new QLabel(u8"1) Алгоритм Уоршалла \n Количество маршрутов из вершины ... в вершину ... : ", this)},
        {"dfs", new QLabel(u8"2) Обход графа по ребрам в глубину из вершины ... : ", this)},
        {"bellman_ford", new QLabel(u8"3) Алгоритм Беллмана-Форда \n Кратчайшее расстояние из вершины ... в вершину ... : ", this)},
        {"iterations", new QLabel(u8"Количество итераций алгоритма 2: ...\nКоличество итераций алгоритма 3: ...", this)},
        {"spanning_tree_count", new QLabel(this)},
        {"min_spanning_tree", new QLabel(this)},
        {"colors", new QLabel(this)},
        {"prufer_code", new QLabel(this)},
        {"flow", new QLabel(u8"Максимальный поток по алгоритму Форда-Фалкерсона : ...", this)},
        {"min_cost_flow", new QLabel(u8"Минимальная стоимость величины потока ... : ...", this)},
        {"eulerian_info", new QLabel(this)}
    };

    for (auto& label : labels) {
        label->setStyleSheet("font-size: 12pt;");
        label->setAlignment(Qt::AlignHCenter);
        label->setWordWrap(true);
    }
}

void MainWindow::updateDynamicLabels() {
    labels["spanning_tree_count"]->setText(u8"1) Число остовных деревьев графа : " +
        QString::number(graph->getSpanningTreesCount()));

    labels["min_spanning_tree"]->setText(u8"2) Минимальный по весу остов : " +
        QString::number(graph->getMSTWeight()));

    labels["colors"]->setText(u8"3) Минимальное количество цветов : " +
        QString::number(graph->getColorsCount()));

    labels["prufer_code"]->setText(updatePruferInfo());

    QString eulerian_text = this->graph->getVertexesCount() == 2 ? u8"Эйлеров цикл невозможен!\n" : this->graph->isEulerian() ? u8"1) Граф уже является эйлеровым\n" : u8"1) Граф не является эйлеровым\n";
    eulerian_text += this->graph->getVertexesCount() == 2 ? u8"Гамильтонов цикл невозможен!" : this->graph->isHamiltonian() ? u8"2) Граф уже является гамильтоновым" : u8"2) Граф не является гамильтоновым";
    labels["eulerian_info"]->setText(eulerian_text);
}

QString MainWindow::updatePruferInfo() {
    auto [pruferCode, weights] = graph->encodePrufer();
    auto decodedEdges = graph->decodePrufer(pruferCode, weights);

    QString pruferInfo;
    pruferInfo += u8"Код Прюфера данного остова: (" + joinQVector(pruferCode) + ")\n";
    pruferInfo += u8"Код Прюфера весов данного остова: (" + joinQVector(weights) + ")\n";
    pruferInfo += u8"Декодирование: " + formatDecodedEdges(decodedEdges);

    return pruferInfo;
}

QString MainWindow::joinQVector(const QVector<int>& vec) const {
    QStringList elements;
    for (int value : vec) {
        elements << QString::number(value);
    }
    return elements.join(", ");
}

QString MainWindow::formatDecodedEdges(const QVector<QPair<QPair<int, int>, int>>& edges) const {
    QStringList edgeStrings;
    for (const auto& [vertices, weight] : edges) {
        edgeStrings << QString("(%1, %2): %3").arg(vertices.first)
            .arg(vertices.second)
            .arg(weight);
    }
    return edgeStrings.join(", ");
}

void MainWindow::createMatricesTab(QWidget* tab) {
    QVBoxLayout* layout = new QVBoxLayout(tab);

    layout->addWidget(tables["adjacency"].label);
    layout->addWidget(tables["adjacency"].widget);

    layout->addWidget(tables["weight"].label);
    layout->addWidget(tables["weight"].widget);
}

void MainWindow::createShimbellTab(QWidget* tab) {
    QVBoxLayout* layout = new QVBoxLayout(tab);

    layout->addWidget(tables["min_path"].label);
    layout->addWidget(tables["min_path"].widget);

    layout->addWidget(tables["max_path"].label);
    layout->addWidget(tables["max_path"].widget);
}

void MainWindow::createWarshallTab(QWidget* tab) {
    QVBoxLayout* layout = new QVBoxLayout(tab);

    layout->addWidget(labels["warshall"]);
    layout->addWidget(labels["dfs"]);
    layout->addWidget(labels["bellman_ford"]);
    layout->addWidget(labels["iterations"]);
}

void MainWindow::createFlowTab(QWidget* tab) {
    QVBoxLayout* layout = new QVBoxLayout(tab);

    layout->addWidget(labels["flow"]);
    layout->addWidget(labels["min_cost_flow"]);

    layout->addWidget(tables["flow"].label);
    layout->addWidget(tables["flow"].widget);

    layout->addWidget(tables["cost"].label);
    layout->addWidget(tables["cost"].widget);
}

void MainWindow::createSpanningTreeTab(QWidget* tab) {
    QVBoxLayout* layout = new QVBoxLayout(tab);

    layout->addWidget(labels["spanning_tree_count"]);
    layout->addWidget(tables["kirchhoff"].label);
    layout->addWidget(tables["kirchhoff"].widget);

    layout->addSpacing(116);
    layout->addWidget(labels["min_spanning_tree"]);
    layout->addWidget(labels["prufer_code"]);
    layout->addSpacing(116);
    layout->addWidget(labels["colors"]);
    layout->addSpacing(116);
}

void MainWindow::createEulerianTab(QWidget* tab) {
    QVBoxLayout* layout = new QVBoxLayout(tab);

    layout->addWidget(labels["eulerian_info"]);

    QPushButton* find_euler_cycle_button = new QPushButton(u8"Найти Эйлеров цикл", this);
    layout->addWidget(find_euler_cycle_button);
    connect(find_euler_cycle_button, &QPushButton::clicked, this, &MainWindow::findEulerianCycle);

    QPushButton* find_hamilton_cycle_button = new QPushButton(u8"Найти Гамильтонов цикл", this);
    layout->addWidget(find_hamilton_cycle_button);
    connect(find_hamilton_cycle_button, &QPushButton::clicked, this, &MainWindow::findHamiltonianCycle);
}

void MainWindow::setupConnections() {
    connect(generate_button, &QPushButton::clicked, this, &MainWindow::generateGraph);
    connect(calculate_button, &QPushButton::clicked, this, &MainWindow::calculateShimbell);
    connect(warshall_calculate_button, &QPushButton::clicked, this, [this]() {
        tab_widget->currentIndex() == 2 ? calculateWarshall() : calculateFlow();
        });
    connect(tab_widget, &QTabWidget::currentChanged, this, &MainWindow::onTabChanged);
}

void MainWindow::setupTabWidget() {
    QWidget* matricesTab = new QWidget(this);
    createMatricesTab(matricesTab);
    tab_widget->addTab(matricesTab, u8"Матрицы графа");

    QWidget* shimbellTab = new QWidget(this);
    createShimbellTab(shimbellTab);
    tab_widget->addTab(shimbellTab, u8"Метод Шимбелла");

    QWidget* warshallTab = new QWidget(this);
    createWarshallTab(warshallTab);
    tab_widget->addTab(warshallTab, u8"Маршрут");

    QWidget* flowTab = new QWidget(this);
    createFlowTab(flowTab);
    tab_widget->addTab(flowTab, u8"Поток");

    QWidget* spanningTreeTab = new QWidget(this);
    createSpanningTreeTab(spanningTreeTab);
    tab_widget->addTab(spanningTreeTab, u8"Остовы");

    QWidget* eulerianTab = new QWidget(this);
    createEulerianTab(eulerianTab);
    tab_widget->addTab(eulerianTab, u8"Циклы");
}

void MainWindow::onTabChanged(int index) {
    graph->clearHighlight();
    graph->clearEulerian();
    graph->clearHamiltonian();

    vertex_input->setVisible(index == 0 || index == 4 || index == 5);
    generate_button->setVisible(index == 0 || index == 4 || index == 5);

    shimbell_edges_input->setVisible(index == 1);
    calculate_button->setVisible(index == 1);

    bool showWarshallInputs = index == 2 || index == 3;
    warshall_first_input->setVisible(showWarshallInputs);
    warshall_second_input->setVisible(showWarshallInputs);
    warshall_calculate_button->setVisible(showWarshallInputs);

    graph->setShowFlowInfo(index == 3);
    graph->setNeorFlag(index == 4 || index == 5);
    graph->setShowMST(index == 4);

    graph->update();
    current_tab = index;
}

void MainWindow::generateGraph() {
    QString vertexInputText = vertex_input->text();
    if (vertexInputText.isEmpty() || vertexInputText.toInt() < 2 || vertexInputText.toInt() > 100) return;

    delete graph;
    graph = new Graph(vertexInputText.toInt());
    this->graph_layout->addWidget(graph);

    QHash<QString, TableData>::iterator it;
    for (it = tables.begin(); it != tables.end(); ++it) {
        if (it.value().getter && it.value().widget) {
            it.value().widget->RegenerateTable((graph->*(it.value().getter))(), this, double_tables_k);
        }
    }

    updateDynamicLabels();
    updateGraphViewSettings();
}

void MainWindow::updateGraphViewSettings() {
    graph->setShowFlowInfo(current_tab == 3);
    graph->setNeorFlag(current_tab == 4 || current_tab == 5);
    graph->setShowMST(current_tab == 4);
    graph->update();
}


void MainWindow::calculateShimbell() {
    QString edgesText = shimbell_edges_input->text();
    if (edgesText.isEmpty() || edgesText.toInt() < 1 || edgesText.toInt() > graph->getVertexesCount()) return;

    int edgesCount = edgesText.toInt();

    QVector<QVector<int>>* minMatrix = graph->getMinWayByShimbell(edgesCount);
    tables["min_path"].widget->RegenerateTable(*minMatrix, this, double_tables_k);
    delete minMatrix;

    QVector<QVector<int>>* maxMatrix = graph->getMaxWayByShimbell(edgesCount);
    tables["max_path"].widget->RegenerateTable(*maxMatrix, this, double_tables_k);
    delete maxMatrix;
}

void MainWindow::calculateWarshall() {
    QString warshall_first_input = this->warshall_first_input->text();
    QString warshall_second_input = this->warshall_second_input->text();

    if (warshall_first_input.isEmpty() || warshall_first_input.toInt() < 0 || warshall_first_input.toInt() > this->graph->getVertexesCount() - 1) return;
    if (warshall_second_input.isEmpty() || warshall_second_input.toInt() < 0 || warshall_second_input.toInt() > this->graph->getVertexesCount() - 1) return;

    int first_vertex_id = this->warshall_first_input->text().toInt();
    int second_vertex_id = this->warshall_second_input->text().toInt();

    if (first_vertex_id == second_vertex_id) {
        labels["warshall"]->setText(
            u8"1) Алгоритм Уоршалла \n Вы уже находитесь в данной вершине!");
    }

    else {
        int pathCount = graph->FindPathCount(first_vertex_id, second_vertex_id);
        labels["warshall"]->setText(
            u8"1) Алгоритм Уоршалла \n Количество маршрутов из вершины " +
            warshall_first_input + u8" в вершину " + warshall_second_input +
            u8" : " + QString::number(pathCount));
    }

    QString way;
    int dfsIterations = graph->DFS(way, first_vertex_id);
    labels["dfs"]->setText(
        u8"2) Обход графа по ребрам в глубину из вершины " + warshall_first_input + " :\n" + way);

    QVector<QString> distances(this->graph->getVertexesCount());
    QVector<int> path;
    int bfIterations = graph->BellmanFord(first_vertex_id, second_vertex_id, distances, path);

    QString bfResult;
    if (bfIterations == -1) {
        bfResult = u8"3) Алгоритм Беллмана-Форда \n Найден отрицательный цикл!";
    }
    else {
        bfResult = u8"3) Алгоритм Беллмана-Форда \n Кратчайшее расстояние из вершины " + warshall_first_input + u8" в вершину " + warshall_second_input + " : " +
            distances[second_vertex_id] + u8"\nКратчайший путь: ";
        for (int i = 0; i < path.size(); ++i) {
            bfResult += QString::number(path[i]);
            if (i < path.size() - 1) bfResult += " -> ";
        }
        bfResult += "\n";

        bfResult += QString::fromUtf8(u8"Вектор расстояний от вершины ") + warshall_first_input + ": (";
        for (int i = 0; i < distances.size() - 1; ++i) {
            bfResult += distances[i] + ", ";
        }
        bfResult += distances[distances.size() - 1] + ")";
    }
    labels["bellman_ford"]->setText(bfResult);
    this->graph->update();

    labels["iterations"]->setText(QString::fromUtf8(u8"Количество итераций алгоритма 2: ") +
        QString::number(dfsIterations) + QString::fromUtf8(u8"\n Количество итераций алгоритма 3: ") + QString::number(bfIterations));
}

void MainWindow::calculateFlow() {
    QString first_input = this->warshall_first_input->text();
    QString second_input = this->warshall_second_input->text();
    if (first_input.isEmpty() || first_input.toInt() < 0 || first_input.toInt() > this->graph->getVertexesCount() - 1) return;
    if (second_input.isEmpty() || second_input.toInt() < 0 || second_input.toInt() > this->graph->getVertexesCount() - 1) return;

    int source = first_input.toInt();
    int sink = second_input.toInt();

    QString base_flow_text = u8"Максимальный поток по алгоритму Форда-Фалкерсона : ";

    if (source == sink) {
        labels["flow"]->setText(base_flow_text + u8"Источник и сток совпадают!");
        return;
    }

    bool is_wrong_sink_or_source = false;
    for (int i = 0; i < this->graph->getVertexesCount(); i++) {
        if (this->graph->getAdjacencyMatrix()[i][source] || this->graph->getAdjacencyMatrix()[sink][i]) {
            is_wrong_sink_or_source = true;
            break;
        }
    }

    if (is_wrong_sink_or_source) {
        labels["flow"]->setText(base_flow_text + u8"Вершина не является истоком или стоком!");
        return;
    }
    
    int maxFlow = graph->FordFulkerson(source, sink);
    labels["flow"]->setText(base_flow_text + QString::number(maxFlow));

    QVector<QVector<int>> flowMatrix(graph->getVertexesCount(), QVector<int>(graph->getVertexesCount(), 0));
    auto [flow, cost] = graph->FindMinCostFlow(source, sink, flowMatrix, maxFlow);

    tables["flow"].widget->RegenerateTable(flowMatrix, this, double_tables_k);
    labels["min_cost_flow"]->setText(u8"Минимальная стоимость величины потока " + QString::number(flow) + " : " + QString::number(cost));
}

void MainWindow::findEulerianCycle() {
    if (this->graph->getVertexesCount() == 2) return;

    graph->clearHamiltonian();
    graph->makeEulerian();
    auto cycle = graph->findEulerianCycle();

    QString cycleText;
    for (const auto& edge : cycle) {
        cycleText += QString("(%1, %2), ").arg(edge.first).arg(edge.second);
    }
    labels["eulerian_info"]->setText(u8"Эйлеров цикл: " + cycleText);

    graph->update();
}

void MainWindow::findHamiltonianCycle() {
    if (this->graph->getVertexesCount() == 2) return;

    graph->clearEulerian();
    
    if (!graph->isHamiltonian()) graph->makeHamiltonian();
    auto result = graph->solveTSP("123.txt");

    QString cycleText;
    int tmp_vertex = -1;
    for (int i = 0; i < result.second.size(); ++i) {
        if (tmp_vertex != -1)
            cycleText += QString("(%1, %2), ").arg(tmp_vertex).arg(result.second[i]);
        tmp_vertex = result.second[i];
    }

    labels["eulerian_info"]->setText(
        u8"Минимальное пройденное расстояние для задачи коммивояжёра: " +
        QString::number(result.first) +
        u8"\nМинимальный цикл: " + cycleText);
    graph->update();
}

TableWidget::TableWidget(const QVector<QVector<int>>& matrix, QWidget* parent, double k, const QColor& highlightColor) 
    : QTableWidget(parent) {

    this->highlightColor = highlightColor;

    RegenerateTable(matrix, parent, k);
    setAutoScroll(false);
    horizontalHeader()->setStretchLastSection(false);
    verticalHeader()->setStretchLastSection(false);
    setEditTriggers(QAbstractItemView::NoEditTriggers);
}

void TableWidget::RegenerateTable(const QVector<QVector<int>>& matrix, QWidget* parent, double k) {
    int n = matrix.size();
    setRowCount(n);
    setColumnCount(n);

    int columnWidth = (parent->width() / 2) / n - (parent->width() * k / n);

    for (int j = 0; j < n; ++j) {
        setColumnWidth(j, columnWidth);
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            QTableWidgetItem* item = new QTableWidgetItem(QString::number(matrix[i][j]));

            if (matrix[i][j] != 0) {
                item->setBackground(QBrush(this->highlightColor));
            }
            setItem(i, j, item);
        }
    }

    QStringList headers = GenerateHeaderLabels(n);
    setHorizontalHeaderLabels(headers);
    setVerticalHeaderLabels(headers);
}

QStringList TableWidget::GenerateHeaderLabels(int size) {
    QStringList headers;
    for (int i = 0; i < size; ++i) {
        headers << QString::number(i);
    }
    return headers;
}
