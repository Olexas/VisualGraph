#include "MainWindow.h"

MainWindow::MainWindow() {
    resize(1920, 1080);
    //QWidget::showFullScreen();

    QWidget* centralWidget = new QWidget(this);
    QHBoxLayout* main_layout = new QHBoxLayout(centralWidget);

    this->inputWidget = new QWidget(this);
    QHBoxLayout* inputLayout = new QHBoxLayout(inputWidget);

    this->vertex_input = new QLineEdit(inputWidget);
    vertex_input->setPlaceholderText(QString::fromUtf8(u8"Введите количество вершин (от 2 до 100)"));
    vertex_input->setStyleSheet("color: grey;");
    inputLayout->addWidget(vertex_input);

    QRegularExpression vertexRegex(QStringLiteral("^[0-9]*$"));
    QRegularExpressionValidator* vertexValidator = new QRegularExpressionValidator(vertexRegex, this);
    vertex_input->setValidator(vertexValidator);

    this->shimbell_edges_input = new QLineEdit(inputWidget);
    shimbell_edges_input->setStyleSheet("color: grey;");
    shimbell_edges_input->setVisible(false);
    inputLayout->addWidget(shimbell_edges_input);

    shimbell_edges_input->setValidator(vertexValidator);

    this->warshall_first_input = new QLineEdit(inputWidget);
    this->warshall_second_input = new QLineEdit(inputWidget);
    warshall_first_input->setPlaceholderText(QString::fromUtf8(u8"Введите номер первой вершины"));
    warshall_second_input->setPlaceholderText(QString::fromUtf8(u8"Введите номер второй вершины"));
    warshall_first_input->setStyleSheet("color: grey;");
    warshall_second_input->setStyleSheet("color: grey;");
    warshall_first_input->setVisible(false);
    warshall_second_input->setVisible(false);
    inputLayout->addWidget(warshall_first_input);
    inputLayout->addWidget(warshall_second_input);

    warshall_first_input->setValidator(vertexValidator);
    warshall_second_input->setValidator(vertexValidator);

    this->generate_button = new QPushButton(QString::fromUtf8(u8"Сгенерировать"), inputWidget);
    inputLayout->addWidget(generate_button);
    connect(generate_button, &QPushButton::clicked, this, &MainWindow::generateGraph);

    this->calculate_button = new QPushButton(QString::fromUtf8(u8"Рассчитать"), inputWidget);
    calculate_button->setVisible(false);
    inputLayout->addWidget(calculate_button);
    connect(calculate_button, &QPushButton::clicked, this, &MainWindow::calculateShimbell);

    this->warshall_calculate_button = new QPushButton(QString::fromUtf8(u8"Рассчитать"), inputWidget);
    warshall_calculate_button->setVisible(false);
    inputLayout->addWidget(warshall_calculate_button);

    this->graph = new Graph(10);

    shimbell_edges_input->setPlaceholderText(QString::fromUtf8(u8"Введите количество рёбер (от 1 до ") + QString::number(this->graph->getVertexesCount()) + ")");

    this->table_widget = new QTabWidget(this);
    connect(table_widget, &QTabWidget::currentChanged, this, &MainWindow::onTabChanged);

    this->adjacency_table = new TableWidget(this->graph->getAdjacencyMatrix(), this);
    this->weight_table = new TableWidget(this->graph->getWeightMatrix(), this);

    QWidget* shimbellTab = new QWidget(this);
    QVBoxLayout* shimbellLayout = new QVBoxLayout(shimbellTab);

    this->min_path_table = new TableWidget(QVector<QVector<int>>(this->graph->getVertexesCount(), QVector<int>(this->graph->getVertexesCount())), this, double_tables_k);
    this->max_path_table = new TableWidget(QVector<QVector<int>>(this->graph->getVertexesCount(), QVector<int>(this->graph->getVertexesCount())), this, double_tables_k);

    QLabel* minPathLabel = new QLabel(QString::fromUtf8(u8"Минмальный путь"), this);
    QLabel* maxPathLabel = new QLabel(QString::fromUtf8(u8"Максимальный путь"), this);

    shimbellLayout->addWidget(minPathLabel);
    shimbellLayout->addWidget(min_path_table);
    shimbellLayout->addWidget(maxPathLabel);
    shimbellLayout->addWidget(max_path_table);

    QWidget* warshallTab = new QWidget(this);
    QVBoxLayout* warshallLayout = new QVBoxLayout(warshallTab);

    this->warshall_label = new QLabel(QString::fromUtf8(u8"1) Алгоритм Уоршалла \n Количество маршрутов из вершины ... в вершину ... : "));
    warshall_label->setStyleSheet("font-size: 12pt;");
    warshall_label->setAlignment(Qt::AlignHCenter);
    warshallLayout->addWidget(warshall_label);
        
    this->DFS_label = new QLabel(QString::fromUtf8(u8"2) Поиск в глубину \n Количество маршрутов из вершины ... в вершину ... : "));
    DFS_label->setStyleSheet("font-size: 12pt;");
    DFS_label->setAlignment(Qt::AlignHCenter);
    warshallLayout->addWidget(DFS_label);
    DFS_label->setWordWrap(true);

    this->bellmanford_label = new QLabel(QString::fromUtf8(u8"3) Алгоритм Беллмана-Форда \n Количество маршрутов из вершины ... в вершину ... : "));
    bellmanford_label->setStyleSheet("font-size: 12pt;");
    bellmanford_label->setAlignment(Qt::AlignHCenter);
    warshallLayout->addWidget(bellmanford_label);
    bellmanford_label->setWordWrap(true);

    this->iterations_label = new QLabel(QString::fromUtf8(u8"Количество итераций алгоритма 2: ... \n Количество итераций алгоритма 3: ..."));
    iterations_label->setStyleSheet("font-size: 12pt;");
    iterations_label->setAlignment(Qt::AlignHCenter);
    warshallLayout->addWidget(iterations_label);

    QWidget* flowTab = new QWidget(this);
    QVBoxLayout* flowLayout = new QVBoxLayout(flowTab);

    this->flow_label = new QLabel(QString::fromUtf8(u8"Максимальный поток по алгоритму Форда-Фалкерсона : ..."), this);
    flow_label->setStyleSheet("font-size: 11pt;");
    flow_label->setAlignment(Qt::AlignHCenter);
    flowLayout->addWidget(flow_label);

    this->min_cost_flow_label = new QLabel(QString::fromUtf8(u8"Минимальная стоимость потока ... : ..."), this);
    min_cost_flow_label->setStyleSheet("font-size: 11pt;");
    min_cost_flow_label->setAlignment(Qt::AlignHCenter);
    flowLayout->addWidget(min_cost_flow_label);

    QLabel* capacityLabel = new QLabel(QString::fromUtf8(u8"Матрица пропускных способностей"), this);
    QLabel* costLabel = new QLabel(QString::fromUtf8(u8"Матрица стоимостей"), this);

    this->capacity_table = new TableWidget(this->graph->getCapacityMatrix(), this, double_tables_k, QColor(255, 255, 200));
    this->cost_table = new TableWidget(this->graph->getCostMatrix(), this, double_tables_k, QColor(144, 238, 144));

    flowLayout->addWidget(capacityLabel);
    flowLayout->addWidget(capacity_table);
    flowLayout->addWidget(costLabel);
    flowLayout->addWidget(cost_table);

    connect(warshall_calculate_button, &QPushButton::clicked, this, [this]() {
        if (table_widget->currentIndex() == 3) {
            calculateWarshall();
        }
        else if (table_widget->currentIndex() == 4) {
            calculateFlow();
        }
        });
  
    table_widget->addTab(adjacency_table, QString::fromUtf8(u8"Матрица смежности"));
    table_widget->addTab(weight_table, QString::fromUtf8(u8"Весовая матрица"));
    table_widget->addTab(shimbellTab, QString::fromUtf8(u8"Метод Шимбелла"));
    table_widget->addTab(warshallTab, QString::fromUtf8(u8"Маршрут"));
    table_widget->addTab(flowTab, QString::fromUtf8(u8"Поток"));

    adjacency_table->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

    this->graph_layout = new QVBoxLayout();
    graph_layout->addWidget(inputWidget);
    graph_layout->addWidget(graph);

    main_layout->addLayout(graph_layout);
    main_layout->addWidget(table_widget);

    setCentralWidget(centralWidget);
}

void MainWindow::generateGraph() {
    QString vertex_input_ = this->vertex_input->text();
    if (vertex_input_.isEmpty() || vertex_input_.toInt() < 2 || vertex_input_.toInt() > 100) return;

    delete this->graph;
    this->graph = new Graph(vertex_input_.toInt());
    this->graph_layout->addWidget(graph);

    this->adjacency_table->RegenarateTable(this->graph->getAdjacencyMatrix(), this);

    this->weight_table->RegenarateTable(this->graph->getWeightMatrix(), this);

    this->min_path_table->RegenarateTable(QVector<QVector<int>>(this->graph->getVertexesCount(), QVector<int>(this->graph->getVertexesCount())), this, double_tables_k);
    this->max_path_table->RegenarateTable(QVector<QVector<int>>(this->graph->getVertexesCount(), QVector<int>(this->graph->getVertexesCount())), this, double_tables_k);

    shimbell_edges_input->setPlaceholderText(QString::fromUtf8(u8"Введите количество ребер (от 1 до ") + QString::number(this->graph->getVertexesCount()) + ")");

    this->warshall_label->setText(QString::fromUtf8(u8"1) Алгоритм Уоршалла \n Количество маршрутов из вершины ... в вершину ... : "));
    this->DFS_label->setText(QString::fromUtf8(u8"2) Поиск в глубину \n Количество маршрутов из вершины ... в вершину ... : "));
    this->bellmanford_label->setText(QString::fromUtf8(u8"3) Алгоритм Беллмана-Форда \n Количество маршрутов из вершины ... в вершину ... : "));
    this->iterations_label->setText(QString::fromUtf8(u8"Количество итераций алгоритма 2: ... \n Количество итераций алгоритма 3: ..."));

    capacity_table->RegenarateTable(this->graph->getCapacityMatrix(), this, double_tables_k);
    cost_table->RegenarateTable(this->graph->getCostMatrix(), this, double_tables_k);

    this->flow_label->setText(QString::fromUtf8(u8"Максимальный поток по алгоритму Форда-Фалкерсона : ..."));
    this->min_cost_flow_label->setText(QString::fromUtf8(u8"Минимальная стоимость потока ... : ..."));
}   

void MainWindow::calculateShimbell() {
    QString shimbell_edges_input_ = this->shimbell_edges_input->text();
    if (shimbell_edges_input_.isEmpty() || shimbell_edges_input_.toInt() < 1 || shimbell_edges_input_.toInt() > this->graph->getVertexesCount()) return;

    int edges_count = shimbell_edges_input_.toInt();

    QVector<QVector<int>>* min_matrix = this->graph->getMinWayByShimbell(edges_count);
    this->min_path_table->RegenarateTable(*min_matrix, this, double_tables_k);
    delete min_matrix;

    QVector<QVector<int>>* max_matrix = this->graph->getMaxWayByShimbell(edges_count);
    this->max_path_table->RegenarateTable(*max_matrix, this, double_tables_k);
    delete max_matrix;
}

void MainWindow::onTabChanged(int index) {
    if (index == 4) {
        this->graph->setShowFlowInfo(true);
        this->graph->update();
    }

    else {
        this->graph->setShowFlowInfo(false);
        this->graph->update();
    }

    if (index == 2) {
        vertex_input->setVisible(false);
        shimbell_edges_input->setVisible(true);
        generate_button->setVisible(false);
        calculate_button->setVisible(true);
        warshall_first_input->setVisible(false);
        warshall_second_input->setVisible(false);
        warshall_calculate_button->setVisible(false);
    }

    else if (index == 3 || index == 4) {
        vertex_input->setVisible(false);
        shimbell_edges_input->setVisible(false);
        generate_button->setVisible(false);
        calculate_button->setVisible(false);
        warshall_first_input->setVisible(true);
        warshall_second_input->setVisible(true);
        warshall_calculate_button->setVisible(true);
    }

    else {
        vertex_input->setVisible(true);
        shimbell_edges_input->setVisible(false);
        generate_button->setVisible(true);
        calculate_button->setVisible(false);
        warshall_first_input->setVisible(false);
        warshall_second_input->setVisible(false);
        warshall_calculate_button->setVisible(false);
    }
}

void MainWindow::calculateWarshall() {
    QString warshall_first_input = this->warshall_first_input->text();
    QString warshall_second_input = this->warshall_second_input->text();
    if (warshall_first_input.isEmpty() || warshall_first_input.toInt() < 0 || warshall_first_input.toInt() > this->graph->getVertexesCount() - 1) return;
    if (warshall_second_input.isEmpty() || warshall_second_input.toInt() < 0 || warshall_second_input.toInt() > this->graph->getVertexesCount() - 1) return;

    int first_vertex_id = this->warshall_first_input->text().toInt();
    int second_vertex_id = this->warshall_second_input->text().toInt();

    int path_count;
    bool firs_second_flag = false;
    if (first_vertex_id == second_vertex_id) {
        firs_second_flag = true;
        this->warshall_label->setText(QString::fromUtf8(u8"1) Алгоритм Уоршалла \n Вы уже находитесь в данной вершине!"));
    }
    else {
        path_count = this->graph->FindPathCount(first_vertex_id, second_vertex_id);

        this->warshall_label->setText(QString::fromUtf8(u8"1) Алгоритм Уоршалла \n Количество маршрутов из вершины ") +
            warshall_first_input + QString::fromUtf8(u8" в вершину ") + warshall_second_input + QString::fromUtf8(u8" : ") + QString::number(path_count));
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    QString way;
    int iterations2 = this->graph->DFS(way, first_vertex_id);
    QString final = QString::fromUtf8(u8"2) Обход графа по ребрам в глубину из вершины ") + warshall_first_input + " : " + "\n" + way;
    this->DFS_label->setText(final);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    QString final2;  int iterations3 = 0;

    if (firs_second_flag) final2 = final2 = QString::fromUtf8(u8"3) Алгоритм Беллмана-Форда \n Вы уже находитесь в данной вершине!");

    else if (path_count == 0) {
        final2 = QString::fromUtf8(u8"3) Алгоритм Беллмана-Форда \n Маршрут из вершины ") +
            warshall_first_input + QString::fromUtf8(u8" в вершину ") + warshall_second_input + QString::fromUtf8(u8" не существует!");
    }

    else {
        QVector<int> distances;
        QVector<int> path;

        iterations3 = this->graph->BellmanFord(first_vertex_id, second_vertex_id, distances, path);
        graph->update();

        if (iterations3 == -1) final2 = QString::fromUtf8(u8"3) Алгоритм Беллмана-Форда \n Найден отрицательный цикл! ");

        else {
            final2 = QString::fromUtf8(u8"3) Алгоритм Беллмана-Форда \n Кратчайшее расстояние из вершины ") +
                warshall_first_input + QString::fromUtf8(u8" в вершину ") + warshall_second_input + QString::fromUtf8(u8" : ") + QString::number(distances[second_vertex_id]) + "\n";

            final2 += QString::fromUtf8(u8"Кратчайший путь: ");
            for (int i = 0; i < path.size(); ++i) {
                final2 += QString::number(path[i]);
                if (i != path.size() - 1) {
                    final2 += " -> ";
                }
            }
            final2 += "\n";

            final2 += QString::fromUtf8(u8"Вектор расстояний от вершины ") + warshall_first_input + ": (";
            for (int i = 0; i < distances.size() - 1; ++i) {
                final2 += QString::number(distances[i]) + ", ";
            }
            final2 += QString::number(distances[distances.size() - 1]) + ")";
        }
    }

    this->bellmanford_label->setText(final2);

    this->iterations_label->setText(QString::fromUtf8(u8"Количество итераций алгоритма 2: ") +
        QString::number(iterations2) + QString::fromUtf8(u8"\n Количество итераций алгоритма 3: ") + QString::number(iterations3));
}

void MainWindow::calculateFlow() {
    QString first_input = this->warshall_first_input->text();
    QString second_input = this->warshall_second_input->text();
    if (first_input.isEmpty() || first_input.toInt() < 0 || first_input.toInt() > this->graph->getVertexesCount() - 1) return;
    if (second_input.isEmpty() || second_input.toInt() < 0 || second_input.toInt() > this->graph->getVertexesCount() - 1) return;

    int source = first_input.toInt();
    int sink = second_input.toInt();

    bool is_wrong_sink_or_source = false;
    for (int i = 0; i < this->graph->getVertexesCount(); i++) {
        if (this->graph->getAdjacencyMatrix()[i][source] || this->graph->getAdjacencyMatrix()[sink][i]) {
            is_wrong_sink_or_source = true;
            break;
        }
    }

    if (is_wrong_sink_or_source) {
        flow_label->setText(QString::fromUtf8(u8"Максимальный поток по алгоритму Форда-Фалкерсона : Вершина не является истоком или стоком!"));
        return;
    }

    if (source == sink) {
        flow_label->setText(QString::fromUtf8(u8"Максимальный поток по алгоритму Форда-Фалкерсона : Источник и сток совпадают!"));
        return;
    }

    int max_flow = this->graph->FordFulkerson(source, sink);
    flow_label->setText(QString::fromUtf8(u8"Максимальный поток по алгоритму Форда-Фалкерсона : ") + QString::number(max_flow));

    auto pair = this->graph->FindMinCostFlow(source, sink, max_flow);
    int flow = pair.first; int cost = pair.second;

    this->min_cost_flow_label->setText(QString::fromUtf8(u8"Минимальная стоимость потока ") + QString::number(flow) + " : " + QString::number(cost));

    this->graph->update();
}


TableWidget::TableWidget(const QVector<QVector<int>>& matrix, QWidget* parent, double k, const QColor& highlightColor) : QTableWidget(parent) {
    this->highlightColor = highlightColor;
    this->RegenarateTable(matrix, parent, k);

    setAutoScroll(false);
    horizontalHeader()->setStretchLastSection(false);
    verticalHeader()->setStretchLastSection(false);
    setEditTriggers(QAbstractItemView::NoEditTriggers);
}

void TableWidget::RegenarateTable(const QVector<QVector<int>>& matrix, QWidget* parent, double k) {
    int n = matrix.size();
    setRowCount(n);
    setColumnCount(n);

    int totalWidth = (*parent).width() / 2;
    int columnWidth = totalWidth / n - ((*parent).width() * k / n);

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

    setHorizontalHeaderLabels(getHeaderLabels(n));
    setVerticalHeaderLabels(getHeaderLabels(n));
}

QStringList TableWidget::getHeaderLabels(int size) {
    QStringList headers;
    for (int i = 0; i < size; ++i) {
        headers << QString::number(i);
    }
    return headers;
}
