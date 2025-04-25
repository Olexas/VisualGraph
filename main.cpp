#include <QApplication>
#include <iostream>
#include <vector>

#include "MainWindow.h"
#include "Graph.h"  

int main(int argc, char* argv[]) {
    QApplication a(argc, argv);
    
    MainWindow mainWindow;
    mainWindow.show();

    return a.exec();
}