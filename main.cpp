#include <iostream>
#include <math.h>
#include <memory.h>
#include <vector>
#include<chrono>
//#include<omp.h>

#include<SFML/Graphics.hpp>
#include "Menu.hpp"
#include "Theme.hpp"
#include "Gui.hpp"

using namespace std;
using namespace std::chrono;

class CrsMatrix{
public:
    int N;                  // размер матрицы NxN
    int NZ;                 // Количество не нулевых элементов
    vector<double> Value;   // Массив значений (len = NZ)
    vector<int> Col;        // Массив номеров столбцов (len = NZ)
    vector<int> RowIndex;   // массив индексов строк (len = N + 1)
};

// скалярное произведение вектора А на вектор В
double MuliplicateVV(vector<double> A,      //вектор А
                    vector<double> B,       //вектор В
                    int size){              //размер векторов
    double result = 0;
//    #pragma omp parallel for reduction (+:result)
    for(int i = 0; i < size; i++){
        result += A[i] * B[i];
    }
    return result;
}

// произведение матрицы M на вектор V
vector<double> MuliplicateMV(CrsMatrix M,           //матрица M
                   vector<double> V){
    int i = 0;
    int j = 0;
    vector<double> result(V.size(), 0);
//    #pragma omp parallel for private(i, j)
    for (i = 0; i < M.N; i++){
        for(j = M.RowIndex[i]+1; j<M.RowIndex[i+1]; j++){
//            #pragma omp atomic
            result[i] += M.Value[j] * V[M.Col[j]];
        }
//        #pragma omp atomic
        result[i] += M.Value[M.RowIndex[i]] * V[M.Col[M.RowIndex[i]]];
    }
    return result;
}

vector<double>  CGMethod (CrsMatrix A,   //матрица
               vector<double> b,         //вектор правой части
               vector<double> x0,        //начальное приблежение
               float eps,                //требуемая точность
               int maxIter)              //максимальное число итераций

{
    vector<double> Xk(x0.size(), 0);     //вектор ответов
    int count = 0;                       //число выполненных итераций

    int size = A.N;
    vector<double> Rk_1(x0.size(), 0);
    vector<double> Rk(x0.size(), 0);
    vector<double> Zk_1(x0.size(), 0);
    vector<double> Ax0(x0.size(), 0);
    vector<double> AZk_1(x0.size(), 0);
    float beta,alpha,check,norm;
    norm = sqrt(MuliplicateVV(b,b,size));       // Вычисляем сумму квадратов элементов вектора F

    Ax0 = MuliplicateMV(A, x0);
    for (int j = 0; j < size; j++){
        Rk_1[j] = b[j] - Ax0[j];                //Задаем начальное значение r0
    }

    for (int j = 0; j < size; j++){
        Zk_1[j] = Rk_1[j];                      //Задаем начальное значение z0.
        Xk[j] = x0[j];                          // Задаем начальное приближение корней
    }

    do                                         //итерации поиска решения
    {
        (count)++;
        AZk_1 = MuliplicateMV(A, Zk_1);
        alpha = MuliplicateVV(Rk_1, Rk_1,size) /  MuliplicateVV(Zk_1,AZk_1,size);

        for(int j = 0; j < size; j++){
            Xk[j] += alpha * Zk_1[j];
            Rk[j] = Rk_1[j] - alpha * AZk_1[j];
        }

        beta = MuliplicateVV(Rk, Rk, size)/MuliplicateVV(Rk_1, Rk_1, size);
        check = sqrt(MuliplicateVV(Rk, Rk, size));
        for(int j = 0; j < size; j++){
            Zk_1[j] = beta * Zk_1[j] + Rk[j];
        }
        for (int j = 0; j < size; j++){
            Rk_1[j] = Rk[j];
        }
    }
    while ((check/norm > eps) && (count <= maxIter));  //прерываем решение если сошлось или достигли лимита
    return Xk;
}

CrsMatrix getStiffnessMatrix(double c, int BindingCount, double NumLoads, bool log){
    std::vector<double> Val;        //вектор ненулевых значений
    std::vector<int> Col;           //индексы столбцов
    std::vector<int> Row;           //количество ненулевых элементов
    Row.push_back(0);               //нуль - заглушка.
    int counter = 0;
    for(int i = 0; i < NumLoads + 1; i++){                          //строки
        int j;
        int max_col;
        if(log){
            j = 0;
            max_col = NumLoads + 1;
        } else {
            j = i - BindingCount;
            max_col = i + BindingCount + 1;
        }
        //for(int j = 0; j < NumLoads + 1; j++){                        //столбцы
        //for(int j = i - BindingCount; j <= i + BindingCount; j++){    //столбцы
        for(; j < max_col; j++){
            if(j >= 0){
                if(i == 0 && j == 0){
                    Val.push_back(c * 1000000000000);       //заделка
                    Col.push_back(j);
                    counter++;
                    if(log){cout<< "  " << "∞";}
                    continue;
                }
                if(i == NumLoads & j == NumLoads){
                    Val.push_back(c * BindingCount);
                    Col.push_back(j);
                    counter++;
                    if(log){cout  << "  " << c * BindingCount;}
                    continue;
                }
                if(i == j)
                {
                    Col.push_back(j);
                    counter++;
                    int count = BindingCount;
                    if(i < BindingCount || i > NumLoads - BindingCount)
                    {
                        int fix;
                        if(i < BindingCount) fix = BindingCount - i;
                        if(i > NumLoads - BindingCount) fix = BindingCount - (NumLoads - i);
                        count = BindingCount - fix;
                    }
                    Val.push_back(c*BindingCount + c*count);
                    if (log) cout<< "  " << c*BindingCount + c*count;
                    continue;
                }

                if(abs(i - j) < BindingCount+1){
                    Val.push_back(-1*c);
                    Col.push_back(j);
                    counter++;
                    if(log) cout << " " << (-1 * c);
                    continue;
                }
                if(log) cout << "  " << 0;
            }
        }
        Row.push_back(counter);
        if(log) cout << endl;
    }
    CrsMatrix matrixA;
    matrixA.N = NumLoads+1;                 // размер матрицы NxN
    matrixA.NZ = Val.size();                // Количество не нулевых элементов
    matrixA.Col = Col;
    matrixA.RowIndex = Row;
    matrixA.Value = Val;
    return matrixA;
}

struct Theme
{
    sf::Color backgroundColor;
    std::string texturePath;
};

sf::Color hex2color(const std::string& hexcolor)
{
    sf::Color color = sf::Color::Black;
    if (hexcolor.size() == 7) // #ffffff
    {
        color.r = strtoul(hexcolor.substr(1, 2).c_str(), NULL, 16);
        color.g = strtoul(hexcolor.substr(3, 2).c_str(), NULL, 16);
        color.b = strtoul(hexcolor.substr(5, 2).c_str(), NULL, 16);
    }
    else if (hexcolor.size() == 4) // #fff
    {
        color.r = strtoul(hexcolor.substr(1, 1).c_str(), NULL, 16) * 17;
        color.g = strtoul(hexcolor.substr(2, 1).c_str(), NULL, 16) * 17;
        color.b = strtoul(hexcolor.substr(3, 1).c_str(), NULL, 16) * 17;
    }
    return color;
}

void inputScreen(double& c, int& bindingCount, double& numLoads, double& endForce)
{
    Theme defaultTheme = {
        hex2color("#dddbde"),
        "demo/texture-default.png"
    };

    sf::RenderWindow app(sf::VideoMode(300, 320), "Cellular automaton", sf::Style::Close);
    app.setPosition(sf::Vector2i(0, 0));
    gui::Menu menu(app);
    menu.setPosition(10, 10);

    gui::Theme::loadFont("demo/tahoma.ttf");
    gui::Theme::loadTexture(defaultTheme.texturePath);
    gui::Theme::textSize = 11;

    gui::Theme::PADDING = 2.f;
    gui::Theme::windowBgColor = defaultTheme.backgroundColor;

    gui::HBoxLayout* hbox = menu.addHBoxLayout();
    gui::FormLayout* form = hbox->addFormLayout();

    gui::TextBox* stiffTextbox = new gui::TextBox();
    stiffTextbox->setText("3");
    form->addRow("Spring Stiffness: ", stiffTextbox);

    gui::TextBox* bindCntTextbox = new gui::TextBox();
    bindCntTextbox->setText("1");
    form->addRow("Number of Bounds: ", bindCntTextbox);

    gui::TextBox* numLoadsTextbox = new gui::TextBox();
    numLoadsTextbox->setText("10");
    form->addRow("Number of Loads: ", numLoadsTextbox);

    gui::TextBox* endForceTextbox = new gui::TextBox();
    endForceTextbox->setText("100");
    form->addRow("End Force Value: ", endForceTextbox);

    bool open_field = false;

    menu.addButton("Run Calculation!", [&]() {
        open_field = true;
        app.close();
    });

    menu.addButton("Quit", [&]() {
        app.close();
    });

    while (app.isOpen())
    {
        // Process events
        sf::Event event;
        while (app.pollEvent(event))
        {
            // Send events to menu
            menu.onEvent(event);
            if (event.type == sf::Event::Closed)
                app.close();
        }

        // Clear screen
        app.clear(gui::Theme::windowBgColor);
        app.draw(menu);
        // Update the window
        app.display();
        }

    if(open_field)
    {
        c = std::stod((string)stiffTextbox->getText());
        bindingCount = std::stoi((string)bindCntTextbox->getText());
        numLoads = std::stod((string)numLoadsTextbox->getText());
        endForce = std::stod((string)endForceTextbox->getText());
    }
}

int main()
{
    double c;
    int BindingCount;
    double NumLoads;
    double force_at_the_end;

    inputScreen(c, BindingCount, NumLoads, force_at_the_end);

    bool log = true;
    if(NumLoads > 20){
        log = false;
    }

    if(log) cout << "Матрица жесткости: "<< endl;
    CrsMatrix StiffnessMatrix = getStiffnessMatrix(c, BindingCount, NumLoads, log);
    cout << endl;
    std::vector<double> force;
    for(int i = 0; i <= NumLoads; i++){
        if (i != NumLoads)
            force.push_back(0);
        else
            force.push_back(force_at_the_end);
    }
    if (log)
    {
        cout << "Вектор сил: "<< endl;
        for(int i = 0; i < force.size(); i++){
            cout << force[i] << " ";
        }
        cout << endl<< endl;
    }


    std::vector<double> First;    //начальное приблежение
    for(int i = 0; i <= NumLoads; i++){
        First.push_back(0);
    }

    auto start = high_resolution_clock::now();          // DEBUG
    auto result = CGMethod (StiffnessMatrix,            //матрица
              force,                                    //вектор правой части
              First,                                    //начальное приблежение
              0.000001,                                   //требуемая точность
              100);                                      //максимальное число итераций

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    cout << "Время исполнения функции: " << duration.count() <<" сек" <<endl;

    if (log)
    {
        cout << "Вектор перемещений: "<< endl;
        for(int i = 0; i < result.size(); i++){
            cout << result[i] << " ";
        }
        cout << endl<< endl;

        cout << "Контроль: "<< endl;
        auto control = MuliplicateMV(StiffnessMatrix, result);
        for(int i = 0; i < control.size(); i++)
        {
            if(control[i] > 1e-5)
            {
                cout << control[i] << " ";
            } else
            {
                cout << 0 << " ";
            }
        }
        cout << endl<< endl;
    }

    if(!log)
    {
        return 0;
    }
    float t = 0;
    sf::RenderWindow window(sf::VideoMode(1500, 500), "OUR loads!");
    while (window.isOpen())
    {
        t += 0.1;
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }
        window.clear();

        for(int i = 0; i < result.size(); i += 1){
            sf::CircleShape circle;
            circle.setRadius(10);
            circle.setOutlineColor(sf::Color::Red);
            circle.setOutlineThickness(5);
            circle.setPosition(100 + result[i]*0.5* ((sin(t * M_PI/180) + 0.5) / 2) + i*30, 250);
            window.draw(circle);
        }

        window.display();
    }

    return 0;
}
