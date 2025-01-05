#include <iostream>
#include <vector>
#include <stdio.h>
#include <fstream>
#include <algorithm>
#include <iomanip>

using namespace std;

// Координаты узлов (r,z)
vector<pair<double, double>> Cn;

// Материал
vector <pair<double, double>> mat;

// Элементы разбиения (нижняя, верхняя грани, номер материала)
vector<vector<int>> elems;

// Первое краевое условие (глобальный номер узла, значение функции)
vector <pair<int, double>> b1;

// Второе краевое условие (номер элемента, его локальная грань, значение)
vector<vector<double>> b2;

int Nn; // Число узлов
int Nel; // Число элементов
int Nt;//Число временных слоев

vector<int> ig;
vector<int> jg;

vector<double> ggl; //нижний треугольник матрицы слау

vector<double> ggu; //верхний треугольник матрицы слау

vector<double> di; //главная диагональ матрицы слау
vector<double> gggl; //нижний треугольник матрицы жесткости
vector<double> gggu; //верхний треугольник матрицы жесткости
vector<double> gdi; //главная диагональ матрицы жесткости
vector<double> sggl; //нижний треугольник матрицы масс первой производной
vector<double> sggu; //верхний треугольник матрицы масс первой производной
vector<double> sdi; //главная диагональ матрицы масс первой производной
vector<double> vec; //вектор правой части
vector<double> P; // Результат
vector<double> t; //вектор временных слоев
vector<double> P0; //вектор весов функции на первом слое
vector<double> P1; //вектор весов функции на втором слое
vector<double> P2; //вектор весов функции на третьем слое
vector<double> nu; //вектор коэффициентов для первой производной

double Function(double r, double z, double t) //функция правой части
{
    return 1;
}
double Tetta(int n, double r, double z)
{
    switch (n)
    {
        case 1:
            return 1;
        case 0:
            return 0;
        case -1:
            return 1;
        case -2:
            return -2;
    }
}

double Sigma(int n, double r, double z)//Функция сигмы
{
    switch (n)
    {
        case 0:
            return 0;
        case 1:
            return 1;
        case 2:
            return r;
        case 3:
            return 3;
        case 4:
            return 4;
    }
}
double U(int n, double r, double z, double t)//функция третьего краевого условия
{
    switch (n)
    {
    case 0:
        return 0;
    case 1:
        return t;
    case 2:
        return pow(t, 4);
    case 3:
        return 3;
    case 4:
        return 4;
    }
}