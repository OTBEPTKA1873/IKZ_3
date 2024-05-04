#include <iostream>
#include <cmath>
using namespace std;

// Отрез кодов друг от друга
void separator()
{
    cout << "\n==================================================================================\n\n";
}
// Функция, которую надо минимизировать
double J(double* X)
{
    return pow(X[0], 2) + pow(X[1], 2);
}
double grad_J(int num_grad, double* X)
{
    switch (num_grad)
    {
    case 0:
        return 2 * X[0];
    case 1:
        return 2 * X[1];
    }
}
// Облегчение написания функции фи
void phi(int deminsion, double* mas_phi, double alpha)
{
    double* new_mas = new double[3];
    for (int i = 0; i < 3; i++)
    {
        new_mas[i] = mas_phi[i] - alpha * grad_J(i, mas_phi);
    }
    for (int i = 0; i < 3; i++)
    {
        mas_phi[i] = new_mas[i];
    }
}
// Функция J для метода деления пополам
double fast_J(int deminsion, double* masX, double alpha)
{
    double* mas_phi = new double[deminsion];
    phi(deminsion, mas_phi, alpha);
    return pow(mas_phi[0], 2) + pow(mas_phi[1], 2);
}
// Метод Половинного деления
double HalfDivision(int deminsion, double A, double B, double* masX, double EPS)
{
    double x1 = A, x2 = B, delta = EPS / 2; // Динамическая граница метода деления пополам
    do
    {
        if (fast_J(deminsion, masX, (x1 + x2 - delta) / 2) <= fast_J(deminsion, masX, (x1 + x2 + delta) / 2))
        {
            x2 = (x1 + x2 + delta) / 2;
        }
        else
        {
            x1 = (x1 + x2 - delta) / 2;
        }
    } while (abs(x2 - x1) > EPS);
    return (x1 + x2) / 2;
}
// Метод Наискорейшего спуска
void Fastes(int deminsion, double* masX, double eps)
{
    double alpha = HalfDivision(deminsion, 0, 100, masX, eps);
    double* masY = new double[3];
    int exit = 0;
    while (alpha >= eps * eps && exit != 5)
    {
        alpha = HalfDivision(deminsion, 0, 1000, masX, eps);
        for (int i = 0; i < deminsion; i++)
        {
            masY[i] = masX[i] - alpha * grad_J(i, masX);
        }
        if (J(masX) > J(masY))
        {
            for (int i = 0; i < 3; i++)
            {
                masX[i] = masY[i];
            }
            exit = 0;
        }
        else
        {
            exit++;
        }
    }
}
// Возвращаем значения равеств
double rav(int type_of_ne_rav, double* X, double q)
{
    switch (type_of_ne_rav)
    {
    case 0:
        return pow(0, 1);
    }
}
// Выбираем необходимое нам неравество в квадрате
double ne_rav(int type_of_ne_rav, double* X)
{
    switch (type_of_ne_rav)
    {
    case 0:
        return -X[0] + 1;
    case 1:
        return X[0] + X[1] - 2;
    }
}
// Поиск максимума у неравенства
double max_ne_rav(int type_of_ne_rav,double* X, double q)
{
    double znach = ne_rav(type_of_ne_rav, X);
    if (znach > 0) // Возвращаем максимальное значение или 0
    {
        return pow(znach, q);
    }
    else
    {
        return 0.;
    }
}
double Ak(int k) // Функция Ak
{
    return k;
}
// Метод штрафных функций
void penalty_func(int deminsion, int count_rav, int count_ne_rav, double q, double* masX, double eps, int& step)
{
    int k = 2;
    double Jk = 1000000;
    double Jk_2 = 0;
    while (abs(Jk - Jk_2) > eps)
    {
        if (k % 2 == 1) // Если у нас нету целого числа от k/2, то мы пропускаем его
        {
            k++;
        }
        else
        {
            Fastes(deminsion, masX, eps);
            cout << masX[0] << " " << masX[1] << endl;
            // Ищем значения J от шага k
            Jk = 0; // Обнуляем значения
            Jk += J(masX); // Добавляем само уравнение
            cout << Jk << " ";
            for (int i = 0; i < count_rav; i++) // Добавляем равенства
            {
                Jk += Ak(k) * rav(i, masX, q);
            }
            cout << Jk << " ";
            for (int i = 0; i < count_ne_rav; i++) // Добавляем неравенства
            {
                Jk += Ak(k) * max_ne_rav(i, masX, q);
            }
            cout << Jk << endl;
            // Ищем значения J от шага k/2
            Jk_2 = 0; // Обнуляем значения
            Jk_2 += J(masX); // Добавляем само уравнение
            cout << Jk_2 << " ";
            for (int i = 0; i < count_rav; i++) // Добавляем равенства
            {
                Jk_2 += Ak(k / 2) * rav(i, masX, q);
            }
            cout << Jk_2 << " ";
            for (int i = 0; i < count_ne_rav; i++) // Добавляем неравенства
            {
                Jk_2 += Ak(k / 2) * max_ne_rav(i, masX, q);
            }
            cout << Jk_2;
            separator();
            k++;
            step++;
        }
    }
}

int main()
{
    double q; // Константа для метода
    double eps = pow(10, -4);
    int deminsion = 2; // Размерность пространства
    int count_ne_rav = 2; // Количество уравнений типа неравенств
    int count_rav = 0; // Количество уравнений типа равенств
    cout << "This program use A(k)=k\n";
    separator();
    cout << "Enter the q=";
    cin >> q;
    q = abs(q); // На всякий случай нормируем
    if (q < 1)
    {
        cout << "You enter wrong q, set q = 2\n";
        q = 2;
    }
    separator();
    int step = 0; // Вывод количества итераций
    double* masX = new double[deminsion]; // Начальная точка и после точка для ответа
    for (int i = 0; i < deminsion; i++) // Узнаем точку для метода наискорейшего спуска
    {
        cout << "Choice x" << i + 1 << " = ";
        cin >> masX[i];
    }
    cout << endl;
    separator();
    penalty_func(deminsion, count_rav, count_ne_rav, q, masX, eps, step);
    cout << "Minimum point x=(";
    for (int i = 0; i < deminsion; i++) // Вывод ответа
    {
        if ( i != deminsion -1)
        {
            cout << masX[i] << " ; ";
        }
        else
        {
            cout << masX[i] << ")  J(x)=" << J(masX) << "  reached with " << step << " step\n";
        }
    }
    return 1;
}