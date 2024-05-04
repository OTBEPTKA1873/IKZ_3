#include <iostream>
#include <cmath>
using namespace std;

int dimension = 2; // Размерность пространства
int count_nerav = 2; // Количество уравнений типа неравенств
int count_rav = 0; // Количество уравнений типа равенств
// Отрез кодов друг от друга
void separator()
{
    cout << "\n==================================================================================\n\n";
}
// Функция для перехода от P(x) к Pk(x)
double Ak(int k)
{
    return k;
}

// Функция, которую минимизируем
double J(double* masX)
{
    return pow(masX[0], 2) + pow(masX[1], 2);
}
// Градиент функции
double grad_J(int grad, double* masX)
{
    switch (grad)
    {
    case 0:
        return 2 * masX[0];
    case 1:
        return 2 * masX[1];
    }
}

// Функции равенств
double rav(int rav, double* masX, double q, int k)
{
    switch (rav)
    {
    case 0:
        return 0;
    }
}
// Градиенты равенств
double grad_rav(int rav, int grad, double* masX, double q, int k)
{
    switch (rav)
    {
    case 0:
        switch (grad)
        {
        case 0:
            return 0;
        }
    }
}

// Функции неравенств
double nerav(int rav, double* masX, double q, int k)
{
    switch (rav)
    {
    case 0:
        if (-masX[0] + 1 > 0)
        {
            return Ak(k) * pow(-masX[0] + 1, q);
        }
        else
        {
            return 0;
        }
    case 1:
        if (masX[0] + masX[1] - 2 > 0)
        {
            return Ak(k) * pow(masX[0] + masX[1] - 2, q);
        }
        else
        {
            return 0;
        }
    }
}
// Градиенты неравенств
double grad_nerav(int rav, int grad, double* masX, double q, int k)
{
    switch (rav)
    {
    case 0:
        switch (grad)
        {
        case 0:
            if (-masX[0] + 1 > 0)
            {
                return Ak(k) * (-q) * pow(-masX[0] + 1, q - 1);
            }
            else
            {
                return 0;
            }
        case 1:
            return 0;
        }
    case 1:
        switch (grad)
        {
        case 0:
            if (masX[0] + masX[1] - 2 > 0)
            {
                return Ak(k) * q * pow(masX[0] + masX[1] - 2, q - 1);
            }
            else
            {
                return 0;
            }
        case 1:
            if (masX[0] + masX[1] - 2 > 0)
            {
                return Ak(k) * q * pow(masX[0] + masX[1] - 2, q - 1);
            }
            else
            {
                return 0;
            }
        }
    }
}

// Ищем Фk, который будем искать
double PH(double* masX, double q, int k)
{
    double Ph = J(masX); // Даем значение J
    for (int j = 0; j < count_rav; j++) // Добавляем равенства
    {
        Ph += Ak(k) * rav(j, masX, q, k);
    }
    for (int j = 0; j < count_nerav; j++) // Добавляем неравенства
    {
        Ph += Ak(k) * nerav(j, masX, q, k); // Здесь уже есть проверка на max
    }
    return Ph;
}
// Найдем производную Фк для метода наискорейшего спуска
double grad_PH(int grad, double* masX, double q, int k)
{
    double grad_Ph = grad_J(grad, masX); // Даем значение J
    for (int j = 0; j < count_rav; j++) // Добавляем равенства
    {
        grad_Ph += Ak(k) * grad_rav(j, grad, masX, q, k);
    }
    for (int j = 0; j < count_nerav; j++) // Добавляем неравенства
    {
        grad_Ph += Ak(k) * grad_nerav(j, grad, masX, q, k);
    }
    return grad_Ph;
}

// Облегчение написания функции фи
void phi(double* mas_phi, double alpha, double q, int k)
{
    double* new_mas = new double[dimension];
    double* grad_Phi = new double[dimension];
    for (int i = 0; i < dimension; i++)
    {
        grad_Phi[i] = grad_PH(i, mas_phi, q, k);
    }
    for (int i = 0; i < dimension; i++)
    {
        new_mas[i] = mas_phi[i] - alpha * grad_Phi[i];
    }
    for (int i = 0; i < dimension; i++)
    {
        mas_phi[i] = new_mas[i];
    }
}
// Функция J для метода деления пополам
double fast_J(double* masX, double alpha, double q, int k)
{
    double* mas_phi = new double[dimension];
    for (int i = 0; i < dimension; i++)
    {
        mas_phi[i] = masX[i]; // Передаем точку
    }
    phi(mas_phi, alpha, q, k);
    return PH(mas_phi, q, k);
}
// Метод Половинного деления
double HalfDivision(double A, double B, double* masX, double EPS, double q, int k)
{
    double x1 = A, x2 = B, delta = EPS / 2; // Динамическая граница метода деления пополам
    do
    {
        if (fast_J(masX, (x1 + x2 - delta) / 2, q, k) <= fast_J(masX, (x1 + x2 + delta) / 2, q, k))
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
void Fastes(double* masX, double eps, double q, int k)
{
    double alpha = HalfDivision(0, 100, masX, eps, q, k);
    double* masY = new double[dimension];
    int exit = 0;
    while (alpha >= eps * eps && exit != 5)
    {
        alpha = HalfDivision(0, 1000, masX, eps, q, k);
        for (int i = 0; i < dimension; i++)
        {
            masY[i] = masX[i] - alpha * grad_PH(i, masX, q, k);
        }
        if (PH(masX, q, k) > PH(masY, q, k))
        {
            for (int i = 0; i < dimension; i++)
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
// Метод штрафных функций
void penalty_func(double* masX, int& step, double eps, double q)
{
    int k = 1;
    int check = 0;
    double* masXk = new double[dimension];
    double* masXk_2 = new double[dimension];
    cout << "Start method\n\n";
    while (check != dimension)
    {
        check = 0;
        if (k % 2 == 1)
        {
            k++;
        }
        else
        {
            for (int i = 0; i < dimension; i++)
            {
                masXk[i] = masX[i];
                masXk_2[i] = masX[i];
            }
            Fastes(masXk, eps, q, k);
            Fastes(masXk_2, eps, q, k / 2);
            if (step % 10 ==0)
            {
                cout << k << " " << masXk[0] << " " << masXk[1] << " " << J(masXk) << "\n";
                cout << k << " " << masXk_2[0] << " " << masXk_2[1] << " " << J(masXk_2);
                separator();
            }
            k++;
            step++;
            for (int i = 0; i < dimension; i++)
            {
                if (abs(masXk[i] - masXk_2[i]) < eps * 10)
                {
                    check++;
                }
            }
        }
    }
    for (int i = 0; i < dimension; i++)
    {
        masX[i] = masXk[i];
    }
}

int main()
{
    double q = 2; // Константа для метода
    double eps = pow(10, -4);
    cout << "This program use A(k)=k and q=2";
    separator();
    int step = 0; // Вывод количества итераций
    double* masX = new double[dimension]; // Начальная точка и после точка для ответа
    for (int i = 0; i < dimension; i++) // Узнаем точку для метода наискорейшего спуска
    {
        cout << "Choice x" << i + 1 << " = ";
        cin >> masX[i];
    }
    cout << endl;
    separator();
    penalty_func(masX, step, eps, q);
    cout << masX[0] << " " << masX[1] << " " << J(masX) << " with step=" << step;
    return 1;
}