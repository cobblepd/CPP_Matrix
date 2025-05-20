# CPP_Matrix
Реализация библиотеки для работы с матрицами на C++ с использованием ООП подхода
Описание проекта

Библиотека s21_matrix_oop предоставляет класс S21Matrix для работы с матрицами, реализующий основные матричные операции с использованием объектно-ориентированного подхода. Проект разработан на C++17 и включает все базовые операции линейной алгебры.
Основные возможности
Конструкторы и деструкторы

    S21Matrix() - конструктор по умолчанию

    S21Matrix(int rows, int cols) - конструктор с параметрами

    S21Matrix(const S21Matrix& other) - конструктор копирования

    S21Matrix(S21Matrix&& other) - конструктор перемещения

    ~S21Matrix() - деструктор

Базовые операции
cpp

bool EqMatrix(const S21Matrix& other);
void SumMatrix(const S21Matrix& other); 
void SubMatrix(const S21Matrix& other);
void MulNumber(const double num);
void MulMatrix(const S21Matrix& other);

Продвинутые операции
cpp

S21Matrix Transpose();
S21Matrix CalcComplements();
double Determinant();
S21Matrix InverseMatrix();

Перегруженные операторы
cpp

// Арифметические
S21Matrix operator+(const S21Matrix& other);
S21Matrix operator-(const S21Matrix& other);
S21Matrix operator*(const S21Matrix& other);
S21Matrix operator*(const double num);

// Присваивания
S21Matrix& operator=(const S21Matrix& other);
S21Matrix& operator+=(const S21Matrix& other);
S21Matrix& operator-=(const S21Matrix& other); 
S21Matrix& operator*=(const S21Matrix& other);
S21Matrix& operator*=(const double num);

// Сравнения
bool operator==(const S21Matrix& other);

// Индексация
double& operator()(int row, int col);

Требования к реализации

    Язык: C++17

    Компилятор: gcc

    Стиль кода: Google Style

    Покрытие тестами: GTest

    Приватные поля: matrix_, rows_, cols_

    Доступ к полям через accessors/mutators

Сборка и использование
Сборка библиотеки
bash

make s21_matrix_oop.a  # сборка статической библиотеки
make all               # сборка библиотеки и тестов

Запуск тестов
bash

make test  # запуск unit-тестов

Очистка
bash

make clean  # удаление временных файлов

Пример использования
cpp

#include "s21_matrix_oop.h"

int main() {
    S21Matrix a(3, 3);
    S21Matrix b(3, 3);
    
    // Заполнение матриц
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            a(i, j) = i + j;
            b(i, j) = i * j; 
        }
    }
    
    S21Matrix c = a + b;  // Сложение матриц
    S21Matrix d = a * b;  // Умножение матриц
    
    return 0;
}

Структура проекта

src/
├── s21_matrix_oop.cc  # реализация методов
├── s21_matrix_oop.h   # заголовочный файл
tests/                 # unit-тесты
Makefile               # файл сборки

Особенности реализации

    Полная поддержка move-семантики

    Обработка исключительных ситуаций

    Оптимизированные алгоритмы вычислений

    Подробное тестирование всех операций

    Соблюдение принципов ООП
