#include<math.h> // pow
#include<iostream>
#include<typeinfo> // typeid
#include <random>
#include "Matrix.h" // Matrix class
#include <map>
using namespace std;


template <typename T>
Matrix<T> compute_grad(Matrix<T>& A, Matrix<T>& B, Matrix<T>& C) {
    Matrix<T> _A = A;
    Matrix<T> _B = B;
    int numIterations = 100000;
    double alpha = 0.0005;
    int m = _A.shape[0];
    Matrix<T> gradient;
    Matrix<T> hypothesis;
    Matrix<T> loss;
    Matrix<T> AT(transpose(_A));
    Matrix<T> theta = C;

    for (int i = 0; i < numIterations; i++) {
        hypothesis = _A * theta;
        loss = hypothesis - _B;
        gradient = AT * loss / m;
        theta = theta - (alpha * gradient);
    }
    return theta;
}

double drand48()
{
    return rand() / (RAND_MAX + 1.);
}

void main()
{
    Matrix<double> A(100, 2);
    Matrix<double> B(100, 1);
    Matrix<double> C(A.shape[1], 1);
    C = C.fill(1);

    double bias = 25;
    double variance = 10;

    for (int i = 0; i < 100; i++) {
        A(i, 0) = 1;
        A(i, 1) = i;

        B(i, 0) = (i + bias) + drand48() * variance;


    }
    C = compute_grad(A, B, C);
    cout << "A\n----------------" << endl;
    print(A);
    cout << "B\n----------------" << endl;
    print(B);
    cout << "C\n----------------" << endl;
    print(C);
}