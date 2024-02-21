//
//  main.cpp
//  Gaussian Elimination
//
//  Created by Behrouz on 10/6/15.
//  Copyright (c) 2015 Behrouz. All rights reserved.
//

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

// Function to print the augmented matrix A
void printMatrix(const vector<vector<double>>& A) {
    int n = (int)A.size(); // Number of rows in matrix A
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n + 1; ++j) {
            cout << A[i][j] << "\t";
            if (j == n - 1) {
                // Print a vertical bar to separate the coefficients from the constants
                cout << "| ";
            }
        }
        cout << "\n"; // Move to the next line after printing each row
    }
    cout << endl; // Add an extra line for better separation
}

// Function to perform Gaussian elimination and solve the system of equations
vector<double> gaussianElimination(vector<vector<double>> A) {
    int n = (int)A.size(); // Number of rows (and variables)

    // Forward elimination to form an upper triangular matrix
    for (int i = 0; i < n; ++i) {
        // Search for maximum element in this column to reduce numerical error
        double maxEl = abs(A[i][i]);
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (abs(A[k][i]) > maxEl) {
                maxEl = abs(A[k][i]);
                maxRow = k;
            }
        }

        // Swap the maximum row with the current row
        swap(A[maxRow], A[i]);

        // Eliminate all rows below the current row
        for (int k = i + 1; k < n; ++k) {
            double coeff = -A[k][i] / A[i][i];
            for (int j = i; j <= n; ++j) {
                if (j == i) {
                    A[k][j] = 0;
                } else {
                    A[k][j] += coeff * A[i][j];
                }
            }
        }
    }

    // Back substitution to solve for variables
    vector<double> solution(n);
    for (int i = n - 1; i >= 0; --i) {
        solution[i] = A[i][n] / A[i][i];
        for (int k = i - 1; k >= 0; --k) {
            A[k][n] -= A[k][i] * solution[i];
        }
    }

    return solution;
}

int main() {
    int n; // Number of variables (and equations)
    cout << "Enter the number of variables: ";
    cin >> n;

    vector<vector<double>> A(n, vector<double>(n + 1, 0)); // Augmented matrix (n x n+1)

    cout << "Enter the coefficients of the equations, row by row:\n";
    for (int i = 0; i < n; ++i) {
        cout << "Row " << i + 1 << ": ";
        for (int j = 0; j < n; ++j) {
            cin >> A[i][j];
        }
    }

    cout << "Enter the constants (right-hand side values):\n";
    for (int i = 0; i < n; ++i) {
        cin >> A[i][n];
    }

    // Print the input matrix
    cout << "Augmented Matrix:\n";
    printMatrix(A);

    // Perform Gaussian elimination and print the solution
    vector<double> solution = gaussianElimination(A);
    cout << "Solution:\n";
    for (int i = 0; i < n; ++i) {
        cout << "Variable " << i + 1 << ": " << solution[i] << "\n";
    }

    return 0;
}
