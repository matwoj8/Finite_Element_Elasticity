#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <fstream>
#include <utility>
#include <sstream>
#include <algorithm>
#include <string>

using namespace std;

string to_comma_format(double number) {
    ostringstream oss;
    oss << fixed << number;
    string result = oss.str();
    replace(result.begin(), result.end(), '.', ',');
    return result;
}

double E(double x) {
    return x <= 1.0 ? 3.0 : 5.0;
}

double e(int n, int i, double x) {
    double h = 2.0 / n;
    return max(0.0, 1.0 - abs(x / h - i));
}

double e_prim(int n, int i, double x) {
    double h = 2.0 / n;
    if (x <= (i - 1) * h || x >= (i + 1) * h) return 0.0;
    return x <= i * h ? 1.0 / h : -1.0 / h;
}

double calculate_integral(int n, int i, int j, function<double(double)> integrand) {
    if (abs(j - i) > 1) return 0.0;
    double h = 2.0 / n, start = max(0.0, h * (i - 1)), end = min(2.0, h * (i + 1));
    double step = (end - start) / 1000;Z
    double result = 0.0;
    for (double x = start; x <= end; x += step) {
        result += integrand(x) * step;
    }
    return result;
}

pair<vector<vector<double>>, vector<double>> fill(int n) {
    vector<vector<double>> B(n, vector<double>(n, 0.0));
    vector<double> L(n, 0.0);
    L[0] = -30.0 * e(n, 0, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            B[i][j] = -3.0 * e(n, i, 0.0) * e(n, j, 0.0) +
                calculate_integral(n, i, j, [&](double x) { return E(x) * e_prim(n, i, x) * e_prim(n, j, x); });
        }
    }
    return { B, L };
}

void show_plot(const vector<double>& solution, int n) {
    ofstream file("elastic_deformation_plot.csv");
    file << "x;u(x)\n";
    double step = 2.0 / n;
    for (int i = 0; i <= n; ++i) {
        file << to_comma_format(i * step) << ";" << to_comma_format(solution[i]) << "\n";
    }
    file.close();
    cout << "Results saved to 'elastic_deformation_plot.csv'.\n";
}

int main() {
    int n;
    cout << "Input n: ";
    cin >> n;

    auto matrices = fill(n); // Tworzenie macierzy B i wektora L
    vector<vector<double>> B = matrices.first;
    vector<double> L = matrices.second;
    vector<double> solution(n + 1, 0.0);

    // Rozwi¹zywanie uk³adu równañ (eliminacja Gaussa)
    for (int i = 0; i < n; ++i) {
        for (int k = i + 1; k < n; ++k) {
            double factor = B[k][i] / B[i][i];
            for (int j = 0; j < n; ++j) {
                B[k][j] -= factor * B[i][j];
            }
            L[k] -= factor * L[i];
        }
    }
    for (int i = n - 1; i >= 0; --i) {
        solution[i] = L[i];
        for (int j = i + 1; j < n; ++j) {
            solution[i] -= B[i][j] * solution[j];
        }
        solution[i] /= B[i][i];
    }

    show_plot(solution, n);
    return 0;
}
