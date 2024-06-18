/*
#include <iostream>
#include <cmath>

using namespace std;

double f(double x, double y) {
    return x*x + 3*y*y + cos(x + y);
}

double minimize_x(double x, double y, double alpha, double tolerance) {
    auto grad = [&](double x) { return 2*x - sin(x + y); };
    double prev_x = x;

    while (true) {
        double grad_val = grad(prev_x);
        double new_x = prev_x - alpha * grad_val;

        if (abs(new_x - prev_x) < tolerance) {
            break;
        }

        prev_x = new_x;
    }

    return prev_x;
}

double minimize_y(double x, double y, double alpha, double tolerance) {
    auto grad = [&](double y) { return 6*y - sin(x + y); };
    double prev_y = y;

    while (true) {
        double grad_val = grad(prev_y);
        double new_y = prev_y - alpha * grad_val;

        if (abs(new_y - prev_y) < tolerance) {
            break;
        }

        prev_y = new_y;
    }

    return prev_y;
}

void coordinate_descent(double initial_x, double initial_y, double alpha, double tolerance, int max_iter) {
    double x = initial_x;
    double y = initial_y;

    for (int iter = 0; iter < max_iter; ++iter) {
        double new_x = minimize_x(x, y, alpha, tolerance);
        double new_y = minimize_y(new_x, y, alpha, tolerance);

        if (abs(f(new_x, new_y) - f(x, y)) < tolerance) {
            break;
        }

        x = new_x;
        y = new_y;
    }

    cout << "Минимум функции достигается в точке (x, y) = (" << x << ", " << y << ")" << endl;
    cout << "Значение функции в минимуме: " << f(x, y) << endl;
}

int main() {
    double initial_x = 1.0;
    double initial_y = 1.0;
    double alpha = 0.1;
    double tolerance = 1e-6;
    int max_iter = 1000;

    coordinate_descent(initial_x, initial_y, alpha, tolerance, max_iter);

    return 0;
}
*/

#include <iostream>
#include <cmath>

using namespace std;

double f(double x, double y) {
    return 7*x*x + 2*x*y + 5*y*y + x - 10*y;
}

void gradient(double x, double y, double& gx, double& gy) {
    gx = 14*x + 2*y + 1;
    gy = 2*x + 10*y - 10;
}

void gradient_descent(double x, double y, double alpha, double tolerance, int max_iter) {
    for (int iter = 0; iter < max_iter; ++iter) {
        double gx, gy;
        gradient(x, y, gx, gy);

        double nx = x - alpha * gx;
        double ny = y - alpha * gy;

        if (sqrt(gx * gx + gy * gy) < tolerance) {
            break;
        }

        x = nx;
        y = ny;
    }

    cout << "Минимум функции достигается в точке (x, y) = (" << x << ", " << y << ")" << endl;
    cout << "Значение функции в минимуме: " << f(x, y) << endl;
}

int main() {
    double x = 1.0;
    double y = 1.0;
    double alpha = 0.1;
    double tolerance = 1e-6;
    int max_iter = 1000;

    gradient_descent(x, y, alpha, tolerance, max_iter);

    return 0;
}
