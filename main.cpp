#include <iostream>
using namespace std;
void forwardEuler(double, double, int, double);
void backwardEuler(double, double, int, double);
void midpoint(double, double, int, double);
void Richardson(double, double, int, double);
void DirectIteration(double x0, double xk, int n, double y0);
void ImprovedEuler(double x0, double xk, int n, double y0);
void RungeKutty(double x0, double xk, int n, double y0);
int main() {
	cout << "Forward auler:" << endl;
	forwardEuler(0, 4, 4, 0);
	cout << "backward auler:" << endl;
	backwardEuler(0, 4, 4, 0);
	cout << "Midpoint:" << endl;
	midpoint(0, 4, 4, 0);	
	cout << "Richardson:" << endl;
	Richardson(0, 4, 4, 0);
	cout << "Direct iteration:" << endl;
	DirectIteration(0, 4, 4, 1);
	cout << "Improved Euler" << endl;
	ImprovedEuler(0, 1, 10, 1);
	cout << "Runge Kutty rz 4" << endl;
	RungeKutty(0, 4, 8, 0.5);
	system("pause");
}
double f(double x,double y) {
	return 2*x;
}
void forwardEuler(double x0,double xn,int n,double y0) {
	double h = abs(xn - x0) / (double)n;
	double xi =x0,yi =y0;
	for (int i = 0; i < n; i++) {
		cout << "x: " << xi << ", yi: " << yi << endl;
		yi = yi + h*(f(xi,yi));
		xi += h;
	}
	cout << "x: " << xi << ", yi: " << yi << endl;
	
}
void backwardEuler(double x0, double xn, int n, double y0) {
	double h = abs(xn - x0) / (double)n;
	double xi = x0, yi = y0;
	for (int i = 0; i < n; i++) {
		cout << "x: " << xi << ", yi: " << yi << endl;
		xi += h;
		yi = yi + h*(f(xi,yi));
	}
	cout << "x: " << xi << ", yi: " << yi << endl;

}
void midpoint(double x0, double xn, int n, double y0) {
	double h = abs(xn - x0) / (double)n;
	double xi = x0, yi = y0;
	for (int i = 0; i < n; i++) {
		cout << "x: " << xi << ", yi: " << yi << endl;		
		double xinext = x0 + (i+1)*h;
		yi = yi + h*((f(xi, yi) + f(xinext, yi))/2.0);
		xi += h;
	}
	cout << "x: " << xi << ", yi: " << yi << endl;
}
void Richardson(double x0, double xn, int n, double y0) {
	double xi = x0, yi = y0;
	auto ForwardEulerLambda = [](double x0, double xn, int n, double y0) -> double {
		double h = abs(xn - x0) / (double)n;
		double xi = x0, yi = y0;
		for (int i = 0; i < n; i++) {
			yi = yi + h*(f(xi, yi));
			xi += h;
		}
		return yi;
	};
	for (int i = 0; i < n; i++) {
		cout << "x: " << xi << ", yi: " << yi << endl;
		yi = (4 * ForwardEulerLambda(x0, i+1, (i+1)*2, y0) - ForwardEulerLambda(x0, i + 1, (i + 1), y0)) / 3.0;
		xi += 1;
	}
	cout << "x: " << xi << ", yi: " << yi << endl;
}
void DirectIteration(double x0, double xk, int n, double y0) {
	auto func = [](double x, double y) -> double {
		return x/y;
	};
	double yn = y0, xn = x0;
	double ynext,xnext;
	double dx = abs(x0 - xk) / n;
	double e = 0.1;
	for (int i = 0; i < n;i++) {
		xnext = x0 + (i+1)*dx;
		double pred = yn + dx*func(xnext, yn);
		double kor = yn + dx*func(xnext, pred);

		//cout << "pred: " << pred << ", kor: " << kor << endl;
		while (abs(kor - pred)/kor > e) {
			pred = kor;
			kor = yn + dx*func(xnext, pred);
			//cout << "pred: " << pred << ", kor: " << kor << endl;
		}
		yn = kor;
		cout << "x: " << (xnext) << ", yi: " << yn << endl;
	}
}
void ImprovedEuler(double x0, double xk, int n, double y0) {
	auto func = [](double x, double y) -> double {
		return 2*x - y;
	};
	double yn = y0, xn = x0, xnext;
	double dx = abs(x0 - xk) / n;
	for (int i = 0; i < n; i++) {
		double pred = yn + dx * func(xn, yn);
		 xnext = x0 + (i + 1)*dx;
		 double kor = yn + dx* func(xnext, pred);
		 yn = (pred + kor) / 2.0;
		 xn = xnext;
		 cout << "x: " << (xnext) << ", yi: " << yn << endl;
	}
}
void RungeKutty(double x0, double xk, int n, double y0) {
	auto func = [](double x, double y) -> double {
		return y - x*x + 1;
	};
	double yn = y0, xn = x0, dx = abs(x0-xk)/n;
	for (int i = 0; i < n;i++) {
		double k1, k2, k3, k4;
		k1 = dx * func(xn, yn);
		k2 = dx * func(xn + 0.5 *dx, yn + 0.5 *k1);
		k3 = dx * func(xn + 0.5 *dx, yn + 0.5 *k2);
		k4 = dx * func(xn + dx, yn + k3);
		yn = yn + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
		xn = x0 + (i + 1)*dx;
		cout << "x: " << xn << ", yi: " << yn << endl;
	}
}