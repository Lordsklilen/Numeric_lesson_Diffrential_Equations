#include <iostream>
using namespace std;
void forwardEuler(double, double, int, double);
void backwardEuler(double, double, int, double);
void midpoint(double, double, int, double);
void Richardson(double,double,int ,double );
double forwardEulerRichardson(double x0, double xn, int n, double y0);
int main() {
	cout << "Forward auler:" << endl;
	forwardEuler(0, 4, 4, 0);
	cout << "backward auler:" << endl;
	backwardEuler(0, 4, 4, 0);
	cout << "Midpoint:" << endl;
	midpoint(0, 4, 4, 0);	
	cout << "Richardson:" << endl;
	Richardson(0, 4, 4, 0);

	
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
double forwardEulerRichardson(double x0, double xn, int n, double y0) {
	double h = abs(xn - x0) / (double)n;
	double xi = x0, yi = y0;
	for (int i = 0; i < n; i++) {
		yi = yi + h*(f(xi, yi));
		xi += h;
	}
	return yi;
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
	for (int i = 0; i < n; i++) {
		cout << "x: " << xi << ", yi: " << yi << endl;
		yi = (4 * forwardEulerRichardson(x0, i+1, (i+1)*2, y0) - forwardEulerRichardson(x0, i + 1, (i + 1), y0)) / 3.0;
		xi += 1;
	}
	cout << "x: " << xi << ", yi: " << yi << endl;
}