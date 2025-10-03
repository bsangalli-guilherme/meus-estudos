#include <iostream>
#include "funcs.h"
#include <locale.h>


using namespace std;


int main() {
	setlocale(LC_ALL, "");
	int op;
	float a = 0, b = 0;
	menu();
	cin >> op;
	switch (op) {
		case 1:
			cout << "\nInforme o primeiro número: ";
			cin >> a;
			cout << "\nInforme o segundo número: ";
			cin >> b;
			cout << a << " + " << b << " = " << soma(a, b) << endl;
			break;
		case 2:
			cout << "\nInforme o primeiro número: ";
			cin >> a;
			cout << "\nInforme o segundo número: ";
			cin >> b;
			cout << a << " - " << b << " = " << subtracao(a, b) << endl;
			break;
	}

}