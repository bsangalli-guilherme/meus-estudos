#include <iostream>
#include <locale.h>
#include "ehprimo.h"

using namespace std;

int main() {
	setlocale(LC_ALL, "");

	int numero;
	cout << "Digite um número inteiro positivo: ";
	cin >> numero;

	if(ehPrimo(numero)){
		cout << numero << " é primo!" << endl;
	}
	else {
		cout << numero << " não é primo!" << endl;
	}
}