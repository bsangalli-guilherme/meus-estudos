#include <iostream>
#include <locale.h>
#include "fatorial.h"

using namespace std;

int main() {
	setlocale(LC_ALL, "");
	int n;
	
	cout << "Informe um numero inteiro positivo para calcular sua fatorial: ";
	cin >> n;
	int resultado = fatorial(n);
	if (resultado == -1) {
		cout << "Não existe fatorial de número negativo!" << resultado << endl;

	}
	else {
		cout << "O fatorial de " << n << " é " << resultado << endl;
	}
	return 0;
}