#include <iostream>
#include <locale.h>
#include "fatorial.h"

using namespace std;

int main() {
	int n;
	
	cout << "Informe um numero nteiro positivo oara calcular sua fatorial: ";
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