#include <iostream>
#include <locale.h>
#include "dobro.h"
using namespace std;

int main() {
	setlocale(LC_ALL, "");
	int numero, dobrado;


	cout << "Digite um numero para saber seu dobro: ";
	cin >> numero;
	dobrado = dobro(numero);
	cout << "O dobro do numero " << numero << " é " << dobrado;



}
