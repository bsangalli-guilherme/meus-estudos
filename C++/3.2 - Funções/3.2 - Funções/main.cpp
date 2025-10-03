#include <iostream>
#include "media.h"
#include <iomanip>
#include <locale.h>


using namespace std;

int main() {
	setlocale(LC_ALL, "");
	const int TAM = 5;
	float vet[TAM];
	int i = 0;

	for (; i < TAM; i++) {
		cout << "Informe a " << i + 1 << "º primeira nota: ";
		cin >> vet[i];
	}
	cout << fixed << setprecision(2);
	cout << "A media da turma foi " << media(vet, TAM) << endl;



}