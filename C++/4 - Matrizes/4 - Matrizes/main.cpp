#include <iostream>
#include <locale.h>

using namespace std;

int main() {
	setlocale(LC_ALL, "");

	int mat[3][3], soma = 0;

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			cout << "Elemento [" << i << "][" << j << "]: ";
			cin >> mat[i][j];
		}
	}

	
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			if (i == j) {
				soma += mat[i][j];
			}
		}
		cout << endl;
	}
	cout << "A soma da diagonal principal é: " << soma << endl;

	return 0;

}


