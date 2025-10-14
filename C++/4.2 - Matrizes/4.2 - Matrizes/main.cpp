#include <iostream>
#include <locale.h>
#include <random>
#include <vector>

using namespace std;

int main() {
	setlocale(LC_ALL, "");
	int m = 0, n = 0;

	cout << "Digite o número de linhas da matriz: ";
	if(!(cin >> m)|| m<=0) return 0;
	cout << "Digite o número de colunas da matriz: ";
	if(!(cin >> n)|| n<=0) return 0;

	vector<vector<int>> mat(m, vector<int>(n));
	vector<int>soma_linha(m);
	vector<int>soma_col(n);

	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<> dist(0, 9);

	cout << "Matriz preenchida: \n";
	for (int i = 0; i < m; i++) {
		
		for (int j = 0; j < n; j++) {
			mat[i][j] = dist(gen);
			cout << mat[i][j] << "\t";
			soma_linha[i] += mat[i][j];
			soma_col[j] += mat[i][j];
		}
		cout << endl;
	}
	for (int i = 0; i < m; i++) {
		cout << "\nA soma da linha " << i + 1 << " é: "  << soma_linha[i];
	}
	for (int i = 0; i < n; i++) {
		cout << "\nA soma da coluna " << i + 1 << " é: " << soma_col[i];
	}


}

/*

Dada uma matriz NxM, calcule e mostre a soma de cada linha e de cada coluna.
Dica: use dois loops — um para linhas, outro para colunas, ou acumule em vetores somaLinha[N], somaColuna[M].
 
*/