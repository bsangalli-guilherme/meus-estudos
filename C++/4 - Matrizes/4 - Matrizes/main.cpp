#include <iostream>
#include <locale.h>
#include <random>
#include <vector>


using namespace std;

int main() {
	setlocale(LC_ALL, "");

	int n;
	cout << "Digite o tamanho da matriz (N): ";
	if(!(cin >> n) || n<=0) return 0;

	//Cria um vetor de vetores dinamico 
	vector<vector<int>> mat(n, vector<int>(n));

	random_device rd; //seed aleatória, coleta entropia real do hardware (se disponível) -> pode usar clock, movimentos do mouse, ruido elétrico, etc...
	mt19937 gen(rd());// mersenne twister, algoritmo que gera uma sequência deterministica de números aleatórios, baseado nesse caso na seed rd
	uniform_int_distribution<> dist(0, 9); //define como os números são gerados 

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			mat[i][j] = dist(gen);
			cout << mat[i][j] << " ";
		}
		cout << "\n";
	}

	int soma = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			soma += mat[i][i];
		}
		cout << "A soma da diagonal principal é: " << soma << endl;

		return 0;
	}
}


