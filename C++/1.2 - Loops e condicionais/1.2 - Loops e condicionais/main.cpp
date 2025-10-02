/*Escreva um programa em C++ que faça o seguinte:

O programa deve pedir ao usuário para digitar vários números inteiros positivos.

A leitura deve continuar até o usuário digitar o número 0 (esse valor não conta na soma nem no cálculo).

Depois disso, o programa deve mostrar:

A quantidade de números digitados (sem contar o 0).

A soma de todos os números digitados.

A média aritmética desses números.

Se o usuário não digitar nenhum número antes do 0, mostre uma mensagem avisando que nenhum valor válido foi inserido.

*/

#include <iostream>
#include <locale.h>
using namespace std;

int main() {
	setlocale(LC_ALL, "");//usa local padrão do PC

	int i = 0, n = -1, soma = 0;
	float media = 0;

	while (n != 0) {
		cout << "Digite um numero (para sair digite 0): ";
		cin >> n;
		
		if (n != 0) {
			i++;
			soma += n;
		}
		
	}
	if (i == 0) {
		cout << "Nenhum valor foi inserido." << endl;
	}else{
		media = (float)soma / i;
		cout << "Foram digitados" << i << " números, sua soma foi de " << soma << ", e a media foi " << media << endl;
	}
	return 0;
}