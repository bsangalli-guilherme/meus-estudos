#include <iostream> // biblioteca de entrada de dados 
#include <locale.h> // configuração de localidade
using namespace std;//faz com eu que não precise usar std::cout e std::cin -> uso só cin e cout

int main() {
	setlocale(LC_ALL, "");//usa local padrão do PC
	
	int n=0, soma = 0, i = 1; //declaração do valor
	cout << "Digite um número inteiro positivo: ";
	cin >> n;


	if (n > 0) {
		
		for (i; i <= n; i++) {
			soma += i;
		}
		cout << "A soma de 1 até " << n << " é de " << soma << endl;
	}
	return 0;
}