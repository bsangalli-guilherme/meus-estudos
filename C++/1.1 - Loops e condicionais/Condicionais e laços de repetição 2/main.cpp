#include <iostream> // biblioteca de entrada de dados 
#include <locale.h> // configuração de localidade

int main() {
	setlocale(LC_ALL, "");//usa local padrão do PC
	
	int n=0, soma = 0, i = 1; //declaração do valor
	std::cout << "Digite um número inteiro positivo: ";
	std::cin >> n;



	if (n > 0) {
		
		for (i; i <= n; i++) {
			soma += i;
		}
		std::cout << "A soma de 1 até " << n << " é de " << soma << std::endl;
	}
	return 0;
}