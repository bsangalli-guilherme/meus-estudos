#include <iostream>
#include <string>
#include <locale.h>

typedef struct {
	std:: string nome;
	int idade;
	std::string CPF;
}Pessoa;

int main() {
	setlocale(LC_ALL, "");
	Pessoa pessoas[10];
	for (int i = 0; i < 10; i++) {
		std::cout << "Digite o nome da " << i + 1<< "ª pessoa "  << ": ";
		std::getline(std::cin, pessoas[i].nome);

		std::cout << "Digite a idade da " << i + 1<< "ª pessoa "  << ": ";
		std::cin >> pessoas[i].idade;

		std::cout << "Digite o CPF da " << i + 1 << "ª pessoa " << ": ";
		std::cin >> pessoas[i].CPF;

		std::cin.ignore();

	}

	std::cout << "\n--- Lista de Pessoas ---\n";
	for (int i = 0; i < 3; i++) {
		std::cout << i + 1 << "ª pessoa " <<  ": " << pessoas[i].nome
			<< ", " << pessoas[i].idade << " anos, CPF " << pessoas[i].CPF << "\n";
	}

	return 0;
}