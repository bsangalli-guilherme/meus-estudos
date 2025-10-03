#include "funcs.h"
#include <iostream>

void menu() {
	int op;
	std :: cout << "\nEscolha umas das opções abaixo\n";
	std::cout << "(1) Soma \n(2) Subtração \n(0) Sair\n" << std::endl;
}

float soma(float a, float b) {
	return a + b;
}

float subtracao(float a, float b) {
	return a - b;
}