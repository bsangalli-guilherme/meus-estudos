#include "funcs.h"
#include <iostream>
#include <cmath>
using namespace std;

void menu() {
    cout << "\nEscolha uma das opções abaixo\n";
    cout << "(1) Soma\n(2) Subtração\n(3) Multiplicação\n(4) Divisão\n(0) Sair\n";
}

float soma(float a, float b) { return a + b; }
float subtracao(float a, float b) { return a - b; }
float multiplicacao(float a, float b) { return a * b; }

float divisao(float a, float b) {
    if (b == 0) {
        return NAN;
    }
    return a / b;
}


void imprimirResultado(float a, float b, char op, float resultado) {
    if (op == '/' && isnan(resultado)) {
        cout << "Impossível dividir por 0 " << endl;
    }
    else {
        cout << a << " " << op << " " << b << " = " << resultado << endl;
    }
}
