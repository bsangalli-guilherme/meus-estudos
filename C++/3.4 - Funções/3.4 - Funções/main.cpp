#include <iostream>
#include <locale.h>
#include "funcs.h"
using namespace std;

int main() {
    setlocale(LC_ALL, "");
    int op;
    float a, b;

    do {
        menu();
        cin >> op;

        switch (op) {
        case 1:
            cout << "Digite dois números: ";
            cin >> a >> b;
            imprimirResultado(a, b, '+', soma(a, b));
            break;
        case 2:
            cout << "Digite dois números: ";
            cin >> a >> b;
            imprimirResultado(a, b, '-', subtracao(a, b));
            break;
        case 3:
            cout << "Digite dois números: ";
            cin >> a >> b;
            imprimirResultado(a, b, 'x', multiplicacao(a, b));
            break;
        case 4:
            cout << "Digite dois números: ";
            cin >> a >> b;
            
            imprimirResultado(a, b, '/', divisao(a, b));
            break;
        case 0:
            cout << "Saindo..." << endl;
            break;
        default:
            cout << "Opção inválida!" << endl;
        }
    } while (op != 0);

    return 0;
}
