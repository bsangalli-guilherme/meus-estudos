#include <iostream>
#include <locale.h>

using namespace std;

int main() {
    setlocale(LC_ALL, "");

    int num[10]; // vetor para armazenar até 10 números
    int i = 0;   // índice do vetor

    cout << "Digite números inteiros não negativos (máx 10 números, 0 para sair):\n";

    // Lê o primeiro número
    while (true) {
        cout << "Digite o 1º número: ";
        cin >> num[0];

        if (num[0] < 0) {
            cout << "Número inválido! Digite um número NÃO NEGATIVO.\n";
            continue;
        }
        break;
    }

    // Se o usuário digitar 0 no primeiro número, encerra
    if (num[0] == 0) {
        cout << "Nenhum número válido foi digitado.\n";
        return 0;
    }

    int maior = num[0]; // inicializa maior com o primeiro número
    int menor = num[0]; // inicializa menor com o primeiro número
    i = 1;              // já temos um número válido

    // Loop para ler os próximos números
    for (; i < 10; i++) {
        cout << "Digite o " << i + 1 << "º número: ";
        cin >> num[i];

        while (num[i] < 0) { // validação de negativo
            cout << "Número inválido! Digite um número NÃO NEGATIVO: ";
            cin >> num[i];
        }

        if (num[i] == 0) { // encerra se digitar 0
            break;
        }

        if (num[i] > maior) maior = num[i]; // atualiza maior
        if (num[i] < menor) menor = num[i]; // atualiza menor
    }

    // Número total de elementos válidos digitados
    int total = i;

    // Mostra os números na ordem inversa
    cout << "\nNúmeros na ordem inversa:\n";
    for (int j = total - 1; j >= 0; j--) {
        cout << num[j] << " ";
    }
    cout << endl;

    // Mostra maior e menor
    cout << "Maior número: " << maior << endl;
    cout << "Menor número: " << menor << endl;

    return 0;
}