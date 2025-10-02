#include <iostream>
#include <locale.h>
using namespace std;

int main() {
    setlocale(LC_ALL, "");

    float vet[5];   // vetor de notas
    float soma = 0; // soma para calcular a média
    float maior, menor;

    // Leitura da primeira nota (separada para inicializar maior/menor)
    cout << "Informe a nota do 1º aluno (0 a 10): ";
    cin >> vet[0];
    while (vet[0] < 0 || vet[0] > 10) {
        cout << "Nota inválida! Informe novamente (0 a 10): ";
        cin >> vet[0];
    }

    soma = vet[0];
    maior = vet[0];
    menor = vet[0];

    // Leitura das próximas 4 notas
    for (int i = 1; i < 5; i++) {
        cout << "Informe a nota do " << i + 1 << "º aluno (0 a 10): ";
        cin >> vet[i];

        while (vet[i] < 0 || vet[i] > 10) {
            cout << "Nota inválida! Informe novamente (0 a 10): ";
            cin >> vet[i];
        }

        soma += vet[i];

        if (vet[i] > maior) maior = vet[i];
        if (vet[i] < menor) menor = vet[i];
    }

    float media = soma / 5;

    // Exibição dos resultados
    cout << "\n===== Resultados =====\n";
    for (int i = 0; i < 5; i++) {
        cout << "Nota do " << i + 1 << "º aluno: " << vet[i];
        if (vet[i] >= media) {
            cout << " -> Acima da média\n";
        }
        else {
            cout << " -> Abaixo da média\n";
        }
    }

    cout << "\nMédia da turma: " << media << endl;
    cout << "Maior nota: " << maior << endl;
    cout << "Menor nota: " << menor << endl;

    return 0;
}
