#include <iostream>
#include <locale.h> // Para configurar localidade

using namespace std;

int main() {
    setlocale(LC_ALL, ""); // Usa a localidade padrão do sistema (Português, se configurado no Windows)

    int N;
    cout << "Digite um número inteiro positivo: ";
    cin >> N;

    int i = 1;
    while (i <= N) {
        if (i % 2 == 0) {
            cout << i << " é par" << endl;
        }
        else {
            cout << i << " é ímpar" << endl;
        }
        i++;
    }

    return 0;
}
