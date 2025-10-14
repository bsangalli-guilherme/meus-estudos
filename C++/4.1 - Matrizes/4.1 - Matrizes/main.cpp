#include <iostream>
#include <vector>
#include <random>
#include <locale.h>
using namespace std;

int main() {
    setlocale(LC_ALL, "");

    int m, n;
    cout << "Digite o número de linhas: ";
    if (!(cin >> m) || m <= 0) return 0;

    cout << "Digite o número de colunas: ";
    if (!(cin >> n) || n <= 0) return 0;

    // Cria matriz m x n usando Array of arrays tendo tamanho variável por causa de vector
    vector<vector<int>> mat(m, vector<int>(n));

    // Gerador de números aleatórios
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dist(0, 9);

    // Preenche e imprime a matriz original
    cout << "\nMatriz original:\n";
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            mat[i][j] = dist(gen);
            cout << mat[i][j] << " ";
        }
        cout << endl;
    }

    // Imprime a matriz transposta diretamente
    cout << "\nMatriz transposta:\n";
    for (int i = 0; i < n; i++) {      
        for (int j = 0; j < m; j++) {  
            cout << mat[j][i] << " ";
        }
        cout << endl;
    }

    return 0;
}
