#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <stdbool.h>
#include <time.h>

typedef struct {
    int proximo_numero;
    int contador_primos;
    int numero_maximo;
    pthread_mutex_t mutex;
} DadosCompartilhados;

typedef struct {
    int id;
    DadosCompartilhados* compartilhado;
} DadosThread;

int ler_input_user_valido(const char* mensagem);
bool verificar_se_eh_primo(int numero);
int obter_proximo_numero(DadosCompartilhados* compartilhado);
void adicionar_primos_encontrados(DadosCompartilhados* compartilhado, int quantidade);
void* processar_numeros(void* arg);
double calcular_tempo_decorrido(clock_t inicio, clock_t fim);

int main() {
    int numero_maximo = ler_input_user_valido("Digite a quantidade de numeros a processar: ");
    int quantidade_threads = ler_input_user_valido("Digite a quantidade de threads a criar: ");

    DadosCompartilhados compartilhado = {
        .proximo_numero = 1,
        .contador_primos = 0,
        .numero_maximo = numero_maximo
    };

    pthread_mutex_init(&compartilhado.mutex, NULL);

    pthread_t* threads = (pthread_t*)malloc(quantidade_threads * sizeof(pthread_t));
    DadosThread* dados_threads = (DadosThread*)malloc(quantidade_threads * sizeof(DadosThread));

    if (threads == NULL || dados_threads == NULL) {
        printf("Erro ao alocar memoria\n");
        free(threads);
        free(dados_threads);
        pthread_mutex_destroy(&compartilhado.mutex);
        return 1;
    }

    
    clock_t inicio = clock();

    for (int i = 0; i < quantidade_threads; i++) {
        dados_threads[i].id = i;
        dados_threads[i].compartilhado = &compartilhado;
        pthread_create(&threads[i], NULL, processar_numeros, &dados_threads[i]);
    }

    for (int i = 0; i < quantidade_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    clock_t fim = clock();
    double tempo_decorrido = calcular_tempo_decorrido(inicio, fim);

    
    printf("Resultados\n\n");
    
    printf("Primos encontrados: %d\n", compartilhado.contador_primos);
    printf("Tempo de execucao: %.3f segundos\n", tempo_decorrido);
  

    pthread_mutex_destroy(&compartilhado.mutex);
    free(threads);
    free(dados_threads);

    return 0;
}

int ler_input_user_valido(const char* mensagem) {
    int valor;
    do {
        printf("%s", mensagem);
        if (scanf("%d", &valor) != 1) {
            while (getchar() != '\n');
            printf("Entrada invalida. Tente novamente.\n");
            continue;
        }
        if (valor <= 0) {
            printf("O valor deve ser positivo.\n");
        }
    } while (valor <= 0);
    return valor;
}

void* processar_numeros(void* arg) {
    DadosThread* dados = (DadosThread*)arg;
    int contador_local = 0;

    while (true) {
        int numero = obter_proximo_numero(dados->compartilhado);

        if (numero == -1) {
            break;
        }

        bool primo = verificar_se_eh_primo(numero);

        if (primo) {
            contador_local++;
            //printf("Thread %d: O numero %d e primo\n", dados->id, numero);
        }
        else {
            //printf("Thread %d: O numero %d nao e primo\n", dados->id, numero);
        }
    }

    adicionar_primos_encontrados(dados->compartilhado, contador_local);

    return NULL;
}

int obter_proximo_numero(DadosCompartilhados* compartilhado) {
    int numero_para_processar = -1;

    pthread_mutex_lock(&compartilhado->mutex);

    if (compartilhado->proximo_numero <= compartilhado->numero_maximo) {
        numero_para_processar = compartilhado->proximo_numero;
        compartilhado->proximo_numero++;
    }

    pthread_mutex_unlock(&compartilhado->mutex);

    return numero_para_processar;
}

bool verificar_se_eh_primo(int numero) {
    if (numero == 1) return false;
    if (numero == 2) return true;
    if (numero % 2 == 0) return false;

    for (int i = 3; i * i <= numero; i += 2) {
        if (numero % i == 0) return false;
    }
    return true;
}

void adicionar_primos_encontrados(DadosCompartilhados* compartilhado, int quantidade) {
    if (quantidade == 0) return;

    pthread_mutex_lock(&compartilhado->mutex);
    compartilhado->contador_primos += quantidade;
    pthread_mutex_unlock(&compartilhado->mutex);
}

double calcular_tempo_decorrido(clock_t inicio, clock_t fim) {
    return ((double)(fim - inicio)) / CLOCKS_PER_SEC;
}