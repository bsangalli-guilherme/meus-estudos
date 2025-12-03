#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define TAMANHO_SEGMENTO 32768
#define TAMANHO_LOTE_CONTAGEM 1000

typedef struct {
    int* primos_base;
    int total_primos_base;
} PrimosBase;

typedef struct {
    long long proximo_segmento_inicio;
    long long numero_maximo;
    PrimosBase* primos_base;
    long long contador_primos_total;
    pthread_mutex_t mutex;
} DadosSegmentados;

typedef struct {
    int id;
    DadosSegmentados* dados;
} DadosThread;

int ler_input_user_valido(const char* mensagem);
void crivo_simples_primos_base(int limite, PrimosBase* primos);
void* thread_crivo_segmentado(void* arg);
double calcular_tempo_decorrido(clock_t inicio, clock_t fim);

int main() {
    printf("=== CRIVO SEGMENTADO PARALELO COM HPC ===\n");
    printf("Otimizacoes: Cache L1 + Wheel Factorization\n\n");

    long long numero_maximo = ler_input_user_valido("Digite a quantidade de numeros a processar: ");
    int quantidade_threads = ler_input_user_valido("Digite a quantidade de threads a criar: ");

    clock_t inicio_total = clock();

    int limite_raiz = (int)sqrt(numero_maximo);

    printf("\nFase 1: Calculando primos base ate %d...\n", limite_raiz);
    PrimosBase primos_base;
    crivo_simples_primos_base(limite_raiz, &primos_base);
    printf("Encontrados %d primos base\n", primos_base.total_primos_base);

    clock_t fim_fase1 = clock();

    printf("\nFase 2: Processando segmentos em paralelo...\n");
    printf("Tamanho do segmento: %d bytes (otimizado para Cache L1)\n", TAMANHO_SEGMENTO);
    printf("Threads: %d\n\n", quantidade_threads);

    DadosSegmentados dados_segmentados = {
        .proximo_segmento_inicio = 0,
        .numero_maximo = numero_maximo,
        .primos_base = &primos_base,
        .contador_primos_total = 0
    };
    pthread_mutex_init(&dados_segmentados.mutex, NULL);

    pthread_t* threads = (pthread_t*)malloc(quantidade_threads * sizeof(pthread_t));
    DadosThread* dados_threads = (DadosThread*)malloc(quantidade_threads * sizeof(DadosThread));

    clock_t inicio_fase2 = clock();

    for (int i = 0; i < quantidade_threads; i++) {
        dados_threads[i].id = i;
        dados_threads[i].dados = &dados_segmentados;
        pthread_create(&threads[i], NULL, thread_crivo_segmentado, &dados_threads[i]);
    }

    for (int i = 0; i < quantidade_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    clock_t fim_total = clock();

    double tempo_fase1 = calcular_tempo_decorrido(inicio_total, fim_fase1);
    double tempo_fase2 = calcular_tempo_decorrido(inicio_fase2, fim_total);
    double tempo_total = calcular_tempo_decorrido(inicio_total, fim_total);

 
    printf("Resultados\n\n");

    
    printf("Primos encontrados:     %lld\n\n", dados_segmentados.contador_primos_total);
    
    printf("Tempo primos base:      %.4f seg\n", tempo_fase1);
    printf("Tempo segmentos:        %.4f seg\n", tempo_fase2);
    printf("Tempo total:            %.4f seg\n", tempo_total);

    pthread_mutex_destroy(&dados_segmentados.mutex);
    free(primos_base.primos_base);
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
            continue;
        }
    } while (valor <= 0);
    return valor;
}

void crivo_simples_primos_base(int limite, PrimosBase* primos) {
    bool* eh_primo = (bool*)malloc((limite + 1) * sizeof(bool));

    for (int i = 0; i <= limite; i++) {
        eh_primo[i] = true;
    }
    eh_primo[0] = eh_primo[1] = false;

    for (int i = 2; i * i <= limite; i++) {
        if (eh_primo[i]) {
            for (int j = i * i; j <= limite; j += i) {
                eh_primo[j] = false;
            }
        }
    }

    int count = 0;
    for (int i = 2; i <= limite; i++) {
        if (eh_primo[i]) count++;
    }

    primos->primos_base = (int*)malloc(count * sizeof(int));
    primos->total_primos_base = 0;

    for (int i = 2; i <= limite; i++) {
        if (eh_primo[i]) {
            primos->primos_base[primos->total_primos_base++] = i;
        }
    }

    free(eh_primo);
}

void* thread_crivo_segmentado(void* arg) {
    DadosThread* dados_thread = (DadosThread*)arg;
    DadosSegmentados* dados = dados_thread->dados;

    bool* segmento = (bool*)malloc(TAMANHO_SEGMENTO * sizeof(bool));
    long long contador_local = 0;

    while (true) {
        long long seg_inicio, seg_fim;

        pthread_mutex_lock(&dados->mutex);

        seg_inicio = dados->proximo_segmento_inicio;
        seg_fim = seg_inicio + TAMANHO_SEGMENTO - 1;

        if (seg_fim > dados->numero_maximo) {
            seg_fim = dados->numero_maximo;
        }

        if (seg_inicio > dados->numero_maximo) {
            pthread_mutex_unlock(&dados->mutex);
            break;
        }

        dados->proximo_segmento_inicio = seg_fim + 1;

        pthread_mutex_unlock(&dados->mutex);

        long long tamanho_seg = seg_fim - seg_inicio + 1;

        for (long long i = 0; i < tamanho_seg; i++) {
            segmento[i] = true;
        }

        for (int p = 0; p < dados->primos_base->total_primos_base; p++) {
            long long primo = dados->primos_base->primos_base[p];

            long long primeiro_multiplo;

            if (seg_inicio <= primo) {
                primeiro_multiplo = primo * primo;
            }
            else {
                primeiro_multiplo = ((seg_inicio + primo - 1) / primo) * primo;
            }

            for (long long j = primeiro_multiplo; j <= seg_fim; j += primo) {
                if (j >= seg_inicio) {
                    segmento[j - seg_inicio] = false;
                }
            }
        }

        for (long long i = 0; i < tamanho_seg; i++) {
            long long numero = seg_inicio + i;

            if (numero < 2) {
                continue;
            }

            if (numero == 2) {
                contador_local++;
                //printf("Thread %d: O numero %lld e primo\n", dados_thread->id, numero);
                continue;
            }

            if (numero % 2 == 0) {
               // printf("Thread %d: O numero %lld nao e primo\n", dados_thread->id, numero);
                continue;
            }

            if (segmento[i]) {
                contador_local++;
                //printf("Thread %d: O numero %lld e primo\n", dados_thread->id, numero);
            }
            else {
                //printf("Thread %d: O numero %lld nao e primo\n", dados_thread->id, numero);
            }
        }
    }

    pthread_mutex_lock(&dados->mutex);
    dados->contador_primos_total += contador_local;
    pthread_mutex_unlock(&dados->mutex);

    free(segmento);
    return NULL;
}

double calcular_tempo_decorrido(clock_t inicio, clock_t fim) {
    return ((double)(fim - inicio)) / CLOCKS_PER_SEC;
}