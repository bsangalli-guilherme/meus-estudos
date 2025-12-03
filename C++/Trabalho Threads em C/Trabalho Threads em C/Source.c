#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>

// Tamanho do bloco que cada thread processa na fase de contagem e impressão
#define TAMANHO_BLOCO_LEITURA 1000


typedef struct {
    int inicio;
    int fim;
} Intervalo;

//struct para a fase 1 - marcação dos múltiplos
typedef struct {
    int proximo_primo_semente;
    int limite_raiz_quadrada;
    int limite_superior_global;
    bool* mapa_de_primos;
    pthread_mutex_t mutex_controle;
} ContextoMarcacao;

//identificadores para as threads na fase 1
typedef struct {
    int id_thread;
    ContextoMarcacao* contexto;
} WorkerMarcacao;

//struct para a fase 2 - verificação se é primo(true) no mapa_primos e contagem da quantidade de primos
typedef struct {
    int proximo_indice_leitura; // Cursor para distribuir blocos de leitura
    int total_primos_encontrados;
    int limite_superior_global;
    bool* mapa_de_primos;       // Vetor já marcado (somente leitura agora)
    pthread_mutex_t mutex_controle;
} ContextoContagem;

//identificadores para as threads na fase 2
typedef struct {
    int id_thread;
    ContextoContagem* contexto;
} WorkerContagem;

// Protótipos
int ler_inteiro_positivo(const char* mensagem);
double obter_tempo_decorrido(clock_t inicio, clock_t fim);

//fase 1
void orquestrar_fase_marcacao(bool* mapa_primos, int limite_superior, int num_threads);
void* executar_tarefa_marcacao(void* arg);

//fase 2
void* executar_tarefa_contagem_impressao(void* arg);
Intervalo reservar_proximo_intervalo(ContextoContagem* contexto);


int main() {

    int limite_superior = ler_inteiro_positivo("Digite a quantidade de numeros a processar: ");
    int num_threads = ler_inteiro_positivo("Digite a quantidade de threads a criar: ");


    bool* mapa_primos = (bool*)malloc((limite_superior + 1) * sizeof(bool));
    if (mapa_primos == NULL) {
        printf("Erro Fatal: Memoria insuficiente para alocar o vetor de %d posicoes.\n", limite_superior);
        return 1;
    }



    clock_t inicio_total = clock();

    orquestrar_fase_marcacao(mapa_primos, limite_superior, num_threads);

    clock_t fim_fase1 = clock();

    ContextoContagem contexto_contagem = {
        .proximo_indice_leitura = 1,
        .total_primos_encontrados = 0,
        .limite_superior_global = limite_superior,
        .mapa_de_primos = mapa_primos
    };
    pthread_mutex_init(&contexto_contagem.mutex_controle, NULL);

    pthread_t* threads = (pthread_t*)malloc(num_threads * sizeof(pthread_t));
    WorkerContagem* workers_contagem = (WorkerContagem*)malloc(num_threads * sizeof(WorkerContagem));


    if (threads == NULL || workers_contagem == NULL) {
        printf("Erro Fatal: Memoria insuficiente para threads da fase 2.\n");
        free(threads);
        free(workers_contagem);
        free(mapa_primos);
        pthread_mutex_destroy(&contexto_contagem.mutex_controle);
        return 1;
    }

    clock_t inicio_fase2 = clock();


    for (int i = 0; i < num_threads; i++) {
        workers_contagem[i].id_thread = i;
        workers_contagem[i].contexto = &contexto_contagem;
        pthread_create(&threads[i], NULL, executar_tarefa_contagem_impressao, &workers_contagem[i]);
    }


    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    clock_t fim_total = clock();

    double tempo_calculo = obter_tempo_decorrido(inicio_total, fim_fase1);
    double tempo_contagem = obter_tempo_decorrido(inicio_fase2, fim_total);
    double tempo_total = obter_tempo_decorrido(inicio_total, fim_total);

    printf("Resultados\n\n");
    printf("Total de Primos:        %d\n", contexto_contagem.total_primos_encontrados);
    printf("Tempo Marcacao (F1):    %.4f seg\n", tempo_calculo);
    printf("Tempo Contagem (F2):    %.4f seg\n", tempo_contagem);
    printf("TEMPO TOTAL:            %.4f seg\n", tempo_total);


    pthread_mutex_destroy(&contexto_contagem.mutex_controle);
    free(threads);
    free(workers_contagem);
    free(mapa_primos);

    return 0;
}


// Implementações
int ler_inteiro_positivo(const char* mensagem) {
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

double obter_tempo_decorrido(clock_t inicio, clock_t fim) {
    return ((double)(fim - inicio)) / CLOCKS_PER_SEC;
}

//fase 1
void* executar_tarefa_marcacao(void* arg) {
    WorkerMarcacao* worker = (WorkerMarcacao*)arg;
    ContextoMarcacao* ctx = worker->contexto;

    while (true) {
        int primo_semente = -1;


        pthread_mutex_lock(&ctx->mutex_controle);

        int candidato = ctx->proximo_primo_semente;

        while (candidato <= ctx->limite_raiz_quadrada) {
            if (ctx->mapa_de_primos[candidato]) {
                primo_semente = candidato;
                ctx->proximo_primo_semente = candidato + 1;
                break;
            }
            candidato++;
        }


        if (primo_semente == -1) {
            ctx->proximo_primo_semente = candidato;
        }

        pthread_mutex_unlock(&ctx->mutex_controle);


        if (primo_semente == -1) {
            break;
        }


        long long inicio_risco = (long long)primo_semente * primo_semente;

        if (inicio_risco <= ctx->limite_superior_global) {

            for (int j = (int)inicio_risco; j <= ctx->limite_superior_global; j += primo_semente) {
                ctx->mapa_de_primos[j] = false;
            }
        }
    }
    return NULL;
}

void orquestrar_fase_marcacao(bool* mapa_primos, int limite_superior, int num_threads) {

    for (int i = 0; i <= limite_superior; i++) mapa_primos[i] = true;
    mapa_primos[0] = false;
    mapa_primos[1] = false;

    int limite_raiz = (int)sqrt(limite_superior);

    ContextoMarcacao contexto = {
        .proximo_primo_semente = 2,
        .limite_raiz_quadrada = limite_raiz,
        .limite_superior_global = limite_superior,
        .mapa_de_primos = mapa_primos
    };
    pthread_mutex_init(&contexto.mutex_controle, NULL);

    pthread_t* threads = (pthread_t*)malloc(num_threads * sizeof(pthread_t));
    WorkerMarcacao* workers = (WorkerMarcacao*)malloc(num_threads * sizeof(WorkerMarcacao));

    
    if (threads == NULL || workers == NULL) {
        printf("Erro Fatal: Memoria insuficiente para threads da fase 1.\n");
        free(threads);
        free(workers);
        pthread_mutex_destroy(&contexto.mutex_controle);
        return;
    }

    
    for (int i = 0; i < num_threads; i++) {
        workers[i].id_thread = i;
        workers[i].contexto = &contexto;
        pthread_create(&threads[i], NULL, executar_tarefa_marcacao, &workers[i]);
    }

    
    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    free(threads);
    free(workers);
    pthread_mutex_destroy(&contexto.mutex_controle);
}

//fase 2
Intervalo reservar_proximo_intervalo(ContextoContagem* contexto) {
    Intervalo intervalo = { .inicio = -1, .fim = -1 };

    
    pthread_mutex_lock(&contexto->mutex_controle);

    if (contexto->proximo_indice_leitura <= contexto->limite_superior_global) {
        intervalo.inicio = contexto->proximo_indice_leitura;
        intervalo.fim = intervalo.inicio + TAMANHO_BLOCO_LEITURA - 1;

        // Ajuste se o bloco passar do final
        if (intervalo.fim > contexto->limite_superior_global) {
            intervalo.fim = contexto->limite_superior_global;
        }

        contexto->proximo_indice_leitura = intervalo.fim + 1;
    }

    // CORREÇÃO AQUI: Usar -> em vez de .
    pthread_mutex_unlock(&contexto->mutex_controle);
    return intervalo;
}

void* executar_tarefa_contagem_impressao(void* arg) {
    WorkerContagem* worker = (WorkerContagem*)arg;
    ContextoContagem* ctx = worker->contexto;

    int contador_local = 0;

    
    
    // otimização pro printf acumula texto pra só depois mandar pra tela diminuindo reduz a disputa pelo "cadeado" do terminal.
    char buffer[65536];
    int pos = 0;

    while (true) {
        Intervalo bloco = reservar_proximo_intervalo(ctx);

        
        if (bloco.inicio == -1) break;

        
        for (int numero = bloco.inicio; numero <= bloco.fim; numero++) {

            // Se o buffer estiver quase cheio (margem de 100 bytes), despeja na tela
            if (pos > 65400) {
                fwrite(buffer, 1, pos, stdout);
                pos = 0;
            }

            if (ctx->mapa_de_primos[numero]) {
                contador_local++;
                //sprintf escreve na memória RAM (muito rápido e paralelo)
                pos += sprintf(buffer + pos, "Thread %d: O numero %d e primo\n", worker->id_thread, numero);
            }
            else {
                pos += sprintf(buffer + pos, "Thread %d: O numero %d nao e primo\n", worker->id_thread, numero);
            }
        }
    }

    
    if (pos > 0) {
        fwrite(buffer, 1, pos, stdout);
    }

    
    if (contador_local > 0) {
        pthread_mutex_lock(&ctx->mutex_controle);
        ctx->total_primos_encontrados += contador_local;
        pthread_mutex_unlock(&ctx->mutex_controle);
    }
    return NULL;
}


/*Explicação código
* 
* Logo que o main começa ele pede o intervalo de números, e a quantidade de threads a ser criada 
* proximo passo do main eu aloco memória pro mapa primos (parte mais pesada do código e o porque dele ser ruim pq usa muita memoria) esse mapa vai ter tamanho n+1 pra facilitar a leitura do código por exemplo mapa_primos[100] é o próprio numero 100
* quando mapa primos for true significa que o numero é primo
* 
* aqui começa a fase 1
* o main chama então osquestrar marcação, ele seta todo o mapa como true, logo em seguida seta 0 e 1 que ja sabemos não serem primos
* inicializa struct ContextoMarcacao que guarda o proximo_primo_semente a ser usado (o primeiro é o 2)
* 
* agora que vem a criação das threads com executar_tarefa_marcacao
* como de praxe uso mutex pra uma thread pgar um primo semente por vez
* exemplo: thread 1 pega o numero 2, e sai riscando todos seus multiplos pq esses não são primos (seta os multplos pra false)
* enquanto isso thread 2 pega o 3 e sai riscando os seus multiplos, e assim até acbaar o trabalho aqui
* nada é impresso nessa parte
* 
* dai começa a fase 2
* mapa_primos ja foi escrito, agora só vai ser lido
* o main inicializa ContextoContagem começando do numero 1
* threads são criadas novamente, agora elas pegam lotes de leitura
* thread 1 pega de 1 a 1000, a 2 de 1001-2000, e assim por diante
* cada thread percorre seu lote e olha se mapa_primos[i] é true -> imprime é primo e conta +1 | se e false imprime não é primo
* quando uma thread termina ela soma quantos primos acho no contador global
* 
*/