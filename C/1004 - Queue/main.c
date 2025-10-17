#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#define MAX_SIZE 10


//struct queue
typedef struct{
    int front; //indice do inicio da queue 
    int back; //indice final da queue
    int data[MAX_SIZE];
} Queue;

//cria uma queue vazia
Queue* createQueue(){
    //sizeof(Fila) devolve quantos bytes são necessários
    //pra armazenar uma var do tipo queue
    //malloc(...) "compra" um espaço o tamanho de sizeof(...)
    //e devolve um ponteiro void* pro inicio do bloco
    //(Queue*) converte (faz um cast) o tipo void* recebido para o 
    //tipo Queue*
    //Queue* f = ...aqui cria o f que é um ponteiro que aponta para a 
    //queue recem criada
    Queue* q = (Queue*) malloc(sizeof(Queue));
    q -> front = -1; //usamos indice -1 pra não armazenar nada ali
    q -> back = -1;
    return q; //retorna a fila
}

bool fullQueue(Queue* q){
    //calcula a próxima posição onde um novo elemento seria inserido. 
    //Se essa posição for igual a q->front (o índice do elemento no 
    //início da queue), então não há espaço para inserir sem 
    //sobrescrever o front, então a queue está cheia.
    //obs.: se o num for menor que o size, por exemplo 4%5 ele retorna o 4 o que indica que a fila não esta cheia

    return(q->back + 1) % MAX_SIZE == q->front;
}

//bem simples, se o indice do front for -1 é porque está vazia
bool emptyQueue(Queue* q){
    return q->front == -1;
}

void enqueue(Queue* q, int valor){
    if(fullQueue(q)){
        printf("Erro: Fila cheia!");
        return;
    }
    if(emptyQueue(q)){
        q->front = 0;
    }


    q->back = (q->back + 1) % MAX_SIZE;
    q->data[q->back] = valor;

}

int dequeue(Queue* q){
    if(emptyQueue(q)){
        printf("Erro: Fila vazia!");
        exit(EXIT_FAILURE);
    }

    int temp = q->data[q->front];
    if(q->front == q->back){
        q->front = -1;
        q->back = -1;
    }else{
        q->front = (q->front+1) % MAX_SIZE;
    }
    return temp;
}

int main(){
    Queue *q = createQueue();
    int valor;

    while(1){
        if(scanf("%d", &valor)== -1) break;
        if(valor == -1) break;
        enqueue(q, valor);
    }

    for (int i = 0; i < 100; i++){
        if(!emptyQueue(q)){
            int atual = dequeue(q);
            atual --;
            if(atual>0){
                enqueue(q,atual);
            }
        }
    }

    while(!emptyQueue(q)){
        printf("\n%d", dequeue(q));
    }

    return 0;
}