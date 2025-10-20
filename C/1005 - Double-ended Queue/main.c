#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#define MAX_SIZE 10

typedef struct{
    int front;
    int back;
    int data[MAX_SIZE];
} Deque;

Deque* createDeque(){
    Deque* d = malloc(sizeof *d);
    if(d==NULL){
        fprintf(stderr, "Erro: malloc falhou\n");
        exit(EXIT_FAILURE);
    }
    d->front = -1;
    d->back = -1;
    return d;
}

bool emptyDeque(Deque* d){
    if(d == NULL) return true;
    return d-> front ==-1;
}

bool fullDeque(Deque* d){
    if(d==NULL) return false;
    return ((d->back + 1) % MAX_SIZE) == d->front;
}

bool pushBack(Deque* d, int value){
    if (d == NULL) return false;
    if(fullDeque(d)){
        return false;
    }

    if(emptyDeque(d)){
        d->front=0;
        d->back=0;
        d->data[0]=value;
    }else{
        d->back=(d->back + 1)% MAX_SIZE;
        d->data[d->back] = value;
    }
    return true;
}

bool pushFront(Deque* d, int value){
    if(d==NULL) return false;
    if(fullDeque(d)){
        return false;
    }
    if(emptyDeque(d)){
        d->front=0;
        d->back=0;
        d->data[0]=value;
    }else{
        d->front = (d->front - 1 + MAX_SIZE)% MAX_SIZE;
        d->data[d->front] = value;
    }
    return true;
}


int popFront(Deque* d){
    if(d==NULL || emptyDeque(d)){
        fprintf(stderr,"Erro: deque vazio\n");
        exit(EXIT_FAILURE);
    }

    int val = d->data[d->front];
    if(d->front == d-> back){
        d->front = -1;
        d->back = -1;
    }else{
        d->front = (d->front + 1) % MAX_SIZE;
    }
    return val;
}

int popBack(Deque* d){
    if(d==NULL || emptyDeque(d)){
        fprintf(stderr, "Erro: deque vazio\n");
        exit(EXIT_FAILURE);
    }

    int val = d->data[d->back];
    if(d->front == d->back){
        d->front = -1;
        d->back = -1;
    }else{
        d->back = (d->back - 1 + MAX_SIZE) % MAX_SIZE;
    }
    return val;
}

int peekFront(Deque* d){
    if(d==NULL || emptyDeque(d)){
        fprintf(stderr, "Erro: deque vazio\n " );
        exit(EXIT_FAILURE);
    }
    return d->data[d->front];
}

int peekBack(Deque* d){
    if (d==NULL || emptyDeque(d)){
        fprintf(stderr, "Erro: Deque vazio\n");
        exit(EXIT_FAILURE);
    }
    return d->data[d->back];
}

void destroyDeque(Deque* d){
    if(d==NULL) return;
    free(d);
}

int main() {
    Deque *d = createDeque();
    int valor;

    while (1) {
        if (scanf("%d", &valor) != 1) break;
        if (valor == -1) break;             
        
        if (!pushBack(d, valor)) {
            printf("Aviso: deque cheio. Ignorando %d\n", valor);
        }
    }

    for (int i = 0; i < 100; i++) {
        if (!emptyDeque(d)) {
     
            int atual = popFront(d);
  
            if (atual > 0) {
       
                if (!pushBack(d, atual)) {
         
                    printf("Aviso: nao foi possivel reenfileirar %d\n", atual);
                }
            }
        } else {
  
            break;
        }
    }

    printf("Conteudo restante:\n");
    while (!emptyDeque(d)) {
        printf("%d ", popFront(d));
    }
    printf("\n");

    destroyDeque(d);
    return 0;
}