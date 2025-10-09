#include <stdio.h>
#include <stdlib.h>
#include "stacks.h"

#define MAX 5

//cria uma pilha vazia
void createStack(struct stack* p){
    p->top = -1; // como nenhum elemento foi inserido ainda top não pode apontar para nada valido por isso -1
}

int fullStack(struct stack* p){
    return p->top == MAX - 1; // apesar de serem 5 elementos o array começa a contar de 0, entaço tem que fazer 5-1 para o indice do ultimo ser 4
}

int emptyStack(struct stack* p){
    return p->top == -1; 
}

void push(struct stack* p, int item){
    if(fullStack(p)){
        printf("Pilha cheia! Não é possível adicionar o item.\n");
    }else{
        p->itens[++(p->top)] = item; // array começa com valor -1
    }
}

int pop(struct stack* p){
    if(emptyStack(p)){
        printf("Pilha vazia! Não existem itens para serem removidos.\n");
        return -1;
    }else{
        int item = p->itens[(p->top)--];
        printf("Item %d removido da pilha com sucesso.\n", item);
        return item;
    }
}

int top (struct stack* p){
    if(emptyStack(p)){
        printf("Pilha vazia! Não há item no topo.\n");
        return -1;
    }else{
        return p->itens[p->top];
    }
}

