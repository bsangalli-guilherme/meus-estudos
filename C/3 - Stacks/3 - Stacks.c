#include <stdio.h>
#include <stdlib.h>

#define MAX 5

struct stack{
    int top;
    int itens[MAX];
};

void createStack(struct stack* p){
    p->top = -1; //pilha vazia
}

int fullStack(struct stack* p){
    return p->top == MAX - 1;
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
        printf("Pilha vazia!Não existem itens para serem removidos.\n");
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

int main(){
    struct stack p;
    createStack(&p);

    push(&p, 10);
    push(&p, 20);
    push(&p, 30);
    
    printf("Topo da pilha: %d\n", top(&p));
    pop(&p);
    printf("Topo da pilha: %d\n", top(&p));
    pop(&p);
    pop(&p);
    pop(&p);

    return 0;

}