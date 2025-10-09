#include <stdio.h>
#include <stdlib.h>

#define MAX 5

typedef struct{
    int top;
    int itens[MAX];
}stack;

void createStack(struct stack* p);

int fullStack(struct stack* p);

int emptyStack(struct stack* p);

void push(struct stack* p, int item);

int pop(struct stack* p);

int top (struct stack* p);
