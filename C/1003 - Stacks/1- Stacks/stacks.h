#include <stdio.h>
#include <stdlib.h>

#define MAX 5

typedef struct{
    int top;
    int itens[MAX];
}stack;

void createStack(stack* p);

int fullStack(stack* p);

int emptyStack(stack* p);

void push(stack* p, int item);

int pop(stack* p);

int top (stack* p);
