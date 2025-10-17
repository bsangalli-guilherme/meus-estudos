#include <stdio.h>
#include <stdlib.h>
#include "stacks.h"


int main(){
    stack p;
    int valor=0, i;
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