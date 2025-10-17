#include <stdio.h>
#include <stdlib.h>
#include "LinkedList.h"

int main() {
    Lista l = CriaLista();
    RemoveInicio(&l);
    RemoveFim(&l);
    RemoveEm(&l, 10);
    AdicionaInicio(&l, 10);
    RemoveInicio(&l);
    AdicionaInicio(&l, 20);
    AdicionaFim(&l, 30);
    AdicionaFim(&l, 40);
    AdicionaEm(&l, 2, 50);
    AdicionaEm(&l, 0, 1);
    AdicionaEm(&l, 5, 99);
    AdicionaEm(&l, -2, 234234);
    AdicionaEm(&l, 999, 9789678);
    RemoveFim(&l);
    RemoveInicio(&l);
    RemoveEm(&l, 1);
    PrintaLista(l);

    return 0;
}