#include <stdio.h>
#include <stdlib.h>
#include "LinkedList.h"

No * CriaNo(int valor) {
    No * n = malloc(sizeof(*n));
    if (!n) {
        perror("malloc");
        exit(EXIT_FAILURE);
    }
    n->valor = valor;
    n->prox = NULL;
    return n;
}

Lista CriaLista() {
    Lista l;
    l.inicio = NULL;
    l.fim = NULL;
    l.tamanho = 0;
    return l;
}

No * NoEm(Lista * l, int pos) {
    if (pos < 0 || pos >= l->tamanho) return NULL;
    No * atual = l->inicio;
    for (int i = 0; i < pos; i++) {
        atual = atual->prox;
    }
    return atual;
}

void AdicionaInicio(Lista * l, int valor) {
    No * n = CriaNo(valor);
    if (l->inicio == NULL) {
        l->inicio = l->fim = n;
    } else {
        n->prox = l->inicio;
        l->inicio = n;
    }
    l->tamanho++;
}

void AdicionaFim(Lista * l, int valor) {
    No * n = CriaNo(valor);
    if (l->inicio == NULL) {
        l->inicio = l->fim = n;
    } else {
        l->fim->prox = n;
        l->fim = n;
    }
    l->tamanho++;
}

void AdicionaEm(Lista * l, int pos, int valor) {
    if (pos < 0 || pos > l->tamanho) return;
    if (pos == 0) { AdicionaInicio(l, valor); return; }
    if (pos == l->tamanho) { AdicionaFim(l, valor); return; }
    No * n = CriaNo(valor);
    No * anterior = NoEm(l, pos-1);
    if (anterior == NULL) { /* proteção extra */ 
        free(n);
        return;
    }
    n->prox = anterior->prox;
    anterior->prox = n;
    l->tamanho++;
}

void RemoveInicio(Lista * l) {
    if (l->inicio == NULL) return;
    if (l->tamanho == 1) {
        free(l->inicio);
        l->inicio = l->fim = NULL;
        l->tamanho--;
        return;
    }
    No * remover = l->inicio;
    l->inicio = remover->prox;
    free(remover);
    l->tamanho--;
}

void RemoveFim(Lista * l) {
    if (l->inicio == NULL) return;
    if (l->tamanho == 1) {
        free(l->inicio);
        l->inicio = l->fim = NULL;
        l->tamanho--;
        return;
    }
    No * penultimo = NoEm(l, l->tamanho-2);
    if (penultimo == NULL) return; /* proteção */
    free(l->fim);
    penultimo->prox = NULL;
    l->fim = penultimo;
    l->tamanho--;
}

void RemoveEm(Lista * l, int pos) {
    if (pos < 0 || pos >= l->tamanho) return;
    if (pos == 0) { RemoveInicio(l); return; }
    if (pos == l->tamanho-1) { RemoveFim(l); return; }
    No * anterior = NoEm(l, pos-1);
    if (anterior == NULL || anterior->prox == NULL) return; /* proteção */
    No * remover = anterior->prox;
    anterior->prox = remover->prox;
    free(remover);
    l->tamanho--;
}

void PrintaLista(Lista l) {
    No * atual = l.inicio;
    while (atual != NULL) {
        printf("%d ", atual->valor);
        atual = atual->prox;
    }
    printf("\n");
}
