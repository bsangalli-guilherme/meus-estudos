#include <stdbool.h>

typedef struct ArrayList{
    int tamanho;
    int *valores;
} ArrayList;

ArrayList CriarLista(int tamanho);

void DesalocaLista(ArrayList l);

void EscreveValorEm(ArrayList l, int pos, int valor);

int LeValorEm(ArrayList l, int pos);

void PrintaLista(ArrayList l);

void AdicionaValorEm(ArrayList *l,int pos, int valor);
void AdicionaValorInicio(ArrayList *l, int valor);
void AdicionaValor(ArrayList *l, int valor);

void RemoveValorEm(ArrayList *l, int pos);
void RemoveValorInicio(ArrayList *l);
void RemoveValorFim(ArrayList * l);

int Tamanho(ArrayList l);

bool EstaVazia(ArrayList l);

int PosicaoDe(ArrayList l, int valor);

bool Contem(ArrayList l, int valor);



