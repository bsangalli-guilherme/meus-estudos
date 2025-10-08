typedef struct No No;
struct No {
    int valor;
    No * prox;
};

typedef struct {
    No * inicio;
    No * fim;
    int tamanho;
} Lista;

No * CriaNo(int valor);
Lista CriaLista();
No * NoEm(Lista * l, int pos);
void AdicionaInicio(Lista * l, int valor);
void AdicionaFim(Lista * l, int valor);
void AdicionaEm(Lista * l, int pos, int valor);
void RemoveInicio(Lista * l);
void RemoveFim(Lista * l);
void RemoveEm(Lista * l, int pos);
void PrintaLista(Lista l);