#include <stdio.h>
#include "ArrayList.h"

int main() {
    int n1, n2, n3, sentinela;

    while(sentinela!=0){
        scanf("%d %d %d", &n1, &n2, &n3);
        
        ArrayList l = CriarLista(0); //cria uma lista vazia
        //adiciona os valores lidos na lista
        for (int i = 0; i < n1; i++) {
            int valor; //valor a ser lido
            scanf("%d", &valor);
            AdicionaValorInicio(&l, valor); //adiciona o valor no início da lista, sempre na posicao 0, empurrando os outros para a direita
        }

        for (int i = 0; i < n2; i++) {
            int valor; //valor a ser lido
            scanf("%d", &valor);
            AdicionaValor(&l, valor); //adiciona o valor no final da lista, sempre na posicao tamanho, que é a ultima posicao + 1
        }

        for (int i = 0; i < n3; i++) {
            int valor; //valor a ser excluido
            scanf("%d", &valor);
            if(Contem(l, valor)){
                RemoveValorEm(&l, PosicaoDe(l, valor)); //remove o valor da lista, se ele existir, na posicao que ele está
            }
        }

        if(EstaVazia(l)){
            printf("Lista vazia\n");
        } else {
            PrintaLista(l); //imprime a lista
        }
        printf("\n Digite 0 para sair ou outro valor para continuar: ");
        scanf("%d", &sentinela); //lê o valor do sentinela
        
        DesalocaLista(l); //desaloca a memoria da lista
    };
    

    return 0;

}