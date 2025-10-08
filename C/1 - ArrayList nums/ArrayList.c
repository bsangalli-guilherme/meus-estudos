#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "ArrayList.h"

/*o que essa função abaixo faz ?  
Cria uma lista do tipo ArrayList, aloca memória para os valores e inicializa o tamanho da lista.
    Como ela faz isso ?
        define o tamanho da lista com o valor passado como parâmetro, aloca memória para um array(vetor) de inteiros
        com tamanho especificado usando calloc(que inicializa todos os elementos com zero) e retorna a lista criada.
*/
ArrayList CriarLista(int tamanho){
    ArrayList l;
    l.tamanho = tamanho;
    l.valores = calloc(tamanho, sizeof(int));
    return l;
}

/*o que essa função abaixo faz ?
Libera a memória alocada para os valores da lista.*/
void DesalocaLista(ArrayList l){
    free(l.valores);
}

/*o que essa função abaixo faz ?
Escreve um valor em uma posição específica da lista.
Exemplo: se a lista tem tamanho 5 e a posição é -1, a função converte -1 para 4 (5 + -1) e escreve o valor na posição 4.
*/
void EscreveValorEm(ArrayList l, int pos, int valor){ //passa a lista por valor pq n vai modificar o tamanho
    if(pos >= l.tamanho || pos < -l.tamanho) //verifica se a posição é válida 
        return;
    if(pos < 0) //se a posição for negativa, converte para positiva para iterar de trás para frente
        l.valores[l.tamanho + pos] = valor;
    else // se a posição for positiva, escreve diretamente na posição
        l.valores[pos] = valor;
}

/*o que essa função abaixo faz ?
Lê e retorna o valor em uma posição específica da lista.
*/
int LeValorEm(ArrayList l, int pos){
    if(pos < 0)//se a posicao for negativa, converte para positiva
        return l.valores[l.tamanho + pos];//retorna o valor na posicao convertida
    else//se a posicao for positiva, retorna diretamente o valor
        return l.valores[pos];
}

/*o que essa função abaixo faz ?
Imprime todos os valores da lista no console.
*/
void PrintaLista(ArrayList l){
    for(int i = 0; i < l.tamanho; i++){//itera sobre todos os elementos da lista
        printf("%d ", LeValorEm(l, i));//imprime o valor na posicao i
    }
}


void AdicionaValorEm(ArrayList *l, int pos, int valor){
    l -> tamanho ++; //aumenta o tamanho da lista
    l -> valores = realloc(l -> valores, l -> tamanho * sizeof(int)); //realoca a memoria para o novo tamanho
    for (int i = l -> tamanho -1; i > pos; i--){ //desloca os elementos para a direita a partir da posicao pos
        l -> valores[i] = l -> valores[i-1]; //atribui o valor do elemento anterior ao elemento atual
    }
    EscreveValorEm(*l, pos, valor); //escreve o valor na posicao pos
}

/*o que essa função abaixo faz ?
Adiciona um valor no início da lista.
*/
void AdicionaValorInicio(ArrayList *l, int valor) { //passa a lista por referencia pq vai modificar o tamanho
    AdicionaValorEm(l, 0, valor); //adiciona o valor na posicao 0, chamando a funcao AdicionaValorEm
}

/* o que essa função abaixo faz ?
Adiciona um valor no final da lista.
*/
void AdicionaValor(ArrayList *l, int valor){ //passa a lista por referencia pq vai modificar o tamanho
    AdicionaValorEm(l, l -> tamanho, valor); //adiciona o valor na ultima posicao
}

/* o que a função abaixo faz?
Remove um valor em uma posição específica da lista.
*/
void RemoveValorEm(ArrayList *l, int pos){ //passa a lista por referencia pq vai modificar o tamanho
    for (int i = pos; i < l -> tamanho; i++){// desloca os elementos para esquerda a partir da posicao pos
        l -> valores[i] = l -> valores[i+1]; //atribui o valor do elemento seguinte ao elemento atual
    }
    l -> tamanho--; //diminui o tamanho da lista
    l -> valores = realloc(l ->valores, l -> tamanho * sizeof(int)); //realoca a memoria para o novo tamanho
}


/* o que a função abaixo faz?
Remove o valor no início da lista.
*/
void RemoveValorInicio(ArrayList *l){
    RemoveValorEm(l,0); //remove o valor na posicao 0, chamando a funcao RemoveValorEm
}

/* o que a função abaixo faz?
Remove o valor no final da lista.
 */
void RemoveValorFim(ArrayList * l){
    RemoveValorEm(l, l -> tamanho - 1);//remove o valor na ultima posicao, chamando a funcao RemoveValorEm
}

/* o que a função abaixo faz?
Retorna o tamanho da lista.
*/
int Tamanho(ArrayList l){
    return l.tamanho; //retorna o tamanho da lista, atraves do campo tamanho da struct
}

/* o que a função abaixo faz?
Verifica se a lista está vazia.
*/
bool EstaVazia(ArrayList l){
    return Tamanho(l) == 0; //chama a função tamanho e retorna true se o tamanho da lista for 0, false caso contrario
}

int PosicaoDe(ArrayList l, int valor){
    for (int i = 0; i < l.tamanho; i++){//itera sobre todos os elementos da lista
        if (LeValorEm(l, i) == valor) //verifica se o valor na posicao i é igual ao valor procurado chamando a funcao LeValorEm
            return i; //retorna a posicao do valor se encontrado
    }
    return -1; //retorna -1 se o valor nao for encontrado
}
/* o que a função abaixo faz?
Verifica se a lista contém um valor específico.
*/
bool Contem(ArrayList l, int valor){
    return PosicaoDe(l, valor) != -1; //chama a funcao PosicaoDe e retorna true se o valor for encontrado, false caso contrario
}