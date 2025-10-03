#include "media.h"


float media(float vet[], int tamanho) {
	int i = 0;
	float total = 0;
	float media = 0;
	for (; i < tamanho; i++) {
		total += vet[i];
	}
	media = total / tamanho;
	return media;
}