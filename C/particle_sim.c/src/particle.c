#include <stddef.h>
#include <stdint.h>


typedef float real_t;

typedef struct {
    real_t x,y;    //posição 
    real_t vx,vy; //velocidade
    real_t mass; //massa != 0
    real_t raio;//raio p colisões
    uint32_t color;
    int alive; // 0 = livre 1 = ativa
    //real_t life //tempo de vida(em segundos)((opcional))
} Particle;

