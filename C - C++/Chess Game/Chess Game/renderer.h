#pragma once

#include <SDL2/SDL.h>
#include "common.h"
#include "game_state.h"

typedef struct {
    SDL_Window* window;
    SDL_Renderer* renderer;
    SDL_Texture* piece_textures[13];  /* Uma para cada tipo de peça + vazio */
    bool initialized;
} Renderer;

/* Inicialização e destruição */
bool renderer_init(Renderer* rend);
void renderer_destroy(Renderer* rend);

/* Renderização */
void renderer_clear(Renderer* rend);
void renderer_draw_board(Renderer* rend);
void renderer_draw_pieces(Renderer* rend, const GameState* state);
void renderer_draw_selection(Renderer* rend, const GameState* state);
void renderer_draw_valid_moves(Renderer* rend, const GameState* state);
void renderer_draw_status_bar(Renderer* rend, const GameState* state);
void renderer_draw_promotion_dialog(Renderer* rend, Player player);
void renderer_present(Renderer* rend);

/* Conversão de coordenadas */
Position renderer_screen_to_board(int mouse_x, int mouse_y);
