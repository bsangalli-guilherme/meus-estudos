#include <SDL2/SDL.h>
#include "renderer.h"
#include <stdbool.h>

bool renderer_init(Renderer* rend) {
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        printf("SDL não pode ser inicializado! SDL_Error: %s\n", SDL_GetError());
        return false;
    }

    rend->window = SDL_CreateWindow(
        "Xadrez Profissional em C",
        SDL_WINDOWPOS_CENTERED,
        SDL_WINDOWPOS_CENTERED,
        WINDOW_WIDTH,
        WINDOW_HEIGHT,
        SDL_WINDOW_SHOWN
    );

    if (!rend->window) {
        printf("Janela não pode ser criada! SDL_Error: %s\n", SDL_GetError());
        return false;
    }

    rend->renderer = SDL_CreateRenderer(rend->window, -1, SDL_RENDERER_ACCELERATED);

    if (!rend->renderer) {
        printf("Renderer não pode ser criado! SDL_Error: %s\n", SDL_GetError());
        return false;
    }

    /* Inicializa texturas (por enquanto NULL, renderiza formas geométricas) */
    for (int i = 0; i < 13; i++) {
        rend->piece_textures[i] = NULL;
    }

    rend->initialized = true;
    return true;
}

void renderer_destroy(Renderer* rend) {
    if (!rend->initialized) return;

    for (int i = 0; i < 13; i++) {
        if (rend->piece_textures[i]) {
            SDL_DestroyTexture(rend->piece_textures[i]);
        }
    }

    if (rend->renderer) {
        SDL_DestroyRenderer(rend->renderer);
    }

    if (rend->window) {
        SDL_DestroyWindow(rend->window);
    }

    SDL_Quit();
    rend->initialized = false;
}

void renderer_clear(Renderer* rend) {
    SDL_SetRenderDrawColor(rend->renderer, COLOR_BACKGROUND, 255);
    SDL_RenderClear(rend->renderer);
}

void renderer_draw_board(Renderer* rend) {
    for (int row = 0; row < BOARD_SIZE; row++) {
        for (int col = 0; col < BOARD_SIZE; col++) {
            bool is_light = (row + col) % 2 == 0;

            SDL_Rect square = {
                col * SQUARE_SIZE,
                row * SQUARE_SIZE,
                SQUARE_SIZE,
                SQUARE_SIZE
            };

            if (is_light) {
                SDL_SetRenderDrawColor(rend->renderer, COLOR_LIGHT_SQUARE, 255);
            }
            else {
                SDL_SetRenderDrawColor(rend->renderer, COLOR_DARK_SQUARE, 255);
            }

            SDL_RenderFillRect(rend->renderer, &square);
        }
    }
}

/* Desenha peça como forma geométrica simples */
static void draw_piece_shape(SDL_Renderer* renderer, PieceType piece, int x, int y, int size) {
    int center_x = x + size / 2;
    int center_y = y + size / 2;
    int radius = size / 3;

    /* Define cor da peça */
    if (is_white_piece(piece)) {
        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
    }
    else {
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    }

    /* Desenha círculo preenchido (aproximado com retângulos) */
    for (int w = 0; w < radius * 2; w++) {
        for (int h = 0; h < radius * 2; h++) {
            int dx = radius - w;
            int dy = radius - h;
            if ((dx * dx + dy * dy) <= (radius * radius)) {
                SDL_RenderDrawPoint(renderer, center_x + dx, center_y + dy);
            }
        }
    }

    /* Desenha símbolo distintivo por tipo */
    SDL_Rect symbol_rect;

    switch (piece) {
    case WHITE_KING:
    case BLACK_KING:
        /* Cruz no topo */
        symbol_rect = (SDL_Rect){ center_x - 2, center_y - radius - 8, 4, 8 };
        SDL_RenderFillRect(renderer, &symbol_rect);
        symbol_rect = (SDL_Rect){ center_x - 6, center_y - radius - 6, 12, 4 };
        SDL_RenderFillRect(renderer, &symbol_rect);
        break;

    case WHITE_QUEEN:
    case BLACK_QUEEN:
        /* Pequena coroa */
        for (int i = -1; i <= 1; i++) {
            symbol_rect = (SDL_Rect){ center_x + i * 6 - 1, center_y - radius - 4, 2, 4 };
            SDL_RenderFillRect(renderer, &symbol_rect);
        }
        break;

    case WHITE_ROOK:
    case BLACK_ROOK:
        /* Forma de torre */
        symbol_rect = (SDL_Rect){ center_x - 8, center_y - radius + 2, 16, 4 };
        SDL_RenderFillRect(renderer, &symbol_rect);
        break;

    case WHITE_KNIGHT:
    case BLACK_KNIGHT:
        /* Forma de L */
        symbol_rect = (SDL_Rect){ center_x - 4, center_y - 6, 4, 12 };
        SDL_RenderFillRect(renderer, &symbol_rect);
        symbol_rect = (SDL_Rect){ center_x - 4, center_y + 4, 8, 4 };
        SDL_RenderFillRect(renderer, &symbol_rect);
        break;

    case WHITE_BISHOP:
    case BLACK_BISHOP:
        /* Ponto no topo */
        for (int w = -2; w <= 2; w++) {
            for (int h = -2; h <= 2; h++) {
                if (w * w + h * h <= 4) {
                    SDL_RenderDrawPoint(renderer, center_x + w, center_y - radius + h);
                }
            }
        }
        break;

    case WHITE_PAWN:
    case BLACK_PAWN:
        /* Círculo menor - já renderizado */
        break;

    default:
        break;
    }

    /* Contorno da peça oposta */
    if (is_white_piece(piece)) {
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    }
    else {
        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
    }

    /* Desenha contorno */
    for (int angle = 0; angle < 360; angle += 5) {
        int px = center_x + (int)(radius * cos(angle * 3.14159 / 180.0));
        int py = center_y + (int)(radius * sin(angle * 3.14159 / 180.0));
        SDL_RenderDrawPoint(renderer, px, py);
    }
}

void renderer_draw_pieces(Renderer* rend, const GameState* state) {
    for (int row = 0; row < BOARD_SIZE; row++) {
        for (int col = 0; col < BOARD_SIZE; col++) {
            Position pos = { row, col };
            PieceType piece = board_get_piece(&state->board, pos);

            if (piece != EMPTY) {
                int x = col * SQUARE_SIZE;
                int y = row * SQUARE_SIZE;
                draw_piece_shape(rend->renderer, piece, x, y, SQUARE_SIZE);
            }
        }
    }
}

void renderer_draw_selection(Renderer* rend, const GameState* state) {
    if (!state->has_selection) return;

    SDL_Rect rect = {
        state->selected_square.col * SQUARE_SIZE,
        state->selected_square.row * SQUARE_SIZE,
        SQUARE_SIZE,
        SQUARE_SIZE
    };

    SDL_SetRenderDrawColor(rend->renderer, COLOR_SELECTED);
    SDL_RenderFillRect(rend->renderer, &rect);
}

void renderer_draw_valid_moves(Renderer* rend, const GameState* state) {
    SDL_SetRenderDrawColor(rend->renderer, COLOR_VALID_MOVE);

    for (int i = 0; i < state->valid_moves_count; i++) {
        Position to = state->valid_moves[i].to;

        SDL_Rect rect = {
            to.col * SQUARE_SIZE + SQUARE_SIZE / 4,
            to.row * SQUARE_SIZE + SQUARE_SIZE / 4,
            SQUARE_SIZE / 2,
            SQUARE_SIZE / 2
        };

        SDL_RenderFillRect(rend->renderer, &rect);
    }
}

void renderer_draw_status_bar(Renderer* rend, const GameState* state) {
    SDL_Rect bar = { 0, BOARD_SIZE * SQUARE_SIZE, WINDOW_WIDTH, 60 };
    SDL_SetRenderDrawColor(rend->renderer, 40, 40, 40, 255);
    SDL_RenderFillRect(rend->renderer, &bar);

    /* Aqui você pode adicionar SDL_ttf para texto */
    /* Por simplicidade, apenas desenhamos indicadores coloridos */

    /* Indicador de turno */
    SDL_Rect turn_indicator = { 20, BOARD_SIZE * SQUARE_SIZE + 20, 20, 20 };
    if (state->current_player == PLAYER_WHITE) {
        SDL_SetRenderDrawColor(rend->renderer, 255, 255, 255, 255);
    }
    else {
        SDL_SetRenderDrawColor(rend->renderer, 0, 0, 0, 255);
    }
    SDL_RenderFillRect(rend->renderer, &turn_indicator);

    /* Indicador de xeque */
    if (state->status == GAME_CHECK) {
        SDL_Rect check_indicator = { 60, BOARD_SIZE * SQUARE_SIZE + 20, 40, 20 };
        SDL_SetRenderDrawColor(rend->renderer, 255, 0, 0, 255);
        SDL_RenderFillRect(rend->renderer, &check_indicator);
    }
}

void renderer_draw_promotion_dialog(Renderer* rend, Player player) {
    /* Desenha caixa de diálogo simples */
    SDL_Rect dialog = {
        WINDOW_WIDTH / 2 - 150,
        WINDOW_HEIGHT / 2 - 50,
        300,
        100
    };

    SDL_SetRenderDrawColor(rend->renderer, 200, 200, 200, 255);
    SDL_RenderFillRect(rend->renderer, &dialog);

    /* Desenha opções de peças (Queen, Rook, Bishop, Knight) */
    PieceType pieces[4] = {
        player == PLAYER_WHITE ? WHITE_QUEEN : BLACK_QUEEN,
        player == PLAYER_WHITE ? WHITE_ROOK : BLACK_ROOK,
        player == PLAYER_WHITE ? WHITE_BISHOP : BLACK_BISHOP,
        player == PLAYER_WHITE ? WHITE_KNIGHT : BLACK_KNIGHT
    };

    for (int i = 0; i < 4; i++) {
        int x = dialog.x + 30 + i * 60;
        int y = dialog.y + 25;
        draw_piece_shape(rend->renderer, pieces[i], x, y, 50);
    }
}

void renderer_present(Renderer* rend) {
    SDL_RenderPresent(rend->renderer);
}

Position renderer_screen_to_board(int mouse_x, int mouse_y) {
    Position pos;
    pos.col = mouse_x / SQUARE_SIZE;
    pos.row = mouse_y / SQUARE_SIZE;

    if (!is_valid_position(pos) || pos.row >= BOARD_SIZE) {
        pos.row = -1;
        pos.col = -1;
    }

    return pos;
}