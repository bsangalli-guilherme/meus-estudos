#define SDL_MAIN_HANDLED


#include <SDL2/SDL.h>
#include "game_state.h"
#include "renderer.h"

int main(int argc, char* argv[]) {
    /* Inicializa sistemas */
    GameState game;
    Renderer renderer;

    if (!renderer_init(&renderer)) {
        printf("Falha ao inicializar renderizador!\n");
        return 1;
    }

    game_state_init(&game);

    /* Loop principal */
    bool running = true;
    SDL_Event event;

    while (running) {
        /* Processa eventos */
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                running = false;
            }

            if (event.type == SDL_MOUSEBUTTONDOWN && event.button.button == SDL_BUTTON_LEFT) {
                int mouse_x = event.button.x;
                int mouse_y = event.button.y;

                /* Verifica se está esperando promoção */
                if (game.waiting_for_promotion) {
                    /* Calcula qual peça foi clicada no diálogo */
                    int dialog_x = WINDOW_WIDTH / 2 - 150;
                    int dialog_y = WINDOW_HEIGHT / 2 - 50;

                    if (mouse_y >= dialog_y && mouse_y <= dialog_y + 100 &&
                        mouse_x >= dialog_x && mouse_x <= dialog_x + 300) {

                        int piece_index = (mouse_x - dialog_x - 30) / 60;
                        if (piece_index >= 0 && piece_index < 4) {
                            PieceType promoted_pieces[4] = {
                                game.current_player == PLAYER_WHITE ? WHITE_QUEEN : BLACK_QUEEN,
                                game.current_player == PLAYER_WHITE ? WHITE_ROOK : BLACK_ROOK,
                                game.current_player == PLAYER_WHITE ? WHITE_BISHOP : BLACK_BISHOP,
                                game.current_player == PLAYER_WHITE ? WHITE_KNIGHT : BLACK_KNIGHT
                            };
                            game_state_promote_pawn(&game, promoted_pieces[piece_index]);
                        }
                    }
                }
                else {
                    Position clicked = renderer_screen_to_board(mouse_x, mouse_y);

                    if (is_valid_position(clicked)) {
                        game_state_select_square(&game, clicked);
                    }
                }
            }
        }

        /* Renderiza */
        renderer_clear(&renderer);
        renderer_draw_board(&renderer);
        renderer_draw_selection(&renderer, &game);
        renderer_draw_valid_moves(&renderer, &game);
        renderer_draw_pieces(&renderer, &game);
        renderer_draw_status_bar(&renderer, &game);

        if (game.waiting_for_promotion) {
            renderer_draw_promotion_dialog(&renderer, game.current_player);
        }

        renderer_present(&renderer);

        /* FPS Control */
        SDL_Delay(16);  /* ~60 FPS */
    }

    /* Limpeza */
    renderer_destroy(&renderer);

    return 0;
}