#pragma once

#include "common.h"
#include "board.h"

#define MAX_MOVES_HISTORY 200

typedef struct {
    Board board;
    Player current_player;
    GameStatus status;

    Position selected_square;
    bool has_selection;
    Move valid_moves[64];
    int valid_moves_count;

    
    Move move_history[MAX_MOVES_HISTORY];
    int move_count;

    
    bool waiting_for_promotion;
    Position promotion_square;

} GameState;

void game_state_init(GameState* state);
void game_state_reset(GameState* state);


void game_state_select_square(GameState* state, Position pos);
bool game_state_try_move(GameState* state, Position to);
void game_state_clear_selection(GameState* state);


void game_state_promote_pawn(GameState* state, PieceType new_piece);


void game_state_update_status(GameState* state);


bool game_state_is_valid_move_target(const GameState* state, Position pos);