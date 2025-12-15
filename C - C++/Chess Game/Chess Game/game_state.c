#include "game_state.h"
#include "rules.h"

void game_state_init(GameState* state) {
    board_init(&state->board);
    state->current_player = PLAYER_WHITE;
    state->status = GAME_PLAYING;
    state->has_selection = false;
    state->selected_square.row = -1;
    state->selected_square.col = -1;
    state->valid_moves_count = 0;
    state->move_count = 0;
    state->waiting_for_promotion = false;
}

void game_state_reset(GameState* state) {
    game_state_init(state);
}

void game_state_select_square(GameState* state, Position pos) {
    if (!is_valid_position(pos)) return;

    PieceType piece = board_get_piece(&state->board, pos);

    /* Se já tem seleção, tenta mover */
    if (state->has_selection) {
        if (game_state_try_move(state, pos)) {
            return;
        }
    }

    /* Seleciona nova peça */
    if (piece != EMPTY && get_piece_owner(piece) == state->current_player) {
        state->selected_square = pos;
        state->has_selection = true;

        /* Calcula movimentos válidos */
        state->valid_moves_count = 0;

        for (int row = 0; row < BOARD_SIZE; row++) {
            for (int col = 0; col < BOARD_SIZE; col++) {
                Position to = { row, col };
                Move move = { pos, to, EMPTY, false, false, false };

                /* Verifica roque */
                if ((piece == WHITE_KING || piece == BLACK_KING) &&
                    pos.row == to.row && abs(to.col - pos.col) == 2) {
                    move.is_castling = true;
                }

                /* Verifica en passant */
                if ((piece == WHITE_PAWN || piece == BLACK_PAWN) &&
                    to.row == state->board.en_passant_target.row &&
                    to.col == state->board.en_passant_target.col) {
                    move.is_en_passant = true;
                }

                if (rules_is_valid_move(&state->board, move, state->current_player)) {
                    if (state->valid_moves_count < 64) {
                        state->valid_moves[state->valid_moves_count++] = move;
                    }
                }
            }
        }
    }
    else {
        game_state_clear_selection(state);
    }
}

bool game_state_try_move(GameState* state, Position to) {
    if (!state->has_selection) return false;

    /* Procura o movimento nos válidos */
    Move* found_move = NULL;
    for (int i = 0; i < state->valid_moves_count; i++) {
        if (state->valid_moves[i].to.row == to.row &&
            state->valid_moves[i].to.col == to.col) {
            found_move = &state->valid_moves[i];
            break;
        }
    }

    if (!found_move) {
        game_state_clear_selection(state);
        return false;
    }

    /* Salva peça capturada */
    found_move->captured_piece = board_get_piece(&state->board, to);

    /* Verifica promoção de peão */
    PieceType piece = board_get_piece(&state->board, state->selected_square);
    if ((piece == WHITE_PAWN && to.row == 0) ||
        (piece == BLACK_PAWN && to.row == 7)) {
        found_move->is_promotion = true;
        state->waiting_for_promotion = true;
        state->promotion_square = to;
    }

    /* Executa o movimento */
    board_make_move(&state->board, *found_move);

    /* Adiciona ao histórico */
    if (state->move_count < MAX_MOVES_HISTORY) {
        state->move_history[state->move_count++] = *found_move;
    }

    /* Limpa seleção */
    game_state_clear_selection(state);

    /* Troca turno se não está esperando promoção */
    if (!state->waiting_for_promotion) {
        state->current_player = (state->current_player == PLAYER_WHITE) ?
            PLAYER_BLACK : PLAYER_WHITE;
        game_state_update_status(state);
    }

    return true;
}

void game_state_clear_selection(GameState* state) {
    state->has_selection = false;
    state->selected_square.row = -1;
    state->selected_square.col = -1;
    state->valid_moves_count = 0;
}

void game_state_promote_pawn(GameState* state, PieceType new_piece) {
    if (!state->waiting_for_promotion) return;

    board_set_piece(&state->board, state->promotion_square, new_piece);
    state->waiting_for_promotion = false;

    /* Troca turno */
    state->current_player = (state->current_player == PLAYER_WHITE) ?
        PLAYER_BLACK : PLAYER_WHITE;
    game_state_update_status(state);
}

void game_state_update_status(GameState* state) {
    if (rules_is_checkmate(&state->board, state->current_player)) {
        state->status = GAME_CHECKMATE;
    }
    else if (rules_is_stalemate(&state->board, state->current_player)) {
        state->status = GAME_STALEMATE;
    }
    else if (rules_is_check(&state->board, state->current_player)) {
        state->status = GAME_CHECK;
    }
    else {
        state->status = GAME_PLAYING;
    }
}

bool game_state_is_valid_move_target(const GameState* state, Position pos) {
    if (!state->has_selection) return false;

    for (int i = 0; i < state->valid_moves_count; i++) {
        if (state->valid_moves[i].to.row == pos.row &&
            state->valid_moves[i].to.col == pos.col) {
            return true;
        }
    }

    return false;
}
