#include "rules.h"
#include "pieces.h"

Position rules_find_king(const Board* board, Player player) {
    PieceType king = (player == PLAYER_WHITE) ? WHITE_KING : BLACK_KING;

    for (int row = 0; row < BOARD_SIZE; row++) {
        for (int col = 0; col < BOARD_SIZE; col++) {
            Position pos = { row, col };
            if (board_get_piece(board, pos) == king) {
                return pos;
            }
        }
    }

    
    Position invalid = { -1, -1 };
    return invalid;
}

bool rules_is_square_attacked(const Board* board, Position pos, Player by_player) {
    for (int row = 0; row < BOARD_SIZE; row++) {
        for (int col = 0; col < BOARD_SIZE; col++) {
            Position attacker_pos = { row, col };
            PieceType attacker = board_get_piece(board, attacker_pos);

            if (attacker == EMPTY) continue;
            if (get_piece_owner(attacker) != by_player) continue;

            if (attacker == WHITE_KING || attacker == BLACK_KING) {
                if (abs(pos.row - row) <= 1 && abs(pos.col - col) <= 1) {
                    return true;
                }
                continue;
            }

            if (piece_can_move_to(board, attacker_pos, pos)) {
                return true;
            }
        }
    }

    return false;
}

bool rules_is_check(const Board* board, Player player) {
    Position king_pos = rules_find_king(board, player);
    if (!is_valid_position(king_pos)) return false;

    Player opponent = (player == PLAYER_WHITE) ? PLAYER_BLACK : PLAYER_WHITE;
    return rules_is_square_attacked(board, king_pos, opponent);
}

bool rules_move_causes_check(const Board* board, Move move, Player player) {

    Board temp_board;
    board_copy(board, &temp_board);
    board_make_move(&temp_board, move);


    return rules_is_check(&temp_board, player);
}

bool rules_is_valid_move(const Board* board, Move move, Player player) {
    PieceType piece = board_get_piece(board, move.from);


    if (piece == EMPTY || get_piece_owner(piece) != player) {
        return false;
    }

 
    if (!piece_can_move_to(board, move.from, move.to)) {
 
        if ((piece == WHITE_KING || piece == BLACK_KING) &&
            move.from.row == move.to.row && abs(move.to.col - move.from.col) == 2) {

            if (move.to.col == 6) { 
                if (!rules_can_castle_kingside(board, player)) return false;
            }
            else if (move.to.col == 2) {
                if (!rules_can_castle_queenside(board, player)) return false;
            }
            else {
                return false;
            }

            move.is_castling = true;
        }
        else {
            return false;
        }
    }

 
    if (rules_move_causes_check(board, move, player)) {
        return false;
    }

    return true;
}

bool rules_can_castle_kingside(const Board* board, Player player) {
    int row = (player == PLAYER_WHITE) ? 7 : 0;

    if (player == PLAYER_WHITE) {
        if (board->white_king_moved || board->white_rook_right_moved) return false;
    }
    else {
        if (board->black_king_moved || board->black_rook_right_moved) return false;
    }

    Position pos1 = { row, 5 };
    Position pos2 = { row, 6 };
    if (board_get_piece(board, pos1) != EMPTY || board_get_piece(board, pos2) != EMPTY) {
        return false;
    }

    if (rules_is_check(board, player)) return false;

    Player opponent = (player == PLAYER_WHITE) ? PLAYER_BLACK : PLAYER_WHITE;
    if (rules_is_square_attacked(board, pos1, opponent) ||
        rules_is_square_attacked(board, pos2, opponent)) {
        return false;
    }

    return true;
}

bool rules_can_castle_queenside(const Board* board, Player player) {
    int row = (player == PLAYER_WHITE) ? 7 : 0;

    if (player == PLAYER_WHITE) {
        if (board->white_king_moved || board->white_rook_left_moved) return false;
    }
    else {
        if (board->black_king_moved || board->black_rook_left_moved) return false;
    }

    Position pos1 = { row, 1 };
    Position pos2 = { row, 2 };
    Position pos3 = { row, 3 };
    if (board_get_piece(board, pos1) != EMPTY ||
        board_get_piece(board, pos2) != EMPTY ||
        board_get_piece(board, pos3) != EMPTY) {
        return false;
    }


    if (rules_is_check(board, player)) return false;


    Player opponent = (player == PLAYER_WHITE) ? PLAYER_BLACK : PLAYER_WHITE;
    if (rules_is_square_attacked(board, pos2, opponent) ||
        rules_is_square_attacked(board, pos3, opponent)) {
        return false;
    }

    return true;
}

bool rules_is_checkmate(const Board* board, Player player) {
  
    if (!rules_is_check(board, player)) {
        return false;
    }

    for (int from_row = 0; from_row < BOARD_SIZE; from_row++) {
        for (int from_col = 0; from_col < BOARD_SIZE; from_col++) {
            Position from = { from_row, from_col };
            PieceType piece = board_get_piece(board, from);

            if (piece == EMPTY || get_piece_owner(piece) != player) continue;

            for (int to_row = 0; to_row < BOARD_SIZE; to_row++) {
                for (int to_col = 0; to_col < BOARD_SIZE; to_col++) {
                    Position to = { to_row, to_col };
                    Move move = { from, to, EMPTY, false, false, false };

                    if (rules_is_valid_move(board, move, player)) {
                        return false;
                    }
                }
            }
        }
    }

    return true;
}

bool rules_is_stalemate(const Board* board, Player player) {
    /* Verifica se NÃO está em xeque */
    if (rules_is_check(board, player)) {
        return false;
    }

    /* Verifica se não tem movimentos válidos */
    for (int from_row = 0; from_row < BOARD_SIZE; from_row++) {
        for (int from_col = 0; from_col < BOARD_SIZE; from_col++) {
            Position from = { from_row, from_col };
            PieceType piece = board_get_piece(board, from);

            if (piece == EMPTY || get_piece_owner(piece) != player) continue;

            for (int to_row = 0; to_row < BOARD_SIZE; to_row++) {
                for (int to_col = 0; to_col < BOARD_SIZE; to_col++) {
                    Position to = { to_row, to_col };
                    Move move = { from, to, EMPTY, false, false, false };

                    if (rules_is_valid_move(board, move, player)) {
                        return false;  /* Encontrou um movimento válido */
                    }
                }
            }
        }
    }

    return true;  /* Afogamento - sem movimentos válidos mas não está em xeque */
}