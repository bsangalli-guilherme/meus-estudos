#include "pieces.h"

static bool is_path_clear(const Board* board, Position from, Position to) {
    int row_dir = (to.row > from.row) ? 1 : (to.row < from.row) ? -1 : 0;
    int col_dir = (to.col > from.col) ? 1 : (to.col < from.col) ? -1 : 0;

    Position current = { from.row + row_dir, from.col + col_dir };

    while (current.row != to.row || current.col != to.col) {
        if (board_get_piece(board, current) != EMPTY) {
            return false;
        }
        current.row += row_dir;
        current.col += col_dir;
    }

    return true;
}

bool is_valid_pawn_move(const Board* board, Position from, Position to) {
    PieceType piece = board_get_piece(board, from);
    PieceType target = board_get_piece(board, to);

    int direction = is_white_piece(piece) ? -1 : 1;
    int start_row = is_white_piece(piece) ? 6 : 1;

    /* Movimento para frente */
    if (from.col == to.col) {
        if (to.row == from.row + direction && target == EMPTY) {
            return true;
        }
        if (from.row == start_row && to.row == from.row + 2 * direction) {
            Position middle = { from.row + direction, from.col };
            if (target == EMPTY && board_get_piece(board, middle) == EMPTY) {
                return true;
            }
        }
        return false;
    }

    /* Captura diagonal */
    if (abs(to.col - from.col) == 1 && to.row == from.row + direction) {
        if (target != EMPTY && get_piece_owner(target) != get_piece_owner(piece)) {
            return true;
        }
        if (to.row == board->en_passant_target.row &&
            to.col == board->en_passant_target.col) {
            return true;
        }
    }

    return false;
}

bool is_valid_rook_move(const Board* board, Position from, Position to) {
    if (from.row != to.row && from.col != to.col) {
        return false;
    }
    return is_path_clear(board, from, to);
}

bool is_valid_knight_move(const Board* board, Position from, Position to) {
    int row_diff = abs(to.row - from.row);
    int col_diff = abs(to.col - from.col);
    return (row_diff == 2 && col_diff == 1) || (row_diff == 1 && col_diff == 2);
}

bool is_valid_bishop_move(const Board* board, Position from, Position to) {
    if (abs(to.row - from.row) != abs(to.col - from.col)) {
        return false;
    }
    return is_path_clear(board, from, to);
}

bool is_valid_queen_move(const Board* board, Position from, Position to) {
    return is_valid_rook_move(board, from, to) || is_valid_bishop_move(board, from, to);
}

bool is_valid_king_move(const Board* board, Position from, Position to) {
    return abs(to.row - from.row) <= 1 && abs(to.col - from.col) <= 1;
}

bool piece_can_move_to(const Board* board, Position from, Position to) {
    if (!is_valid_position(from) || !is_valid_position(to)) {
        return false;
    }

    if (from.row == to.row && from.col == to.col) {
        return false;
    }

    PieceType piece = board_get_piece(board, from);
    PieceType target = board_get_piece(board, to);

    if (piece == EMPTY) {
        return false;
    }

    if (target != EMPTY && get_piece_owner(piece) == get_piece_owner(target)) {
        return false;
    }

    switch (piece) {
    case WHITE_PAWN:
    case BLACK_PAWN:
        return is_valid_pawn_move(board, from, to);
    case WHITE_ROOK:
    case BLACK_ROOK:
        return is_valid_rook_move(board, from, to);
    case WHITE_KNIGHT:
    case BLACK_KNIGHT:
        return is_valid_knight_move(board, from, to);
    case WHITE_BISHOP:
    case BLACK_BISHOP:
        return is_valid_bishop_move(board, from, to);
    case WHITE_QUEEN:
    case BLACK_QUEEN:
        return is_valid_queen_move(board, from, to);
    case WHITE_KING:
    case BLACK_KING:
        return is_valid_king_move(board, from, to);
    default:
        return false;
    }
}

int piece_get_possible_moves(const Board* board, Position from, Move* moves, int max_moves) {
    int count = 0;

    for (int row = 0; row < BOARD_SIZE && count < max_moves; row++) {
        for (int col = 0; col < BOARD_SIZE && count < max_moves; col++) {
            Position to = { row, col };
            if (piece_can_move_to(board, from, to)) {
                moves[count].from = from;
                moves[count].to = to;
                moves[count].captured_piece = board_get_piece(board, to);
                moves[count].is_castling = false;
                moves[count].is_en_passant = false;
                moves[count].is_promotion = false;
                count++;
            }
        }
    }

    return count;
}