#pragma once

#include "common.h"

typedef struct Board {
	PieceType squares[BOARD_SIZE][BOARD_SIZE];
    bool white_king_moved;
    bool black_king_moved;
    bool white_rook_left_moved;
    bool white_rook_right_moved;
    bool black_rook_left_moved;
    bool black_rook_right_moved;
    Position en_passant_target;
}Board;

void board_init(Board* board);
void board_reset(Board* board);

PieceType board_get_piece(const Board* board, Position pos);
void board_set_piece(Board* board, Position pos, PieceType piece);


void board_make_move(Board* board, Move move);
void board_copy(const Board* source, Board* dest);

void board_print(const Board* board);