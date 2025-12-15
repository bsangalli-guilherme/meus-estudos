#pragma once

#include "common.h"
#include "board.h"


bool piece_can_move_to(const Board* board, Position from, Position to);

int piece_get_possible_moves(const Board* board, Position from, Move* moves, int max_moves);

bool is_valid_pawn_move(const Board* board, Position from, Position to);
bool is_valid_rook_move(const Board* board, Position from, Position to);
bool is_valid_knight_move(const Board* board, Position from, Position to);
bool is_valid_bishop_move(const Board* board, Position from, Position to);
bool is_valid_queen_move(const Board* board, Position from, Position to);
bool is_valid_king_move(const Board* board, Position from, Position to);