#pragma once

#include "common.h"
#include "board.h"

Position rules_find_king(const Board* board, Player player);


bool rules_is_square_attacked(const Board* board, Position pos, Player by_player);


bool rules_is_check(const Board* board, Player player);

bool rules_move_causes_check(const Board* board, Move move, Player player);


bool rules_is_valid_move(const Board* board, Move move, Player player);


bool rules_can_castle_kingside(const Board* board, Player player);
bool rules_can_castle_queenside(const Board* board, Player player);

bool rules_is_checkmate(const Board* board, Player player);

bool rules_is_stalemate(const Board* board, Player player);