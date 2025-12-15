#pragma once


#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>


#define BOARD_SIZE 8
#define SQUARE_SIZE 80
#define WINDOW_WIDTH (BOARD_SIZE * SQUARE_SIZE)
#define WINDOW_HEIGHT (BOARD_SIZE * SQUARE_SIZE + 60)

#define COLOR_LIGHT_SQUARE 240, 217, 181
#define COLOR_DARK_SQUARE 181, 136, 99
#define COLOR_SELECTED 255, 255, 0, 128
#define COLOR_VALID_MOVE 0, 255, 0, 128
#define COLOR_CHECK 255, 0, 0, 128
#define COLOR_BACKGROUND 50, 50, 50

typedef enum PieceType {
	EMPTY = 0,
	WHITE_PAWN, WHITE_ROOK, WHITE_KNIGHT, WHITE_BISHOP, WHITE_QUEEN, WHITE_KING,
	BLACK_PAWN, BLACK_ROOK, BLACK_KNIGHT, BLACK_BISHOP, BLACK_QUEEN, BLACK_KING
} PieceType;

typedef enum Player {
	PLAYER_WHITE = 0,
	PLAYER_BLACK = 1,
}Player;

typedef struct Position {
	int row;
	int col;
} Position;

typedef struct Move {
	Position from;
	Position to;
	PieceType captured_piece;
	bool is_castling;
	bool is_en_passant;
	bool is_promotion;
} Move;

typedef enum GameStatus {
	GAME_PLAYING,
	GAME_CHECK,
	GAME_CHECKMATE,
	GAME_STALEMATE,
	GAME_DRAW
}GameStatus;

static inline bool is_white_piece(PieceType piece) {
	return piece >= WHITE_PAWN && piece <= WHITE_KING;
}

static inline bool is_black_piece(PieceType piece) {
	return piece >= BLACK_PAWN && piece <= BLACK_KING;
}

static inline bool is_valid_position(Position pos) {
	return pos.row >= 0 && pos.row < BOARD_SIZE && pos.col >= 0 && pos.col < BOARD_SIZE;
}

static inline Player get_piece_owner(PieceType piece) {
	return is_white_piece(piece) ? PLAYER_WHITE : PLAYER_BLACK;
}
