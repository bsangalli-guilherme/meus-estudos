#define _CRT_SECURE_NO_WARNINGS 1
#include "board.h"



void board_init(Board* board) {
	memset(board->squares, EMPTY, sizeof(board->squares));

	board->squares[0][0] = BLACK_ROOK;
	board->squares[0][1] = BLACK_KNIGHT;
	board->squares[0][2] = BLACK_BISHOP;
	board->squares[0][3] = BLACK_QUEEN;
	board->squares[0][4] = BLACK_KING;
	board->squares[0][5] = BLACK_BISHOP;
	board->squares[0][6] = BLACK_KNIGHT;
	board->squares[0][7] = BLACK_ROOK;
	
	for (int col = 0; col < BOARD_SIZE; col++) {
		board->squares[1][col] = BLACK_PAWN;
	}
	
	for (int col = 0; col < BOARD_SIZE; col++) {
		board->squares[6][col] = WHITE_PAWN;
	}

	board->squares[7][0] = WHITE_ROOK;
	board->squares[7][1] = WHITE_KNIGHT;
	board->squares[7][2] = WHITE_BISHOP;
	board->squares[7][3] = WHITE_QUEEN;
	board->squares[7][4] = WHITE_KING;
	board->squares[7][5] = WHITE_BISHOP;
	board->squares[7][6] = WHITE_KNIGHT;
	board->squares[7][7] = WHITE_ROOK;


	board->white_king_moved = false;
	board->black_king_moved = false;
	board->white_rook_left_moved = false;
	board->white_rook_right_moved = false;
	board->black_rook_left_moved = false;
	board->black_rook_right_moved = false;
	board->en_passant_target.row = -1;
	board->en_passant_target.col = -1;

}

void board_reset(Board* board) {
	board_init(board);
}

PieceType board_get_piece(const Board* board, Position pos) {
	if (!is_valid_position(pos)) return EMPTY;
	return board->squares[pos.row][pos.col];
}

void board_set_piece(Board* board, Position pos, PieceType piece) {
	if (!is_valid_position(pos)) return;
	board->squares[pos.row][pos.col] = piece;
}

void board_make_move(Board* board, Move move) {
	PieceType piece = board_get_piece(board, move.from);

	
	if (piece == WHITE_KING) board->white_king_moved = true;
	if (piece == BLACK_KING) board->black_king_moved = true;
	if (piece == WHITE_ROOK && move.from.row == 7 && move.from.col == 0)
		board->white_rook_left_moved = true;
	if (piece == WHITE_ROOK && move.from.row == 7 && move.from.col == 7)
		board->white_rook_right_moved = true;
	if (piece == BLACK_ROOK && move.from.row == 0 && move.from.col == 0)
		board->black_rook_left_moved = true;
	if (piece == BLACK_ROOK && move.from.row == 0 && move.from.col == 7)
		board->black_rook_right_moved = true;

	
	if (move.is_castling) {
		board_set_piece(board, move.to, piece);
		board_set_piece(board, move.from, EMPTY);

		
		if (move.to.col == 6) {  
			Position rook_from = { move.from.row, 7 };
			Position rook_to = { move.from.row, 5 };
			PieceType rook = board_get_piece(board, rook_from);
			board_set_piece(board, rook_to, rook);
			board_set_piece(board, rook_from, EMPTY);
		}
		else {
			Position rook_from = { move.from.row, 0 };
			Position rook_to = { move.from.row, 3 };
			PieceType rook = board_get_piece(board, rook_from);
			board_set_piece(board, rook_to, rook);
			board_set_piece(board, rook_from, EMPTY);
		}
		return;
	}

	if (move.is_en_passant) {
		board_set_piece(board, move.to, piece);
		board_set_piece(board, move.from, EMPTY);
		Position captured_pawn = { move.from.row, move.to.col };
		board_set_piece(board, captured_pawn, EMPTY);
		return;
	}

	
	board_set_piece(board, move.to, piece);
	board_set_piece(board, move.from, EMPTY);

	
	board->en_passant_target.row = -1;
	board->en_passant_target.col = -1;

	if ((piece == WHITE_PAWN || piece == BLACK_PAWN) &&
		abs(move.to.row - move.from.row) == 2) {
		board->en_passant_target.row = (move.from.row + move.to.row) / 2;
		board->en_passant_target.col = move.from.col;
	}
}

void board_copy(const Board* source, Board* dest) {
	memcpy(dest, source, sizeof(Board));
}

void board_print(const Board* board) {
	printf("\n a b c d e f g h\n");
	for (int row = 0; row < BOARD_SIZE; row++) {
		printf("%d ", 8 - row);
		for (int col = 0; col < BOARD_SIZE; col++) {
			char symbol = '.';
			PieceType piece = board->squares[row][col];
			switch (piece) {
				case WHITE_PAWN: symbol = 'P'; break;
				case WHITE_ROOK: symbol = 'R'; break;
				case WHITE_KNIGHT: symbol = 'N'; break;
				case WHITE_BISHOP: symbol = 'B'; break;
				case WHITE_QUEEN: symbol = 'Q'; break;
				case WHITE_KING: symbol = 'K'; break;
				case BLACK_PAWN: symbol = 'p'; break;
				case BLACK_ROOK: symbol = 'r'; break;
				case BLACK_KNIGHT: symbol = 'n'; break;
				case BLACK_BISHOP: symbol = 'b'; break;
				case BLACK_QUEEN: symbol = 'q'; break;
				case BLACK_KING: symbol = 'k'; break;
				default: symbol = '.'; break;
			}
			printf("%c ", symbol);
		}
		printf("%d\n", 8 - row);
	}
	printf("  a b c d e f g h\n\n");
}