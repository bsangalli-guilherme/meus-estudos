#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <stdbool.h>


#define MAX_WORD_LENGTH 50
#define MAX_TRIES 6

struct WordWithHint{
    char word[MAX_WORD_LENGTH];
    char hint[MAX_WORD_LENGTH];
};


void displayWord(const char [], const bool guessed[]);

void drawHangman(int tries);

