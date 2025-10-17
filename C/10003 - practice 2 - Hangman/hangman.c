#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <stdbool.h>
#include "hangman.h"


void displayWord(const char word[], const bool guessed[]){
    printf("Word: ");
    for (int i = 0; word[i]!= '\0'; i++){
        if(guessed[word[i]] - 'a'){
            printf("%c ", word[i]);
        }else{
            printf("_ ");
        }

    }
    printf("\n");
}

void drawHangman(int tries)
{
    const char* hangmanParts[]
        = { "     _________",    "    |         |",
            "    |         O",   "    |        /|\\",
            "    |        / \\", "    |" };

    printf("\n");
    for (int i = 0; i <= tries; i++) {
        printf("%s\n", hangmanParts[i]);
    }
}