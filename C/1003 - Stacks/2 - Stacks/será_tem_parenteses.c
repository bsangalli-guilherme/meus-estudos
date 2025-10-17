// check_brackets.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>


//stack definition
typedef struct{
    int head; //top of the stack
    char* data; //char array to store the information (the stack)
} arrayStackChar;

//create an empty stack
arrayStackChar* createStackChar(){
    
    arrayStackChar* s = (arrayStackChar*) malloc(sizeof(arrayStackChar));
    s->head = -1;
    return s;
}

//looks at the top value in the stack
int peek(arrayStackChar* s){
    return s->data[s->head];
}

//checks if the stack is empty if head indez is -1 the stack is empty
bool emptyStack(arrayStackChar* s){
    return s->head == -1;
}


void push(arrayStackChar* s, char value){
    s->head++;
    s->data = (char*) realloc(s->data, (s->head + 1)*sizeof(char));
    s->data[s->head] = value;
}

char pop(arrayStackChar* s){
    char value;
    if(!emptyStack(s)){
        value = peek(s);
        s->head--;
        s->data = (char*) realloc(s->data, (s->head+1)*sizeof(char));
        return value;
    }
    exit(EXIT_FAILURE);
}

int main(void) {
    
    char expr[200];

    fgets(expr, sizeof(expr), stdin);

    arrayStackChar* s = createStackChar();

    for (int i = 0; expr[i]!='\0'; i++){
        char c = expr[i];
        if(c == '(' || c == '[' || '{'){
            push(s, c);
        }else if(c == ')' || c == '}' || '}'){
            if(emptyStack){
                printf("0\n");
                return 0;
            }

            char head = peek(s);
            pop(s);
            
            if((c==')' && head != '(') 
            || (c==']' && head != '[') 
            || (c=='}' && head != '{')){
                printf("0\n");
                return 0;
            }
        }
    }

    printf("%d\n", emptyStack);
    return 0;
}
