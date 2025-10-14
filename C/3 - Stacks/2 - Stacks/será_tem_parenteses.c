// check_brackets.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

/* Pilha de chars baseada em array dinâmico */
typedef struct {
    char *data;
    size_t capacity; // capacidade alocada
    size_t top;      // índice do próximo slot livre (== número de elementos)
} CharStack;

bool stack_init(CharStack *s, size_t capacity) {
    s->data = malloc(capacity * sizeof(char));
    if (!s->data && capacity > 0) return false;
    s->capacity = capacity;
    s->top = 0;
    return true;
}

void stack_free(CharStack *s) {
    free(s->data);
    s->data = NULL;
    s->capacity = s->top = 0;
}

bool stack_is_empty(const CharStack *s) {
    return s->top == 0;
}

bool stack_push(CharStack *s, char c) {
    if (s->top >= s->capacity) {
        size_t newcap = (s->capacity == 0) ? 4 : s->capacity * 2;
        char *tmp = realloc(s->data, newcap * sizeof(char));
        if (!tmp) return false;
        s->data = tmp;
        s->capacity = newcap;
    }
    s->data[s->top++] = c;
    return true;
}

bool stack_pop(CharStack *s, char *out) {
    if (stack_is_empty(s)) return false;
    *out = s->data[--s->top];
    return true;
}

char stack_peek(const CharStack *s) {
    if (stack_is_empty(s)) return '\0';
    return s->data[s->top - 1];
}

/* Retorna true se o fechamento 'close' corresponde a abertura 'open' */
bool matches(char open, char close) {
    return (open == '(' && close == ')') ||
           (open == '[' && close == ']') ||
           (open == '{' && close == '}');
}

/* Função principal que verifica balanceamento */
bool is_balanced(const char *s) {
    if (!s) return true; // tratar NULL como "vazio" (balanceado)
    size_t n = strlen(s);

    CharStack st;
    if (!stack_init(&st, n > 0 ? n : 4)) return false; // falha alocação -> considerar inválido

    for (size_t i = 0; i < n; ++i) {
        char ch = s[i];
        if (ch == '(' || ch == '[' || ch == '{') {
            if (!stack_push(&st, ch)) { stack_free(&st); return false; } // falha memória
        } else if (ch == ')' || ch == ']' || ch == '}') {
            char topc;
            if (!stack_pop(&st, &topc)) {
                // fechamento sem abertura correspondente
                stack_free(&st);
                return false;
            }
            if (!matches(topc, ch)) {
                // tipo diferente: exemplo "(]" ou "{)"
                stack_free(&st);
                return false;
            }
        } else {
            // Ignora outros caracteres; se quiser considerar inválidos, trate aqui.
        }
    }

    bool balanced = stack_is_empty(&st);
    stack_free(&st);
    return balanced;
}

/* Programa de teste simples */
int main(void) {
    const char *tests[] = {
        "()", "()[]{}", "(]", "([{}])", "((())", "{[()()]}", "[(])", "", "a(b[c]{d}e)f", NULL
    };

    for (int i = 0; tests[i] != NULL; ++i) {
        printf("Teste: \"%s\" -> %s\n",
               tests[i],
               is_balanced(tests[i]) ? "BALANCEADO" : "INVALIDO");
    }

    /* Exemplo de leitura da entrada (opcional):
    char line[1024];
    printf("Digite uma expressão: ");
    if (fgets(line, sizeof line, stdin)) {
        // remover newline
        line[strcspn(line, "\n")] = '\0';
        printf("%s\n", is_balanced(line) ? "Balanceado" : "Inválido");
    }
    */

    return 0;
}
