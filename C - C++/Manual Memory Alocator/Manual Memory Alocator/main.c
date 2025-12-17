#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stddef.h> 
#include <string.h>




typedef struct BlockHeader {
    struct BlockHeader* next;
    int payload_size;
    int is_free;
} BlockHeader;

#define POOL_SIZE 10000
char MEMORY_POOL[POOL_SIZE];

BlockHeader* global_memory_head = (BlockHeader*)&MEMORY_POOL[0];


void init_memory_pool() {
    global_memory_head->next = NULL;
    global_memory_head->payload_size = POOL_SIZE - sizeof(BlockHeader);
    global_memory_head->is_free = 1;
}


void* memory_allocation(size_t requested_size) {
    BlockHeader* current = global_memory_head;

    while (current != NULL) {

        if (current->is_free == 1 && current->payload_size >= requested_size) {

            
            if (current->payload_size > requested_size + sizeof(BlockHeader)) {

                BlockHeader* next_block = (BlockHeader*)((char*)current + sizeof(BlockHeader) + requested_size);

                
                next_block->payload_size = current->payload_size - sizeof(BlockHeader) - requested_size;
                next_block->is_free = 1;
                next_block->next = current->next;

                
                current->payload_size = requested_size;
                current->next = next_block;
            }

            current->is_free = 0;

            return (void*)(current + 1);
        }

        current = current->next;
    }

    return NULL;
}

void memory_free(void* ptr) {
    if (ptr == NULL) return;

    BlockHeader* header = (BlockHeader*)ptr - 1;

    header->is_free = 1;
}

void visualize_heap() {
    BlockHeader* current = global_memory_head;
    printf("Overhead do Header: %zu bytes", sizeof(BlockHeader));
    printf("| %-14s | %-10s | %-10s | %s\n\n", "Endereco", "Tamanho", "Estado", "Prox");
    int total_free = 0;
    int total_used = 0;

    while (current != NULL) {
        printf("| %p | %8d B | %-10s | %p |\n",
            (void*)current,
            current->payload_size,
            current->is_free ? "LIVRE" : "OCUPADO",
            (void*)current->next);
        if (current->is_free) total_free += current->payload_size;
        else total_used += current->payload_size;

        current = current->next;
    }
    printf("\nTotal Usado: %d B | Total Livre: %d B\n\n", total_used, total_free);
}

int main() {
    init_memory_pool();

    void* pointers[10] = { 0 };
    int command = 0;
    int size = 0;
    int index = 0;
    

    while (1) {
        printf("\n1. Alocar (malloc)\n");
        printf("2. Liberar (free)\n");
        printf("3. Mostrar Mapa (visualize)\n");
        printf("4. Sair\n");
        printf("Escolha: ");

        if (scanf("%d", &command) != 1) break; // Sai se não for número

        if (command == 1) {
            printf("Tamanho (bytes): ");
            scanf("%d", &size);

            // Procura um slot vazio no nosso array de testes
            int found_slot = -1;
            for (int i = 0; i < 10; i++) {
                if (pointers[i] == NULL) { found_slot = i; break; }
            }

            if (found_slot != -1) {
                void* p = memory_allocation(size);
                if (p != NULL) {
                    pointers[found_slot] = p;
                    printf(">> Sucesso! Alocado no slot [%d] -> Endereco Payload: %p\n", found_slot, p);
                }
                else {
                    printf(">> Erro: Sem memoria suficiente (ou fragmentada demais)!\n");
                }
            }
            else {
                printf(">> Erro: Slots de teste cheios (libere algo do array pointers[])\n");
            }

        }
        else if (command == 2) {
            printf("Qual slot liberar [0-9]? ");
            scanf("%d", &index);

            if (index >= 0 && index < 10 && pointers[index] != NULL) {
                memory_free(pointers[index]);
                pointers[index] = NULL;
            }
            else {
                printf(">> Slot invalido ou ja vazio.\n");
            }

        }
        else if (command == 3) {
            visualize_heap();

        }
        else if (command == 4) {
            break;
        }
    }

    return 0;
}