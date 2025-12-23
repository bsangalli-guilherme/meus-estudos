#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stddef.h> 
#include <string.h>

typedef struct BlockHeader {
    struct BlockHeader* next;
    struct BlockHeader* previous;
    size_t payload_size;
    int is_free;
} BlockHeader;

#define POOL_SIZE 10000
char memory_pool[POOL_SIZE];

BlockHeader* heap_head = (BlockHeader*)&memory_pool[0];


void log_heap_operation(const char* operation_type, void* ptr, size_t size);

void initialize_memory_pool() {
    heap_head->next = NULL;
    heap_head->previous = NULL;
    heap_head->payload_size = POOL_SIZE - sizeof(BlockHeader);
    heap_head->is_free = 1;
    FILE* log_file = fopen("log_heap.csv", "w");
}

void* allocate_memory(size_t requested_size) {
    BlockHeader* current_block = heap_head;

    while (current_block != NULL) {
        if (current_block->is_free && current_block->payload_size >= requested_size) {

            if (current_block->payload_size > requested_size + sizeof(BlockHeader)) {
                BlockHeader* new_block = (BlockHeader*)((char*)current_block + sizeof(BlockHeader) + requested_size);

                new_block->payload_size = current_block->payload_size - sizeof(BlockHeader) - requested_size;
                new_block->is_free = 1;

                new_block->next = current_block->next;
                new_block->previous = current_block;

                current_block->payload_size = requested_size;
                current_block->next = new_block;
                if (new_block->next != NULL) {
                    new_block->next->previous = new_block;
                }
            }

            current_block->is_free = 0;
            log_heap_operation("ALLOC", current_block + 1, requested_size);

            return (void*)(current_block + 1);
        }

        current_block = current_block->next;
    }

    return NULL;
}
void coalesce_memory(BlockHeader* block) {
    if (block->next != NULL && block->next->is_free == 1) {
        block->payload_size = block->payload_size + block->next->payload_size + sizeof(BlockHeader);
        block->next = block->next->next;
        if (block->next != NULL) block->next->previous = block;
    }
    if (block->previous != NULL && block->previous->is_free == 1) {
        block->previous->payload_size = block->payload_size + block->previous->payload_size + sizeof(BlockHeader);
        block->previous->next = block->next;
        if (block->next != NULL) block->next->previous = block->previous;
        block = block->previous;
    }

}

void free_memory(void* ptr) {
    if (ptr == NULL) {
        return;
    }

    BlockHeader* header = (BlockHeader*)ptr - 1;
    log_heap_operation("FREE", ptr, header->payload_size);

    header->is_free = 1;

    coalesce_memory(header);

    
}



void visualize_heap() {
    BlockHeader* current_block = heap_head;

    printf("Header Overhead: %zu bytes\n", sizeof(BlockHeader));
    printf("| %-14s | %-10s | %-10s | %s\n\n", "Address", "Size", "Status", "Next");

    size_t total_free_memory = 0;
    size_t total_used_memory = 0;

    while (current_block != NULL) {
        printf("| %p | %8zu B | %-10s | %p |\n",
            (void*)current_block,
            current_block->payload_size,
            current_block->is_free ? "FREE" : "USED",
            (void*)current_block->next);

        if (current_block->is_free) {
            total_free_memory += current_block->payload_size;
        }
        else {
            total_used_memory += current_block->payload_size;
        }

        current_block = current_block->next;
    }

    printf("\nTotal Used: %zu B | Total Free: %zu B\n\n", total_used_memory, total_free_memory);
}

void log_heap_operation(const char* operation_type, void* ptr, size_t size) {
    FILE* log_file = fopen("log_heap.csv", "a");
    if (log_file == NULL) {
        return;
    }

    fprintf(log_file, "%s,%p,%zu\n", operation_type, ptr, size);
    fclose(log_file);
}

int main() {
    initialize_memory_pool();

    void* allocated_pointers[10] = { 0 };
    int user_command = 0;
    size_t allocation_size = 0;
    int slot_index = 0;

    while (1) {
        printf("\n1. Allocate (malloc)\n");
        printf("2. Free\n");
        printf("3. Show Heap Map (visualize)\n");
        printf("4. Exit\n");
        printf("Choose: ");

        if (scanf("%d", &user_command) != 1) {
            while (getchar() != '\n'); // Clear input buffer
            printf(">> Invalid input. Please enter a number.\n");
            continue;
        }

        if (user_command == 1) {
            printf("Size (bytes): ");
            if (scanf("%zu", &allocation_size) != 1) {
                while (getchar() != '\n'); // Clear input buffer
                printf(">> Invalid size.\n");
                continue;
            }

            int available_slot = -1;
            for (int i = 0; i < 10; i++) {
                if (allocated_pointers[i] == NULL) {
                    available_slot = i;
                    break;
                }
            }

            if (available_slot != -1) {
                void* allocated_ptr = allocate_memory(allocation_size);
                if (allocated_ptr != NULL) {
                    allocated_pointers[available_slot] = allocated_ptr;
                    printf(">> Success! Allocated at slot [%d] -> Payload Address: %p\n",
                        available_slot, allocated_ptr);
                }
                else {
                    printf(">> Error: Insufficient memory (or too fragmented)!\n");
                }
            }
            else {
                printf(">> Error: Test slots full (free something from allocated_pointers[])\n");
            }

        }
        else if (user_command == 2) {
            printf("Which slot to free [0-9]? ");
            if (scanf("%d", &slot_index) != 1) {
                while (getchar() != '\n'); // Clear input buffer
                printf(">> Invalid input.\n");
                continue;
            }

            if (slot_index >= 0 && slot_index < 10 && allocated_pointers[slot_index] != NULL) {
                free_memory(allocated_pointers[slot_index]);
                allocated_pointers[slot_index] = NULL;
                printf(">> Slot [%d] freed successfully.\n", slot_index);
            }
            else {
                printf(">> Invalid slot or already empty.\n");
            }

        }
        else if (user_command == 3) {
            visualize_heap();

        }
        else if (user_command == 4) {
            printf("Exiting...\n");
            break;
        }
        else {
            printf(">> Invalid command.\n");
        }
    }

    return 0;
}