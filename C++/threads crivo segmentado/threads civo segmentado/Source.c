#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define SEGMENT_SIZE 32768
#define COUNT_BATCH_SIZE 1000

typedef struct {
    int* base_primes;
    int total_base_primes;
} BasePrimes;

typedef struct {
    long long next_segment_start;
    long long max_number;
    BasePrimes* base_primes;
    long long total_prime_count;
    pthread_mutex_t mutex;
} SegmentedData;

typedef struct {
    int id;
    SegmentedData* data;
} ThreadData;

int read_valid_user_input(const char* message);
void simple_sieve_base_primes(int limit, BasePrimes* primes);
void* segmented_sieve_thread(void* arg);
double calculate_elapsed_time(clock_t start, clock_t end);

int main() {

    long long max_number = read_valid_user_input("Enter how many numbers to process: ");
    int thread_count = read_valid_user_input("Enter how many threads to create: ");

    clock_t start_total = clock();

    int root_limit = (int)sqrt(max_number);

    //printf("\nPhase 1: Calculating base primes up to %d...\n", root_limit);
    BasePrimes base_primes;
    simple_sieve_base_primes(root_limit, &base_primes);
    //printf("Found %d base primes\n", base_primes.total_base_primes);

    clock_t end_phase1 = clock();

    /*printf("\nPhase 2: Processing segments in parallel...\n");
    printf("Segment size: %d bytes (optimized for L1 Cache)\n", SEGMENT_SIZE);
    printf("Threads: %d\n\n", thread_count);*/

    SegmentedData segmented_data = {
        .next_segment_start = 0,
        .max_number = max_number,
        .base_primes = &base_primes,
        .total_prime_count = 0
    };
    pthread_mutex_init(&segmented_data.mutex, NULL);

    pthread_t* threads = (pthread_t*)malloc(thread_count * sizeof(pthread_t));
    ThreadData* thread_data = (ThreadData*)malloc(thread_count * sizeof(ThreadData));

    clock_t start_phase2 = clock();

    for (int i = 0; i < thread_count; i++) {
        thread_data[i].id = i;
        thread_data[i].data = &segmented_data;
        pthread_create(&threads[i], NULL, segmented_sieve_thread, &thread_data[i]);
    }

    for (int i = 0; i < thread_count; i++) {
        pthread_join(threads[i], NULL);
    }

    clock_t end_total = clock();

    double time_phase1 = calculate_elapsed_time(start_total, end_phase1);
    double time_phase2 = calculate_elapsed_time(start_phase2, end_total);
    double time_total = calculate_elapsed_time(start_total, end_total);

    //printf("Results\n\n");

    printf("Primes found:           %lld\n\n", segmented_data.total_prime_count);
    //printf("Base primes time:       %.4f sec\n", time_phase1);
    //printf("Segment processing time: %.4f sec\n", time_phase2);
    printf("Total execution time:   %.4f sec\n", time_total);

    pthread_mutex_destroy(&segmented_data.mutex);
    free(base_primes.base_primes);
    free(threads);
    free(thread_data);

    return 0;
}

int read_valid_user_input(const char* message) {
    int value;
    do {
        printf("%s", message);
        if (scanf("%d", &value) != 1) {
            while (getchar() != '\n');
            continue;
        }
    } while (value <= 0);
    return value;
}

void simple_sieve_base_primes(int limit, BasePrimes* primes) {
    bool* is_prime = (bool*)malloc((limit + 1) * sizeof(bool));

    for (int i = 0; i <= limit; i++) {
        is_prime[i] = true;
    }
    is_prime[0] = is_prime[1] = false;

    for (int i = 2; i * i <= limit; i++) {
        if (is_prime[i]) {
            for (int j = i * i; j <= limit; j += i) {
                is_prime[j] = false;
            }
        }
    }

    int count = 0;
    for (int i = 2; i <= limit; i++) {
        if (is_prime[i]) count++;
    }

    primes->base_primes = (int*)malloc(count * sizeof(int));
    primes->total_base_primes = 0;

    for (int i = 2; i <= limit; i++) {
        if (is_prime[i]) {
            primes->base_primes[primes->total_base_primes++] = i;
        }
    }

    free(is_prime);
}

void* segmented_sieve_thread(void* arg) {
    ThreadData* thread_data = (ThreadData*)arg;
    SegmentedData* data = thread_data->data;

    bool* segment = (bool*)malloc(SEGMENT_SIZE * sizeof(bool));
    long long local_count = 0;

    while (true) {
        long long seg_start, seg_end;

        pthread_mutex_lock(&data->mutex);

        seg_start = data->next_segment_start;
        seg_end = seg_start + SEGMENT_SIZE - 1;

        if (seg_end > data->max_number) {
            seg_end = data->max_number;
        }

        if (seg_start > data->max_number) {
            pthread_mutex_unlock(&data->mutex);
            break;
        }

        data->next_segment_start = seg_end + 1;

        pthread_mutex_unlock(&data->mutex);

        long long seg_size = seg_end - seg_start + 1;

        for (long long i = 0; i < seg_size; i++) {
            segment[i] = true;
        }

        for (int p = 0; p < data->base_primes->total_base_primes; p++) {
            long long prime = data->base_primes->base_primes[p];

            long long first_multiple;

            if (seg_start <= prime) {
                first_multiple = prime * prime;
            }
            else {
                first_multiple = ((seg_start + prime - 1) / prime) * prime;
            }

            for (long long j = first_multiple; j <= seg_end; j += prime) {
                if (j >= seg_start) {
                    segment[j - seg_start] = false;
                }
            }
        }

        for (long long i = 0; i < seg_size; i++) {
            long long number = seg_start + i;

            if (number < 2) {
                continue;
            }

            if (number == 2) {
                local_count++;
                continue;
            }

            if (number % 2 == 0) {
                continue;
            }

            if (segment[i]) {
                local_count++;
            }
        }
    }

    pthread_mutex_lock(&data->mutex);
    data->total_prime_count += local_count;
    pthread_mutex_unlock(&data->mutex);

    free(segment);
    return NULL;
}

double calculate_elapsed_time(clock_t start, clock_t end) {
    return ((double)(end - start)) / CLOCKS_PER_SEC;
}
