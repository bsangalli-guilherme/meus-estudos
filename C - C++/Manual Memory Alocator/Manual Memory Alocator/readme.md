# Simple Memory Manager (Heap Simulator)

A low-level implementation of a memory allocator in C, simulating how functions like malloc and free work under the hood. This project manages a fixed memory pool using a linked list strategy to track free and allocated blocks.

## Overview

This project serves as an educational tool to understand:

- Pointer Arithmetic: Navigating raw memory bytes manually.
- Memory Fragmentation: Visualizing how "holes" appear in memory over time.
- Allocation Strategies: Implementation of the "First-Fit" algorithm with block splitting.

## Features

- Custom Malloc (memory_allocation): Searches for available blocks using a First-Fit strategy and splits large blocks to minimize wasted space.

- Custom Free (memory_free): Marks blocks as free for future reuse.

- Heap Visualization: A dashboard that prints the current state of the memory map, showing addresses, block sizes, and status (FREE/OCCUPIED).


## How It Works

- The Arena: The program allocates a static global array MEMORY_POOL (10KB) acting as the simulated RAM.

- Block Headers: Every allocation is preceded by a hidden BlockHeader structure containing:

- next: Pointer to the next block.

- payload_size: Size of the usable memory.

- is_free: Flag indicating availability.

- Allocation Logic:

- Traverses the linked list to find a free block large enough.

- If the block is significantly larger than requested, it splits the block into two: one for the user and one remaining free block.

- Returns a pointer to the payload (skipping the header).

## TODOs
This project is currently in the Learning Phase. Below is the list of planned features to reach a production-ready simulation:


- [ ] Block Coalescing: Automatically merging adjacent free blocks to reduce external fragmentation upon free().

- [ ] Bitwise Header Optimization:

	- Use the Least Significant Bit (LSB) of the size field to store the is_free flag.

	- This will reduce the header overhead by removing the int is_free field (saving 4 bytes per block).

- [ ] Memory Alignment: Ensure all returned pointers are aligned to 8 bytes (crucial for CPU performance and complex data types).

- [ ] Fragmentation Metric: Add a function to calculate and display the percentage of fragmented memory.