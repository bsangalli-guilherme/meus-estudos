# 🔍 Contains Duplicate (LeetCode 217)

## Problem Description

Given an integer array `nums`, return `true` if any value appears at least twice, and return `false` if every element is distinct.

---

## Examples

**Example 1**

```
Input: nums = [1,2,3,1]
Output: true
Explanation: The number 1 appears twice.
```

**Example 2**

```
Input: nums = [1,2,3,4]
Output: false
Explanation: All numbers are unique.
```

---

## Solution Walkthrough

### 1st Try: Brute Force (Naive)

Compare each number with every other number using nested loops.

```cpp
bool containsDuplicate_bruteforce(const std::vector<int>& nums) {
    for (int i = 0; i < (int)nums.size(); i++) {
        for (int j = i + 1; j < (int)nums.size(); j++) {
            if (nums[i] == nums[j]) return true;
        }
    }
    return false;
}
```

**Drawback:** Extremely slow for large inputs — O(n²) time.

---

### 2nd Try: Sorting

Sort the array and then check if any adjacent values are equal.

```cpp
bool containsDuplicate_sort(std::vector<int> nums) {
    std::sort(nums.begin(), nums.end());
    for (int i = 1; i < (int)nums.size(); i++) {
        if (nums[i] == nums[i - 1]) return true;
    }
    return false;
}
```

**Trade-offs:** O(n log n) time. Modifies the original array unless passed by copy.

---

### 3rd Try: Hash Set (Optimal)

Use an `unordered_set` to store previously seen values.

```cpp
#include <unordered_set>

bool containsDuplicate_hashset(const std::vector<int>& nums) {
    std::unordered_set<int> seen;
    for (int num : nums) {
        if (seen.count(num)) return true; // duplicate found
        seen.insert(num);
    }
    return false;
}
```

**Benefits:** O(n) average time, does not alter the original array.

---