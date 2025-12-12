# Two Sum (LeetCode 1)

## Problem Description

Given an array of integers `nums` and an integer `target`, return the **indices of the two numbers** such that they add up to `target`.

You may assume:

* Each input has **exactly one solution**.
* You **cannot** use the same element twice.
* The answer can be returned in any order.

---

## Examples

### Example 1

```
Input: nums = [2,7,11,15], target = 9
Output: [0,1]
Explanation: nums[0] + nums[1] = 9
```

### Example 2

```
Input: nums = [3,2,4], target = 6
Output: [1,2]
```

### Example 3

```
Input: nums = [3,3], target = 6
Output: [0,1]
```

---

## Constraints

* `2 <= nums.length <= 10^4`
* `-10^9 <= nums[i] <= 10^9`
* `-10^9 <= target <= 10^9`
* Only **one valid answer** exists

### Follow-up

Can you design an algorithm better than **O(n²)**?

---

# Solution process

## 1st Try — Brute Force 

Check every pair of numbers using nested loops.

```cpp
class Solution {
public:
    vector<int> twoSum(vector<int>& nums, int target) {
        for (int i = 0; i < nums.size(); i++) {
            for (int j = i + 1; j < nums.size(); j++) {
                if (nums[i] + nums[j] == target) {
                    return {i, j};
                }
            }
        }
        return {}; // Not expected due to problem constraints
    }
};
```

### ❌ Drawback

* **O(n²)** time — very slow for large arrays.

---

## 2nd Try — Hash Map


Use an `unordered_map` to store each number and its index. 
While iterating, compute the complement (`target - nums[i]`) and check if it already exists in the map.

```cpp
#include <unordered_map>

class Solution {
public:
    vector<int> twoSum(vector<int>& nums, int target) {
        unordered_map<int, int> seen; // value -> index

        for (int i = 0; i < nums.size(); i++) {
            int complement = target - nums[i];

            if (seen.count(complement)) {
                return {seen[complement], i};
            }

            seen[nums[i]] = i;
        }

        return {}; // Not expected due to constraints
    }
};
```

### ✔️ Benefits

* **O(n)** time
* **O(n)** space
* Most efficient and widely used solution

---

If you want, I can format this into a full LeetCode solution README or create a version matching your previous templates!
