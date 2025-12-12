# Valid Anagram (LeetCode 242)

## Problem Description

Given two strings `s` and `t`, return `true` if `t` is an **anagram** of `s`, and `false` otherwise.

An **anagram** uses all letters of another string exactly once.

---

## Examples

### Example 1

```
Input: s = "anagram", t = "nagaram"
Output: true
```

### Example 2

```
Input: s = "rat", t = "car"
Output: false
```

---

## Constraints

* `1 <= s.length, t.length <= 5 * 10^4`
* Strings contain only lowercase English letters (`a`–`z`)

---

# Solution Walkthrough

## 1st try — Sorting

First, check if both strings have the same length. If not, they cannot be anagrams.

Then, sort both strings and compare them.

```cpp
class Solution {
public:
    bool isAnagram(string s, string t) {
        if (s.size() != t.size()) {
            return false;
        }

        sort(s.begin(), s.end());
        sort(t.begin(), t.end());

        return s == t;
    }
};
```

### Pros

* Very easy to implement

### Drawback

* Sorting takes **O(n log n)** time
* Not the most efficient for large inputs

---

## 2nd try — Frequency Count (Optimal Solution)

Check length equality again. Then use a **frequency array** of size 26 to count how many times each letter appears.

Increment for characters in `s` and decrement for characters in `t`. If the strings are anagrams, all frequencies will return to zero.

```cpp
class Solution {
public:
    bool isAnagram(string s, string t) {
        if (s.size() != t.size()) {
            return false;
        }
        
        int count[26] = {};

        for (int i = 0; i < s.size(); i++) {
            count[s[i] - 'a']++;
            count[t[i] - 'a']--;
        }

        for (int i = 0; i < 26; i++) {
            if (count[i] != 0) {
                return false;
            }
        }

        return true;
    }
};
```

### Benefits

* **O(n)** time complexity (linear)
* **O(1)** extra space (only 26 integers)
* Most efficient method for this problem

---
## Bonus try - frequency count in C

the same logic as the above

```c

bool isAnagram(char* s, char* t){
    int sizeS = strlen(s);
    int sizeT = strlen(t);

    if(sizeS != sizeT) return false;

    int frequencyMap[26] = {};

    for (int i = 0 ; i < sizeS; i++){
        frequencyMap[s[i] - 'a']++
        frequencyMap[t[i] - 'a']--
    }
    
    for(int i = 0; i < 26; i++){
        if(frequencyMap[i] != 0){
            return false;
        }
    }
    return true;
}

```