/*

Given an integer array nums, return true if any value appears at least twice in the array, and return false if every element is distinct.



Example 1:

Input: nums = [1,2,3,1]

Output: true

Explanation:

The element 1 occurs at the indices 0 and 3.

Example 2:

Input: nums = [1,2,3,4]

Output: false

Explanation:

All elements are distinct.

Example 3:

Input: nums = [1,1,1,3,3,4,3,2,4,2]

Output: true



Constraints:

1 <= nums.length <= 10^5
-10^9 <= nums[i] <= 10^9
* 
* 
* 
*/




/* How to solve?
 *
 * 1st try BRUTE FORCE
 * I can try a naive approach first using nested loops to compare every element.
 *
 * for(int i = 0; i < nums.size(); i++) {
 * for (int j = i + 1; j < nums.size(); j++) {
 * if(nums[i] == nums[j]) {
 * return true;
 * }
 * }
 * }
 * return false;
 *
 * PROBLEM: While this works logically, the time complexity is O(N^2).
 * For large datasets, this will result in a "Time Limit Exceeded" (TLE).
 * So, I need something more sophisticated.
 *
 *
 * 2nd try SORTING
 * First, I sort the vector in ascending order.
 * This places duplicate elements next to each other (adjacent).
 * sort(nums.begin(), nums.end());
 *
 * Then, I iterate through the sorted array.
 * I start from index 1 to avoid "out of bounds" when checking the previous element.
 *
 * if(nums[i - 1] == nums[i])
 * return true; // Duplicate found
 *
 * This is fast (O(N log N)) and uses little extra memory (O(1)).
 * However, there is a drawback: it mutates (modifies) the original order of the numbers.
 * Useful when i have less memory
 *
 * 3rd try HASH SET
 * I used an unordered_set.
 * This acts as a memory of "seen" numbers.
 *
 * 1. Create `unordered_set<int> seen`.
 * 2. Use a range-based for loop to iterate through the vector: `for (int num : nums)`
 * 3. Inside the loop:
 * - Check if the number exists: `if (seen.count(num))` -> Return true.
 * - If not, insert the number: `seen.insert(num)`.
 *
 * 4. Return false if the loop finishes without finding a duplicate.
 */

#include <iostream>
#include <vector>
//#include <algorithm> need this one for the sort solution
#include <unordered_set>

using namespace std;

class Solution {
public:
    bool containsDuplicate(vector<int>& nums) {
        unordered_set<int> seen;

        for (int num : nums) {
            if (seen.count(num)) {
                return true;
            }
            seen.insert(num);
        }
        return false;
    }
};



int main() {
    Solution solution;

    vector<int> nums1 = { 3, 4, 3, 1, 6, 8, 9 };

    vector<int> nums2 = { 100, 99, 20, 0, 3 };

    bool resultNums1 = solution.containsDuplicate(nums1);

    bool resultNums2 = solution.containsDuplicate(nums2);

    cout << "Test case 1: " << resultNums1 << endl;
    cout << "Test case 2: " << resultNums2 << endl;

}