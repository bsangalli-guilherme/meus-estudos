#include <iostream>
#include <vector>
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