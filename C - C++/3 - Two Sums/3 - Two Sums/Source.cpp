#include <iostream>
#include <vector>
#include <unordered_map>

class Solution {
public:
    std::vector<int> twoSum(std::vector<int>& nums, int target) {
        std::unordered_map<int, int> seen;

        for (int i = 0; i < nums.size(); i++) {
            int complement = target - nums[i];

            if (seen.count(complement)) {
                return { seen[complement], i };
            }
            seen[nums[i]] = i;
        }

        
        return {};
    }
};

int main() {
    Solution sol;

    
    std::vector<int> nums = { 2, 7, 11, 15, 7, 3, 6 };
    int target = 9;

    std::vector<int> result = sol.twoSum(nums, target);

    if (!result.empty()) {
        std::cout << "Indices encontrados: "
            << result[0] << " e " << result[1] << std::endl;
        std::cout << "Valores: "
            << nums[result[0]] << " + " << nums[result[1]]
            << " = " << target << std::endl;
    }
    else {
        std::cout << "Nenhuma resposta encontrada!" << std::endl;
    }

    return 0;
}
