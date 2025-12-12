#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>

using namespace std;
class Solution {
public:
    vector<vector<string>> groupAnagrams(vector<string>& strs) {
        unordered_map<string, vector<string>> map;

        for (string s : strs) {
            string key(26, 0);
            
            for (char c : s) {
                key[c - 'a']++;
            }
            map[key].push_back(s);
        }

        vector<vector<string>> result;
        
        for (auto& pair : map) {
            result.push_back(pair.second);
        }
        return result;
    }
};