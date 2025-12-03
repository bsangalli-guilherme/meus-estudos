#include <iostream>
#include <string>
#include <algorithm>
#include <vector>

using namespace std;

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
		return false;
	}
};


int main() {
	Solution solution;

	string string1 = "anagram";
	string string2 = "maagran";
	string string3 = "maagran";
	string string4 = "esteran";

	cout << boolalpha;

	cout << string1 << " ," << string4 << " : " << solution.isAnagram(string1, string4) << endl;

}