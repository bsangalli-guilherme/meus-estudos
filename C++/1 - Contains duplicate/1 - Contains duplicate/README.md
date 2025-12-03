# 🔍 Contains Duplicate (LeetCode 217)

## Descrição do problema

Dado um array de inteiros `nums`, retorne `true` se algum valor aparecer pelo menos duas vezes no array, e retorne `false` se todos os elementos forem distintos.

---

## Exemplos

**Exemplo 1**

```
Input: nums = [1,2,3,1]
Output: true
Explicação: O número 1 aparece duas vezes.
```

**Exemplo 2**

```
Input: nums = [1,2,3,4]
Output: false
Explicação: Todos os números são únicos.
```

---

## Solução — Passo a passo

### 1ª tentativa: Força bruta (naïve)

Comparar cada elemento com todos os outros usando laços aninhados.

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

**Desvantagem:** Muito lento para entradas grandes — complexidade O(n²).

---

### 2ª tentativa: Ordenação

Ordenar o vetor e verificar elementos adjacentes — duplicatas aparecem lado a lado.

```cpp
bool containsDuplicate_sort(std::vector<int> nums) {
    std::sort(nums.begin(), nums.end());
    for (int i = 1; i < (int)nums.size(); i++) {
        if (nums[i] == nums[i - 1]) return true;
    }
    return false;
}
```

**Trade-offs:** Modifica a ordem do array original (a menos que você passe por cópia). Complexidade O(n log n).

---

### 3ª tentativa: Hash Set (ideal/ótimo)

Usar um `unordered_set` para guardar os valores vistos.

```cpp
#include <unordered_set>

bool containsDuplicate_hashset(const std::vector<int>& nums) {
    std::unordered_set<int> seen;
    for (int num : nums) {
        if (seen.count(num)) return true;
        seen.insert(num);
    }
    return false;
}
```

**Vantagens:** Tempo O(n) em média, não altera o array original.

---

## Como executar (Visual Studio)

1. Clone o repositório.
2. Abra o arquivo `.sln` no Visual Studio 2022.
3. Configure o projeto como *Startup Project*.
4. Pressione **F5** para rodar.

```cpp
#include <iostream>
#include <vector>

int main() {
    std::vector<int> nums1 = {1,2,3,1};
    std::vector<int> nums2 = {1,2,3,4};
    std::cout << std::boolalpha;
    std::cout << "nums1 tem duplicata? " << containsDuplicate_hashset(nums1) << "\n";
    std::cout << "nums2 tem duplicata? " << containsDuplicate_hashset(nums2) << "\n";
    return 0;
}
```

---
