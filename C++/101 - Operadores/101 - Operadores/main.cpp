#include <iostream>

main() {

	//arithmetic operations
	int a = 19, b = 7;
	int sum = a + b;
	int difference = a - b;
	int multiplication = a * b;
	int division = a / b;
	float division = float(a) / float(b);
	int remainder = a % b;
	
	
	int c = 0, d, e;
	d = ++c;// increments the value before returning
	e = c++;//returns the value,then increments



	//relational operators
	int f = 0;
	int g = 0;
	g == f; //comparision between 2 values, true if equal

	g != a;///comparision between 2 values, true if notequal

	a > b;
	b < a;
	b <= a;
	a >= b;


	//logical operators

	true && false; //and -> returns true when both are true, false whenever it is something else
	true || true; //or -> returns false only when both are false
	!true && false; //not -> changes true to false, for example !true == not true == false


	//bitwise operations
	//& compares the numbers bit by bit, returns a new where each bit is set to 1 if the bits in both numbers are 1, otherwise the bit is 0
	int result = 5 & 3; // 0000 0101 & 0000 0011 = 0000 0001 == 1

	//| compares the numbers bit by bit, returns a new where each bit is set to 1 if the bits in one of the numbers is 1, otherwise the bit is 0
	int result = 5 | 3; // 0000 0101 | 0000 0011 = 0000 0111 = 7


	//^ compares the numbers bit by bit, returns a new where each bit is set to 1 if the bits in the input numbers are different, otherwise the bit is 0
	int result = 5 ^ 3;// 0000 0101 ^ 0000 0011 = 0000 0110 = 6

	//>> right shifts the first number by the amout of the second  
	int result = 5 >> 1; // 0000 0101 >> 1 = 0000 0010 = 2

	//<< left shifts the first number by the amout of the second 
	int result = 5 << 1; //0000 0101 << 1 = 0000 1010 = 10

	//~ inverts each bit
	int result = ~5; // ~0000 0101 = 1111 1010 = -6
}