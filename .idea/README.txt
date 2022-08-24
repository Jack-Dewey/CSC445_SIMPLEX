The simplex solving program is located in main.py
To run the code, call the program through "py main.py <INPUT TEST FILE.txt>"
The first line of the input file will contain values c1 c2 ... cn
Each remaining line will encode a constraint ai1 ai2 ... ain bi
Example:
Max 2x1 + 3x2
st x1 + 2x2 <= 5

Would be:
2 3
1 2 5


The program will then solve the problem provided by the test file and print out the solution, if one exists.

The implementation uses bland's rule to avoid cycling. It does so by always selecting the lowest eligible index.
We solve initially infeasible problems by solving an auxiliary problem and then returning a auxiliary-solved, equal,
new feasible problem.
We use floating point arithmetic. The solution rounds to 9 significant decimal points.
