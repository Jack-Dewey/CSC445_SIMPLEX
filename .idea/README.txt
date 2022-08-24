The simplex solving program is locating in the sole submitted file: main.py
To run the code, call the program through "py main.py <INPUT TEST FILE.txt>"
The program will then solve the problem provided by the test file and print out the solution, if one exists.

The implementation uses bland's rule to avoid cycling. It does so by always selecting the lowest eligible index.
We solve initially infeasible problems by solving an auxiliary problem and then returning a auxiliary-solved, equal,
new feasible problem.
We use floating point arithmetic. The solution rounds to 9 significant decimal points.
