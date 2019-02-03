# Bachelorarbeit
Project code for my bachelor thesis which solves the binary knapsack problem.

The bachelor thesis is written in German and can be displayed by opening `Bachelorarbeit_Hilbinger.pdf`.

To run the project code, execute the file `KnapsackProblem.exe` in the `Data` folder via a commandline.
The invocation needs to obey the following usage format:

`KnapsackProblem.exe <input.csv> <output.csv> <GA|DP|GR> <parameters.txt>`

An example invocation:

`./KnapsackProblem.exe input/i-7#_170c_1735o.csv out7.csv GA input/p-7#_170c_1735o.txt`

The following algorithms are supported for solving the binary knapsack problem:
- Greedy Algorithm (GR)
- Dynamic Programming (DP)
- Genetic Algorithms (GA)
