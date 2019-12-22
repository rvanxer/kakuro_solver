# Kakuro Solver
Kakuro (AKA "cross sum" puzzle) is a Japanese logic puzzle similar to crossword but uses numbers (see [Wikipedia](https://en.wikipedia.org/wiki/Kakuro)). Given a grid of blank cells with barred cells at the left/top of each row/column, the goal of this puzzle is to fill numbers, typically 1 through 9, in these cells without repetition across any row or column such that these numbers add up to the number shown in the barred cells. Kakuro is an NP-complete problem [1], but runs quite quickly for standard, small problems.

![Sample Kakuro problem](https://upload.wikimedia.org/wikipedia/commons/thumb/c/c8/Kakuro_black_box.svg/375px-Kakuro_black_box.svg.png)

This is a single-file Kakuro solver program that uses a host of numerical tricks to shrink the possible "pool" of values each cell can attain. If the problem is difficult enough to be unable to be solved by the pool reduction techniques, the program resorts to a pseudo-brute force approach to solve the puzzle. However, since this is not a complete brute force, correct solution is not always guaranteed (though is quite rare for standard 16x16 puzzles).

## Usage 
The program requires a CSV file as input of the problem in a specific way. For reference, see the `Problems` folder, where a few 10x10 sample puzzles are provided. These problems and their solutions are organized in a more legible way in `Problems/Sample Problems with Solution.xlsx`. To run the solver, instantiate a `Problem` object in `main.py` by providing the path of the CSV file of the puzzle. Alternatively, uncomment the last line of `main.py` and run the file.

Requires Python 3.x.

## References
[1]: Akahiro, Seta (February 5, 2002). "The complexities of puzzles, cross sum and their another solution problems (ASP)" (see [PDF](http://www-imai.is.s.u-tokyo.ac.jp/~seta/paper/senior_thesis/seniorthesis.pdf)).
