## Gaussian Elimination row-oriented algorithm
--- This algorithm first transforms given matrix(system of linear equations) into upper triangular matrix and then do back substitution to get final answer. 
--- Results obtained from this algorithm are accurate.
--- I used pivotting method for converting the augmented matrix to upper triangular.
| Number of Processes | Time (in milliseconds) |
|---------------------|------------------------|
|1                     |1900                        |
|2                     |1024                        |
|5                     |1689                        |
|10                     |1510                      |

## Conjugate Gradient Method
--- It is an iterative method to calculate solution to system of linear equations.
--- The results were accurate upto a certain decimal point, depending on given epsilon value. Epsilon is the tolerance.
--- Results with tolerance of 0.0001(10^-4)
| Number of Processes | Time (in milliseconds) |Performance Increase(w.r.t. Gaussian)|
|---------------------|------------------------|-----------|
|1                     |568                        |3.3|
|2                     |678                        |1.5|
|5                     |669                        |2.5|
|10                     |564                        |2.7|

--- Results with tolerance of 0.0000001(10^-7)
| Number of Processes | Time (in milliseconds) |Performance Increase(w.r.t. Gaussian)|
|---------------------|------------------------|-----------|
|1                     |780                        |2.4|
|2                     |944                        |1.1|
|5                     |867                        |1.9|
|10                     |783                        |1.9|

--- All above time computations are for 1000*1000 matrix.
