# Normal Equations
The normal equations are the closed form solution to least squares problems. They are a system of linear equations used to obtain the least squares solution to an overdetermined linear system. This report will cover the background of linear least squares and data fitting, how the normal equation are derived and used, the drawbacks of the normal equations, applications for least squares problems and the normal equations, and the history of the normal equations.
$$A^TAx = A^Tb$$

## Linear Least Squares
The normal equations were born from the linear least squares problem. This linear algebra problem arises when the solution of an overdetermined linear system is desired. That is a system of the standard $Ax = b$ form where A is n x m and n > m. In this type of system, the vector b is not contained in the image of A. So, it is impossible to find a solution x that exactly satisfies b. Therefore, finding an x which makes Ax as close as possible to b is the "best" solution. This least squares solution can be denoted as x*.
![Linear least squares concept](https://user-images.githubusercontent.com/119821953/205551356-4b81ba66-4e14-450c-8d5d-39e1878b9b42.PNG)
This figure is a representation of the least squares problem. In this case, the image of A is the subspace of A is in R^2 represented by a plane, but the b vector is in R^3. This problem can be expanded to higher dimensions, but the issue remains the same.
