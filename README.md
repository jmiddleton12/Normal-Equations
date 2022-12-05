# Normal Equations
The normal equations are the closed form solution to least squares problems. They are a system of linear equations used to obtain the least squares solution to an overdetermined linear system. This report will cover the background of linear least squares and data fitting, how the normal equation are derived and used, the drawbacks of the normal equations, applications for least squares problems and the normal equations, and the history of the normal equations.
$$A^TAx = A^Tb$$

## Linear Least Squares
The normal equations were born from the linear least squares problem. This linear algebra problem arises when the solution of an overdetermined linear system is desired. That is a system of the standard $Ax = b$ form where A is n x m and n > m. In this type of system, the vector b is not contained in the image of A. So, it is impossible to find a solution x that exactly satisfies b. Therefore, finding an x which makes Ax as close as possible to b is the "best" solution. This least squares solution can be denoted as x*.
![Linear least squares concept](https://user-images.githubusercontent.com/119821953/205551356-4b81ba66-4e14-450c-8d5d-39e1878b9b42.PNG)
This figure is a representation of the least squares problem. In this case, the image of A is the subspace in R^2 represented by a plane, but the b vector is in R^3. This problem can be expanded to higher dimensions, but the issue remains the same.

A good way to intuitively understand the problem is to think of the linear regression of a line given a set of data points in 2D space. We know our solution vector x* will contain 2 parameters to describe a line, and the A matrix and b vector contain the information about our given points (The exact setup of this system is shown ahead). Assuming all the data points aren't the same and don't lie on the same line, it is of course not possible to find an equation of a line that goes through all the data points. So, the least squares solution we want to obtain will be the solution x* which describes the equation of the line which minimizes the distance between itself and all the data points.
![scatterplotpoints](https://user-images.githubusercontent.com/119821953/205553729-20b1a550-9eb4-46cb-80f5-312175b207d0.PNG)
This figure shows a linear least squares solution to a linear regression problem involving 3 data points. The solution line does not pass through any individual point, but the distance between the line and all the points is as small as possible.

Setting up linear least squares problems is done simply by writing a system of equations where each equation is written in terms of the x_i parameters for each data point depending on hwat type of data is being fit. In general, each equation in the system is of the form: 
$$f(t,x) = x_1 \phi_1 + x_2 \phi_2 + ... + x_i \phi_i$$
