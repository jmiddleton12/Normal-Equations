---
Name: Joshua Middleton
Topic: [Topic 9]
Title: Normal Equations
----

# Normal Equations
The normal equations are the closed form solution to least squares problems. They are a system of linear equations used to obtain the least squares solution to an overdetermined linear system. This report will cover the background of linear least squares and data fitting, how the normal equation are derived and used, the drawbacks of the normal equations, applications for least squares problems and the normal equations, and the history of the normal equations.
$$A^TAx = A^Tb$$

## Linear Least Squares
The normal equations were born from the linear least squares problem. This linear algebra problem arises when the solution of an overdetermined linear system is desired. That is a system of the standard $Ax = b$ form where A is n x m and n > m. In this type of system, the vector b is not contained in the image of A. So, it is impossible to find a solution x that exactly satisfies b. Therefore, finding an x which makes Ax as close as possible to b is the "best" solution. This least squares solution can be denoted as x*.

![Linear least squares concept](https://user-images.githubusercontent.com/119821953/205551356-4b81ba66-4e14-450c-8d5d-39e1878b9b42.PNG)

This figure is a representation of the least squares problem. In this case, the image of A is the subspace in R^2 represented by a plane, but the b vector is in R^3. This problem can be expanded to higher dimensions, but the issue remains the same.

A good way to intuitively understand the problem is to think of the linear regression of a line given a set of data points in 2D space. We know our solution vector x* will contain 2 parameters to describe a line, and the A matrix and b vector contain the information about our given points (The exact setup of this system is shown below). Assuming all the data points aren't the same and don't lie on the same line, it is of course not possible to find an equation of a line that goes through all the data points. So, the least squares solution we want to obtain will be the solution x* which describes the equation of the line which minimizes the distance between itself and all the data points.

![scatterplotpoints](https://user-images.githubusercontent.com/119821953/205553729-20b1a550-9eb4-46cb-80f5-312175b207d0.PNG)

This figure shows a linear least squares solution to a linear regression problem involving 3 data points. The solution line does not pass through any individual point, but the distance between the line and all the points is as small as possible.

Setting up linear least squares problems is done simply by writing a system of equations where each equation is written in terms of the x_i parameters for each data point depending on what type of data is being fit (linear, quadratic, cubic, etc.). t will be used as the independent variable of each data point and y is the dependent variable of each data point. In general, each equation in the system is of the form: 
$$f(t,x) = x_1 \phi_1(t) + x_2 \phi_2(t) + ... + x_i \phi_i(t)$$
Where $x_i$ are the parameters of the fit equation and $\phi_i$ are functions of t depending on the type of fit and are therefore known for each data point. Because these equations are linear in x, the system can be represented in matrix form of $Ax = b$. Of course, in a case where the number of data points being analyzed is equal to the number of parameters (degree) of the fit funtion, the A matrix is n x n and the problem is not over determined. That is not a least squares problem and will have a unique solution which passes through all the data points in the given set. In a case where the number of data points being analyzed is less than the number of parameters, infinite solutions exist whihc exactly pass through all the data points. Finally, the case where the number of data points being analyzed is greater than the number of parameters is the least squares problem. There are zero exact solutions as described previously.

As promised, the setup of a least squares problem for a linear regression is as follows. There will be 2 x parameters and therefore 2 $\phi(t)$ functions. Conceptually, $x_1$ will be the y intercept of the solution line and $x_2$ is the slope. Therefore, intuitively, each equation in our system is written as $y = x_2t + x_1$. It is now clear what the $\phi(t)$ functions are. $\phi_1(t) = 1$ and $\phi_2(t) = t$. This equation is written for each of the given data points in the set $(t_n,y_n)$ using the matrix form $Ax = b$ where A contains the $\phi_i(t)$ functions given t, b is the vector of corresponding $y_n's$, and x is the vector of unknown parameters.

![LinearRegressionsetup](https://user-images.githubusercontent.com/119821953/205559521-fdca5a1f-d5fe-40d3-b5f8-c071b6a30318.PNG)

This results in the system shown in the figure where A is n x 2 because there are n data points and 2 x parameters for a line. This same process can be used for any type of fit funtion, but the $\phi_i(t)'s$ will change and there may be more of them. For example, a quadratic fit function uses 3 $\phi_i(t)$ functions and 3 x parameters.

## Deriving the Normal Equations
Solving the system $Ax = b$ where A is n x m and n > m requires a least squares approach because there is no exact solution. So, it is not possible to use standard linear algebra approaches such as inverting A (A is not square) or gaussian elimination. That is where the normal equations come in. They will allow the "best" solution to be found given all the information organized into the A matrix and b vector. Two derivation methods will be discussed, the first is more direct and is based on minimizing the residual vector of the system. The second uses geometric reasoning from linear algebra and gives a better visualization of what is happening.

The residual vector is just the difference between $Ax$ and $b$, %r = b - Ax$. The square of the residual vector is to be minimized in order to find the "best" or least squares solution. The magnitude of r is what is important. So, direction and sign doesn't matter. Therefore, taking the square of the euclidean norm of the residual vector and minimizing that is the focus of this approach. Recall that the euclidean norm squared can be written in the form of an inner product of a vector with itself. 
$$||r||^2 = r^Tr = (b - Ax)^T(b - Ax)$$
$$= b^Tb - 2x^TA^Tb + x^TA^TAx$$
To minimize this equation, take the derivative and set it equal to zero:
$$2A^TAX - 2A^Tb = 0$$
Which simplifies to our normal equations:
$$A^TAx = A^Tb$$

The alternative derivation still seeks to minimize the residual vector $b- Ax$. As stated above, b is not in the image of A. Therefore, neither are any of the possible residual vectors for any x in $R^m$. Geometrically, it can be seen that the residual vector will be minimized when it is orthogonal to the image of A:

![GeometricDerivation](https://user-images.githubusercontent.com/119821953/205566000-c1fc6b63-8e32-4cb4-aa7f-6a5cb576276a.PNG)

Knowing that $r^* = b - Ax* $ is orthogonal to the image of A is helpful because that means r* is contained in the orthogonal complement of the image of A. The orthogonal complement of the image of a matrix is equal to the kernal of the transpose of that matrix. Therefore, $r* $ must be contained in the kernal of $A^T$ whihc means $A^Tr* = 0$ leading to this equation:
$$A^T(b - Ax) = 0$$ 
Which can be simplified to again give the normal equations:
$$A^TAx = A^Tb$$

Finding the x which satisifes this equation will be the one which minimizes the square of the magnitude of the residual vector, the least square solution, x*. And because $A^TA$ is square, this system has a unique solution if A has full column rank (rank(A) = m). In the standard cases where the data points are all unique and the $\phi_i(t)'s$ are formulated properly, the column rank will be full and the normal equations will have a unique solution.

## Solving the Normal Equations
The normal equations can be solved using matrix inversion because $A^TA$ will be full column rank and m x m square (ie. invertible). This is how most classical introductory linear algebra textbooks teach the topic because the concepts of orthogonality, geometry and rank are more of a focus than utilizing the normal equations to solve linear least squares problems themselves. $(A^TA)^{-1}A^T$ is known as the Moore-Penrose pseudoinverse denoted by $A^{+}$.
$$x* = (A^TA)^{-1}A^Tb = A^{+}b$$

In the context of numerical analysis, matrix inversion is almost never used because it is inefficient and unstable. A useful property of $A^TA$ is that it is symmetric positive definite. This makes it a candidate for Cholesky factorization which is the preffered solution method for the normal equations. 



