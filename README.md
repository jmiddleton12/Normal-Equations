---
Name: Joshua Middleton
Topic: [Topic 9]
Title: Normal Equations
----

# Normal Equations
The normal equations are the closed form solution to least squares problems. They are a system of linear equations used to obtain the least squares solution to an overdetermined linear system. This report will cover the background of linear least squares and data fitting, how the normal equations are derived and used, the drawbacks of the normal equations, applications for least squares problems and the normal equations, and the history of the normal equations.
$$A^TAx = A^Tb$$

***
1. [Linear Least Squares](#linear-least-squares)
2. [Deriving the Normal Equations](#deriving-the-normal-equations)
3. [Solving the Normal Equations](#solving-the-normal-equations)
4. [Conditionging and Error](#conditioning-and-error)
5. [History and Applications](#history-and-applications)
6. [References](#references)
***

## Linear Least Squares
The normal equations were born from the linear least squares problem. This linear algebra problem arises when the solution of an overdetermined linear system is desired. That is a system of the standard $Ax = b$ form where A is n x m and n > m. In this type of system, the vector b is not contained in the image of A. So, it is impossible to find a solution x that exactly satisfies b. Therefore, finding an x which makes Ax as close as possible to b is the "best" solution. This least squares solution can be denoted as x*.

![Linear least squares concept](https://user-images.githubusercontent.com/119821953/205551356-4b81ba66-4e14-450c-8d5d-39e1878b9b42.PNG)
Figure 1: Linear Least Squares Explanation

This figure is a representation of the least squares problem. In this case, the image of A is the subspace in R^2 represented by a plane, but the b vector is in R^3. This problem can be expanded to higher dimensions, but the issue remains the same.

A good way to intuitively understand the problem is to think of the linear regression of a line given a set of data points in 2D space. We know our solution vector x* will contain 2 parameters to describe a line, and the A matrix and b vector contain the information about our given points (The exact setup of this system is shown below). Assuming all the data points aren't the same and don't lie on the same line, it is of course not possible to find an equation of a line that goes through all the data points. So, the least squares solution we want to obtain will be the solution x* which describes the equation of the line which minimizes the distance between itself and all the data points.

![scatterplotpoints](https://user-images.githubusercontent.com/119821953/205553729-20b1a550-9eb4-46cb-80f5-312175b207d0.PNG)
Figure 2: Linear Regression Example

This figure shows a linear least squares solution to a linear regression problem involving 3 data points. The solution line does not pass through any individual point, but the distance between the line and all the points is as small as possible.

Setting up linear least squares problems is done simply by writing a system of equations where each equation is written in terms of the x_i parameters for each data point depending on what type of data is being fit (linear, quadratic, cubic, etc.). t will be used as the independent variable of each data point and y is the dependent variable of each data point. In general, each equation in the system is of the form: 
$$f(t,x) = x_1 \phi_1(t) + x_2 \phi_2(t) + ... + x_i \phi_i(t)$$
Where $x_i$ are the parameters of the fit equation and $\phi_i$ are functions of t depending on the type of fit and are therefore known for each data point. Because these equations are linear in x, the system can be represented in matrix form of $Ax = b$. Of course, in a case where the number of data points being analyzed is equal to the number of parameters (degree) of the fit function, the A matrix is n x n and the problem is not over determined. That is not a least squares problem and will have a unique solution which passes through all the data points in the given set. In a case where the number of data points being analyzed is less than the number of parameters, infinite solutions exist which exactly pass through all the data points. Finally, the case where the number of data points being analyzed is greater than the number of parameters is the least squares problem. There are zero exact solutions as described previously.

As promised, the setup of a least squares problem for a linear regression is as follows. There will be 2 x parameters and therefore 2 $\phi(t)$ functions. Conceptually, $x_1$ will be the y intercept of the solution line and $x_2$ is the slope. Therefore, intuitively, each equation in our system is written as $y = x_2t + x_1$. It is now clear what the $\phi(t)$ functions are. $\phi_1(t) = 1$ and $\phi_2(t) = t$. This equation is written for each of the given data points in the set $(t_n,y_n)$ using the matrix form $Ax = b$ where A contains the $\phi_i(t)$ functions given t, b is the vector of corresponding $y_n's$, and x is the vector of unknown parameters.

![LinearRegressionsetup](https://user-images.githubusercontent.com/119821953/205559521-fdca5a1f-d5fe-40d3-b5f8-c071b6a30318.PNG)
Figure 3: Linear Least Squares Geometry

This results in the system shown in the figure where A is n x 2 because there are n data points and 2 x parameters for a line. This same process can be used for any type of fit function, but the $\phi_i(t)'s$ will change and there may be more of them. For example, a quadratic fit function uses 3 $\phi_i(t)$ functions and 3 x parameters.

## Deriving the Normal Equations
Solving the system $Ax = b$ where A is n x m and n > m requires a least squares approach because there is no exact solution. So, it is not possible to use standard linear algebra approaches such as inverting A (A is not square) or gaussian elimination. That is where the normal equations come in. They will allow the "best" solution to be found given all the information organized into the A matrix and b vector. Two derivation methods will be discussed, the first is more direct and is based on minimizing the residual vector of the system. The second uses geometric reasoning from linear algebra and gives a better visualization of what is happening.

The residual vector is just the difference between $Ax$ and $b$, $r = b - Ax$. The square of the residual vector is to be minimized in order to find the "best" or least squares solution. The magnitude of r is what is important. So, direction and sign doesn't matter. Therefore, taking the square of the euclidean norm of the residual vector and minimizing that is the focus of this approach. Recall that the euclidean norm squared can be written in the form of an inner product of a vector with itself. 
$$||r||^2 = r^Tr = (b - Ax)^T(b - Ax)$$
$$= b^Tb - 2x^TA^Tb + x^TA^TAx$$
To minimize this equation, take the derivative and set it equal to zero:
$$2A^TAX - 2A^Tb = 0$$
Which simplifies to our normal equations:
$$A^TAx = A^Tb$$

The alternative derivation still seeks to minimize the residual vector $b- Ax$. As stated above, b is not in the image of A. Therefore, neither are any of the possible residual vectors for any x in $R^m$. Geometrically, it can be seen that the residual vector will be minimized when it is orthogonal to the image of A:

![GeometricDerivation](https://user-images.githubusercontent.com/119821953/205566000-c1fc6b63-8e32-4cb4-aa7f-6a5cb576276a.PNG)

Knowing that $r^* = b - Ax* $ is orthogonal to the image of A is helpful because that means r* is contained in the orthogonal complement of the image of A. The orthogonal complement of the image of a matrix is equal to the kernel of the transpose of that matrix. Therefore, $r* $ must be contained in the kernel of $A^T$ which means $A^Tr* = 0$ leading to this equation:
$$A^T(b - Ax) = 0$$ 
Which can be simplified to again give the normal equations:
$$A^TAx = A^Tb$$

Finding the x which satisfies this equation will be the one which minimizes the square of the magnitude of the residual vector, the least square solution, x*. And because $A^TA$ is square, this system has a unique solution if A has full column rank (rank(A) = m). In the standard cases where the data points are all unique and the $\phi_i(t)'s$ are formulated properly, the column rank will be full and the normal equations will have a unique solution.

## Solving the Normal Equations
The normal equations can be solved using matrix inversion because $A^TA$ will be full column rank and m x m square (ie. invertible). This is how most classical introductory linear algebra textbooks teach the topic because the concepts of orthogonality, geometry and rank are more of a focus than utilizing the normal equations to solve linear least squares problems themselves. $(A^TA)^{-1}A^T$ is known as the Moore-Penrose pseudoinverse denoted by $A^{+}$.
$$x* = (A^TA)^{-1}A^Tb = A^{+}b$$

In the context of numerical analysis, matrix inversion is almost never used because it is inefficient and unstable. A useful property of $A^TA$ is that it is symmetric positive definite. This makes it a candidate for Cholesky factorization which is the prefered solution method for the normal equations. Cholesky factorization is an algorithm for reducing a symmetric, positive definite matrix into a lower triangular matrix times its conjugate transpose:
$$P = LL^{*}$$

In the case of linear least squares with data and fitting functions that are real numbers, this conjugate transpose becomes simply the transpose. The resulting factorization being a lower triangular real matrix times its transpose which is an upper triangular real matrix:
$$A^TA = LL^T $$

Once A^TA is factored using Cholesky factorization, the system can be solved in two steps using forward and backward substitution such as with other matrix decomposition methods:

$$Lz = A^Tb$$ where z is an intermediate result vector.

$$L^Tx = z$$ where x is our x* that corresponds to the least squares solution.

The Cholesky algorithm can be written as a sequence of modifying the A matrix and forming L_n matrices which are multiplied successively to ultimately compute L after n steps. Keep in mind that in this general form, A is an arbitrary symmetric, positive definite matrix. In the case of the normal equations, the A matrix we are factoring is actually $A^TA$. Another note about the notation below is that $b_i$ refers to the vector which is a portion of the $A^i$ matrix rather than the b vector which is on the right side of the norm al equations, and B^i is a block matrix which is a sub matrix of $A^i$ reffering to the remaining part of the $A^i$ matrix where $A^i_{(j,k)}; j>i , k>i$. In an n x n matrix, i = 1, 2,...,n.

![Cholesky](https://user-images.githubusercontent.com/119821953/205804011-b7727767-d608-4b6d-81e2-01e55dcb6f33.PNG)

The algorithm starts with an A matrix, at each step, $L_i$ is computed and then A is modified before going to the next step. The new A is used to form the next $L_i$ for a total of n steps. This algorithm is written in a way where only the A and L matrix need to be stored to be more memory efficient. It takes about $n^3/3$ floating point operations to compute the Cholesky factorization which is about half the cost of LU factorization which is used for arbitrary square matrices. 

To summarize, solving the normal equations using Cholesky decomposition is done by first computing $A^TA$ and $A^Tb$. Then, the Cholesky algorithm is used to factor $A^TA$ into a form which can be solved using forward and backward substitution resulting in a solution vector, $x*$ which is the least squares solution.

## Conditioning and Error
The discussion of conditioning and error are where the drawbacks of the normal equations are revealed. Up until now it seems like the normal equations method is great. It is closed form and relatively easy to code or even compute by hand for a small n dimension problem, and it makes sense geometrically to directly use the properties of orthogonality in linear algebra to compute a solution directly. However, conditioning and error associated with the normal equations show why it is rarely used in practice for solving least squares problems, and gradient descent or other optimization methods are more practical ways of approaching least squares problems.

The biggest flaw with the normal equations is that their conditioning is inherently poor. $A^TA$ is essentially squaring a matrix which is not ideal for conditioning:
$$cond(A^TA) = [cond(A)]^2$$
So, before even using any numerical method to solve the $A^TA$ matrix, precision is lost simply by forming $A^TA$. This makes intuitive sense if you think about the amount of information that is lost. For example, in the case of a linear regression using 10 data points, the A matrix will have a column of 1's and then a column of 10 numbers corresponding to the independent variable of each data point. Computing $A^TA$ from that 10 x 2 matrix will result in a 2 x 2 matrix totaling 4 parameters. So, the information which once was expressed using 10 independent parameters is now only 4 parameters. It is obvious that precision of the final answer will be reduced if the first step is already this detrimental.

As a general rule, for $K = cond(A)$ about $log_{10}(K) digits of precision are lost. So, squaring the condition number means the number of digits lost is doubled. That could result in having almost no precision in your final answer even if the method used to solve is very accurate. 

Another major issue with normal equations are their susceptibility to rounding error in finite precision. If entries in A are below the machine precision $\epsilon_{mach}$, then computing $A^TA$ can result in a singular matrix. As discussed above, even when numbers that small aren't being dealt with, computing $A^TA$ still makes the problem more ill-conditioned and rounding error in finite precision further increases the error.

## History and Applications
The history of the normal equations follows closely with the history of least squares methods. Before more advanced numerical methods were developed and computers became available, the normal equations were a practical way to solve least squares problems. Least squares problems are used mainly for curve fitting in a variety of fields. Statistics, astronomy, and engineering to name a few. The type of fitting can be single or multi variable and linear and nonlinear formulations exist (linear was covered in this report). The development of least squares analysis and the normal equations was born in the field of astronomy and geodesy.

For centuries, humans desired to fit curves to the shape of the Earth, Moon, and other celestial bodies, and they wanted to find functions to describe the orbital trajectories of the Moon, Earth, and other celestial bodies. Simpler curve fitting methods were used before Marquis Pierre Simon de Laplace started experimenting with more advanced curve fitting methods he created during the late 18th century. He wished to fit an ellipse to a polar cross section of the Earth to more accurately describe the Earth's shape. The main issues with his methods were that the possible error exceeded the accuracy of his instruments, and for overdetermined systems, his methods generated multiple solutions instead of a single "best" solution.

There is some dispute over whether it was Adrien-Marie Legendre or Carl Friedrich Gauss invented linear least squares and the normal equations. However, it was Legendre who first published them in 1805. It is suspected that Guass had already been using least squares methods between 1800 and 1805 during his study of the dwarf planet Ceres. Guass wanted to fit a curve to the orbit of Ceres based on his observations in order to study how the orbital dynamics of one celestial body is influenced by the gravity of another. In his case, he wanted to know how Jupiter's gravity perturbed the orbital trajectory of Ceres. Guass later published his Theoria Motus in 1809 which contained more detailed types of Least Squares analysis

## References

Bretscher, O. (2019). Linear algebra with applications. Pearson India Education Services.

Forbes, E. G. (1978). The astronomical work of Carl Friedrich Gauss (1777–1855). Historia Mathematica, 5(2), 167–181. https://doi.org/10.1016/0315-0860(78)90048-4

Kwiatkowski, R. (2020, December 4). Performing linear regression using the normal equation. Medium. Retrieved December 4, 2022, from https://towardsdatascience.com/performing-linear-regression-using-the-normal-equation-6372ed3c57

Lawson, C. L., &amp; Hanson, R. J. (1974). Solving least squares problems. Prentice-Hall.

Taboga, Marco (2021). "Normal equations", Lectures on probability theory and mathematical statistics. Kindle Direct Publishing. Online appendix. https://www.statlect.com/glossary/normal-equations.

Nievergelt, Y. (2000). A tutorial history of least squares with applications to astronomy and geodesy. Journal of Computational and Applied Mathematics, 121(1-2), 37–72. https://doi.org/10.1016/s0377-0427(00)00343-5

