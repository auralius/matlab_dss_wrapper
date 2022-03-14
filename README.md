# matlab_dss_wrapper

Wrapper functions for easier dynamic optimization implementations by using direct single shooting method.

There are two solvers applied in the dss_solve.m: SQP and pattern search. 

For SQP, dss_solve.m automatically computes the derivative of the objective function by using the complex-step derivative approximation (CSDA). If we don't provide the derivative of the objective function, MATLAB actually approximates the derivative numerically. However, this process often fails when the objective function is a function of the state variables and not the decission variables. 

When the objective function is a function of the state variables, the system's solution must be calculated first. For a complex system, such as a system modeled with PDE, SQP will almost certainly fail if the derivative function is not explicitly provided by the user.

CSDA signifficantly contributes to SQP's improved performances as it reduces  the number of iterations. In some cases, SQP even works for a longer time horizon when using CSDA.

----------------
The results below is taken from ex1.m. The example is about a moving mass that moves from x = 0 to x = 0.5 within exactly 1 second. The objective is to minimize the effort. At the destination, the mass must stop moving. The solver is SQP.

In these two figures, we can see that with CSDA, SQP can handle long horizon (N=101).

#### Without CSDA - SQP fails 
![](https://github.com/auralius/matlab_dss_wrapper/blob/main/docs/long_horizon_no_derivative.png)

#### With CSDA - SQP still works 
![](https://github.com/auralius/matlab_dss_wrapper/blob/main/docs/long_horizon_with_derivative.png)

----------------

In these two figures, we can see that with CSDA (N=10), SQP takes less iterations to converge.

#### Without CSDA - SQP converges after 27 iterations
<img src="https://github.com/auralius/matlab_dss_wrapper/blob/main/docs/no_derivative.png" width="514" height="493">

#### With CSDA - SQP converges after 9 iterations
<img src="https://github.com/auralius/matlab_dss_wrapper/blob/main/docs/with_derivative.png" width="518" height="267">


