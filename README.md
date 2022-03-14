# matlab_dss_wrapper

Wrapper functions for easier dynamic optimization implementations by using direct single shooting method.

the dss_solve.m automatically computes the derivative of the objective function numerically by using the complex-step derivative approximation (CSDA). If we don't provide the derivative of the objective function, MATLAB actually approximates the derivative numerically. However, this process often fails when the objective function has a high complexity.

CSDA is signifficantly better. It reduces  the number of  iterations for the optimization process to complete. In some cases it even works for longer time horizon.

----------------
The results below is taken from ex1.m. The example is about a moving mass that moves from x = 0 to x = 0.5 within exactly 1 second. The objective is to minimize the effort. At the destination, the mass must stop moving. The solver is SQP.

In these two figures, we can see that with CSDA, SQP can handle long horizon (N=101).

#### Without CSDA - SQP fails 
![](https://github.com/auralius/matlab_dss_wrapper/blob/main/docs/long_horizon_no_derivative.png)

#### With CSDA - SQP still works 
![](https://github.com/auralius/matlab_dss_wrapper/blob/main/docs/long_horizon_with_derivative.png)

----------------

In these two figures, we can see that with CSDA, SQP converges paster (N=10).

#### Without CSDA - SQP converges after 27 iterations
<img src="https://github.com/auralius/matlab_dss_wrapper/blob/main/docs/no_derivative.png" width="514" height="493">

#### With CSDA - SQP converges after 9 iterations
<img src="https://github.com/auralius/matlab_dss_wrapper/blob/main/docs/with_derivative.png" width="518" height="267">


