# matlab_dss_wrapper

Wrapper functions for easier dynamic optimization implementations by using direct single shooting method.

the dss_solve.m automatically computes the derivative of the objective function numerically by using the complex-step derivative approximation (CSDA). If we don't provide the derivative of the objective function, MATLAB actually approximates the derivative numerically. However, this process often fails when the objective function has a high complexity.

CSDA is signifficantly better. It reduces  the number of  iterations for the optimization process to complete. In some cases it even works for longer time horizon.
