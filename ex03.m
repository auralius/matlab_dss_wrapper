% The dynamics model here describes the longitudinal motion of an F-8
% aircraft, taken from the following paper:
% Garrard, W. L., & Jordan, J. M. (1977). Design of nonlinear automatic 
% flight control systems. Automatica, 13(5), 497–505. 
% https://doi.org/10.1016/0005-1098(77)90070-X
%
% or this paper:
% Banks, S. P., & Mhana, K. J. (1992). Optimal control and stabilization 
% for nonlinear systems. IMA Journal of Mathematical Control and 
% Information, 9(2), 179–196. https://doi.org/10.1093/imamci/9.2.179
%
% or this paper:
% Kaya, C. Y., & Noakes, J. L. (2003). Computational Method for Time-
% Optimal Switching Control. Journal of Optimization Theory and 
% Applications, 117(1), 69–92. https://doi.org/10.1023/A:1023600422807
%
% or this paper:
% Teo, K. L., Lee, H. W. J., & Rehbock, V. (1998). Control parametrization 
% enhancing technique for time optimal control and optimal three-valued 
% control problems. Dynamics of Continuous, Discrete and Impulsive Systems 
% Series B: Application and Algorithm, 4(4), 617–631.
%
% -------------------------------------------------------------------------
% Pattern search works better than SQP+CSDA
% -------------------------------------------------------------------------
%
% Auralius Manurung, ME, Universitas Pertamina, Indonesia
%

clear;
close all;
clc;

% Setup the horizon
Tf    = 10;              % 1 second
T_ocp = 0.5;            % Temporal discretization step
t     = 0 : T_ocp : Tf;
N     = length(t);

% Mandatory fields --------------------------------------------------------
dss.n_horizon        = N;
dss.T_ocp            = T_ocp;       % optimal control problem's period

dss.n_inputs         = 1;
dss.n_states         = 3;
dss.lb               = deg2rad(-3)*ones(1,N);
dss.ub               = deg2rad(3)*ones(1,N);
dss.intial_guesses   = zeros(1,N);
dss.T_dyn            = 0.01;        % dynamic simulation's period

dss.obj_fn           = @obj_fn;
dss.state_update_fn  = @state_update_fn;
dss.ic               = [deg2rad(26.7) 0 0];

dss.input_type      = 'zoh'; % zoh or foh?

% Optional fields ---------------------------------------------------------
dss.parallel         = true;
dss.display          = 'iter';
dss.optsolver        = 'ps';
dss.odesolver        = 'ode45';

% Run the solver ----------------------------------------------------------
tic
dss = dss_solve(dss);
toc
dss = dss_resimulate(dss);

%%
function J = obj_fn(U, X, dt)
x_bar = 0.01;
Tf = 10;

J = 1/(Tf/dt+1) * sum( dt*((X(1,:)/x_bar).^4 + (X(2,:)/x_bar).^4 + ...
    (X(3,:)/x_bar).^4) );
end

%% The state update funtion 
function dXdt = state_update_fn(U, X, t)
dXdt = [(-0.877*X(1) + X(3) - 0.088*X(1)*X(3) + 0.47*X(1)^2 - ...
        0.019*X(2)^2 - X(1)^2*X(3) + 3.846*X(1)^3 - 0.215*U + ...      
        0.28*X(1)^2*U + 0.47*X(1)*U^2 + 0.63*U^3);
        X(3);
        (-4.208*X(1) - 0.396*X(3) - 0.47*(X(1)^2) - 3.564*(X(1)^3) - ...
        20.967*U + 6.265*X(1)^2*U + 46*X(1)*U^2 + ...
        61.4*U^3)];
end

