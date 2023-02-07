% A mass-damper system that can move freely in a 2D-space.
% The input is a (2x1) force vector.

clear;
close all;
clc;


% Setup the horizon
Tf    = 1;              % 1 second
T_ocp = 0.1;            % Temporal discretization step
t     = 0 : T_ocp : Tf;
N     = length(t);

n_inputs = 2;
n_states = 4;

% Mandatory fields --------------------------------------------------------
dss.n_horizon        = N;
dss.T_ocp            = T_ocp;       % optimal control problem's period

dss.n_inputs         = n_inputs;
dss.n_states         = n_states;
dss.lb               = -4*ones(n_inputs,N);
dss.ub               = 4*ones(n_inputs,N);
dss.intial_guesses   = 4*ones(n_inputs, N);
dss.T_dyn            = 0.01;        % dynamic simulation's period

dss.obj_fn           = @obj_fn;
dss.state_update_fn  = @state_update_fn;
dss.ic               = zeros(1, n_states);

dss.input_type      = 'foh'; % zoh or foh?

% Optional fields ---------------------------------------------------------
dss.parallel         = true;
dss.display          = 'iter';
dss.optsolver        = 'sqp';
dss.odesolver        = 'ode45';

% Run the solver ----------------------------------------------------------
tic
dss = dss_solve(dss);
toc
dss = dss_resimulate(dss);

%%
function J = obj_fn(U, X, dt)
% Weighting factors for the terminal states
r1 = 100;
r2 = 100;

r3 = 100;
r4 = 100;

% Final state
xf = [0.5; 0; 1.0 ; 0];

J = r1*sum(X(1,end)-xf(1)).^2 + r2*sum(X(2,end)-xf(2)).^2 ...
    + r3*sum(X(3,end)-xf(3)).^2 + r4*sum(X(4,end)-xf(4)).^2 ...
    + dt * (sum(U(1,:).^2) + sum(U(2,:).^2)); 
end

%% The state update funtion 
function dXdt = state_update_fn(U, X, t)

m = 1;   % Mass
b = 0.1; % Damping coefficient

dXdt = [0 1 0 0; 0 -b/m  0 0; 0 0 0 1; 0 0 0 -b/m]*X + ...
    transpose([0 1/m 0 0; 0 0 0 1/m])*U;


end

