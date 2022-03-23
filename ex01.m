clear;
close all;
clc;


% Setup the horizon
Tf    = 1;              % 1 second
T_ocp = 0.1;            % Temporal discretization step
t     = 0 : T_ocp : Tf;
N     = length(t);

% Mandatory fields --------------------------------------------------------
dss.n_horizon        = N;
dss.T_ocp            = T_ocp;       % optimal control problem's period

dss.n_inputs         = 1;
dss.n_states         = 2;
dss.lb               = -4*ones(1,N);
dss.ub               = 4*ones(1,N);
dss.intial_guesses   = 2*ones(1,N);
dss.T_dyn            = 0.01;        % dynamic simulation's period

dss.obj_fn           = @obj_fn;
dss.state_update_fn  = @state_update_fn;
dss.ic               = [0 0];

dss.input_type      = 'zoh'; % previous or foh?

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
r1 = 250;
r2 = 100;

% Final state
xf = [0.5; 0];

terminal_cost = r1*(X(1,end)-xf(1)).^2 + r2*(X(2,end)-xf(2)).^2; 
J = dt*sum(U.^2) + terminal_cost;
end

%% The state update funtion 
function dXdt = state_update_fn(U, X, t)

m = 1;   % Mass
b = 0.1; % Damping coefficient

dXdt = [0 1; 0 -b/m]*X + [0; 1/m]*U;

end

