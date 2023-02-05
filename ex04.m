clear;
close all;
clc;


% Setup the horizon
Tf    = 2;              % 1 second
T_ocp = 0.1;            % Temporal discretization step
t     = 0 : T_ocp : Tf;
N     = length(t);

% Mandatory fields --------------------------------------------------------
dss.n_horizon        = N;
dss.T_ocp            = T_ocp;       % optimal control problem's period

dss.n_inputs         = 1;
dss.n_states         = 4;
dss.lb               = -10*ones(1,N);
dss.ub               = 10*ones(1,N);
dss.intial_guesses   = zeros(1,N);
dss.T_dyn            = 0.01;        % dynamic simulation's period

dss.obj_fn           = @obj_fn;
dss.state_update_fn  = @state_update_fn;
dss.ic               = [0 0 0 0];

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
r1 = 0;
r2 = 0;
r3 = 100;
r4 = 100;

% Final state
xf = [0; 0; 1.0; 0];

terminal_cost = r1*(X(1,end)-xf(1)).^2 + ...
                r2*(X(2,end)-xf(2)).^2 + ...
                r3*(X(3,end)-xf(3)).^2 + ...
                r4*(X(4,end)-xf(4)).^2; 
J = 0*sum(U.^2) + sum(((X(3,:)-xf(3))/1e-3).^4) + ...
    sum(((X(4,:)-xf(4))/1e-2).^4); 
%J = 10*dt*(sum(U.^2) + sum(X(1,:).^2) + sum(X(2,:).^2)) + ...
%    terminal_cost;

end

%% The state update funtion 
function Xdot = state_update_fn(U, X, t)

l = 1.0;
g = 9.8;
B = 2.0;
M = 1.0;

Xdot = zeros(4,1);

theta = X(1,1);
thetad = X(2,1);
x = X(3,1);
xd = X(4,1);

Xdot(1,1) = thetad;
Xdot(2,1) = -B*thetad/(M*l^2) - g*sin(theta)/l - cos(theta)*U/l;
Xdot(3,1) = xd;
Xdot(4,1) = U;

end

