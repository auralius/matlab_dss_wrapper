% Water hammer Effect Minimization
%
% Based on the following paper:
% Chen, T., Xu, C., Lin, Q., Loxton, R., & Teo, K. L. (2015). Water hammer 
% mitigation via PDE-constrained optimization. Control Engineering 
% Practice, 45, 54–63. https://doi.org/10.1016/j.conengprac.2015.08.008
%
% Chen, T., Ren, Z., Xu, C., & Loxton, R. (2015). Optimal boundary control 
% for water hammer suppression in fluid transmission pipelines. Computers 
% & Mathematics with Applications, 69(4), 275–290. 
%
% -------------------------------------------------------------------------
%
% Auralius Manurung, ME, Universitas Pertamina, Indonesia
%

clear;
close all;
clc;

pipe = create_pipe();
% Setup the horizon
Tf    = 10;              % 10 seconds
T_ocp = 1;               % Temporal discretization step
t     = 0 : T_ocp : Tf;
N     = length(t);

% Mandatory fields --------------------------------------------------------
dss.n_horizon           = N;
dss.T_ocp               = T_ocp;       % optimal control problem's period

dss.n_inputs            = 1;           % valve closing rate
dss.n_states            = 2*pipe.M;    % pressure and velocity at each node

dss.lb                  = 0:1/(N-1):1; % diagonal line
dss.ub                  = ones(1,N);   % fully-close
dss.intial_guesses      = dss.lb;      % initializa with the lower bounds

dss.lb(1)               = 0;           % fully closed at t=0
dss.ub(1)               = 0;
dss.intial_guesses(1)   = 0;

dss.lb(end)             = 1;           % fully open at t = 10
dss.ub(end)             = 1;
dss.intial_guesses(end) = 1;

dss.T_dyn               = 0.001;       % dynamic simulation's period
dss.obj_fn              = @obj_fn;
dss.state_update_fn     = @state_update_fn;
dss.ic                  = [pipe.p0 pipe.v0];

% Optional fields ---------------------------------------------------------
dss.input_type          = 'foh'; % zoh or foh?
dss.parallel            = 'always'; % Only works for pattern search
dss.display             = 'iter';
dss.optsolver           = 'sqp';
dss.odesolver           = 'ode23';  % ODE23 as proposed by the paper

% Run the solver ----------------------------------------------------------
tic
dss = dss_solve(dss);
toc
dss = dss_resimulate(dss);

close all; % Too many states for the automatic plotting!

figure
plot(dss.hires_tvect, dss.hires_states(pipe.M,:))
xlabel('Time');
ylabel('Pressure');

figure
plot(dss.hires_tvect, dss.hires_sol)
xlabel('Time');
ylabel('Valve closing rate');

%%
function J = obj_fn(~, X, dt)

pipe = create_pipe();

gamma   = 2;
T       = 10;
p_ref   = 2e5;
p_toll  = 1e5;

p = X(1:pipe.M,:);

% Terminal
delta1 = ((p(end, :) - p_ref) ./ p_toll) .^ (2*gamma);
Sum1 = 1/T*sum(delta1)*dt;

% Stage
delta2 = ((p - p_ref) ./ p_toll).^ (2*gamma);
Sum2 = 1/(pipe.L*T)*sum(sum(delta2))*pipe.dl*dt;

J = Sum1 + Sum2;

end

%% The state update funtion 
function dxdt = state_update_fn(U, X, ~)

pipe = create_pipe();

p_dot = zeros(pipe.M, 1);
v_dot = zeros(pipe.M, 1);

p = X(1:pipe.M);
v = X(pipe.M+1:end);

% Boundary conditions
p(1) = pipe.P;
v(end) = pipe.u_max - pipe.u_max * U;                   

p_dot(end) = pipe.rho*pipe.csqrd/pipe.dl * (v(end-1)-v(end));
v_dot(1) = 1/(pipe.rho*pipe.dl)*(pipe.P-p(2)) - pipe.f/(2*pipe.D)*v(1)*abs(v(1));

% Dynamics update                
k = 2 : pipe.M-1; % Pipe segments
v_dot(k) = -1/pipe.rho * (p(k+1)-p(k)) / (pipe.dl) - pipe.f/(2*pipe.D) .* ...
    abs(v(k)).*v(k);
p_dot(k) = -pipe.rho*pipe.csqrd * (v(k)-v(k-1)) / (pipe.dl);
         
dxdt = [p_dot; v_dot];

end

%%
function pipe = create_pipe()
% The parameters
pipe.L = 200;    % m
pipe.D = 100e-3; % m
pipe.rho = 1000; % kg/m3
pipe.c = 1200;   % m/s
pipe.csqrd = pipe.c^2;
pipe.f = 0.03;   % friction factor
pipe.P = 2e5;    % Pa

pipe.m = 16;             % Integer, number of the pipeline segments
pipe.M = pipe.m + 1;    % Number of the pipeline nodes
pipe.dl = pipe.L / pipe.m;

% maximum velocity of the water inside the pipe
pipe.u_max = 2;

% Initial values for p and v
l = linspace(0, pipe.L, pipe.M);
pipe.p0 = pipe.P - 2 * pipe.rho * pipe.f / pipe.D * l;
pipe.v0 = pipe.u_max * ones(1, pipe.M);  % Sec. 4 of the referenced paper
end

