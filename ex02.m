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
dss.n_horizon        = N;
dss.T_ocp            = T_ocp;         % optimal control problem's period

dss.n_inputs         = 1;             % valve closing rate
dss.n_states         = 2*pipe.M;      % pressure and velocity at every pipeline node
dss.lb               = zeros(N,1);    % fully-open
dss.ub               = ones(N,1);     % fully-close
dss.intial_guesses   = 0.5*ones(N,1); % half-open / half-close
dss.T_dyn            = 0.00001;       % dynamic simulation's period

dss.obj_fn           = @obj_fn;
dss.state_update_fn  = @state_update_fn;
dss.ic               = [pipe.p0 pipe.v0];
dss.parameterization = 'foh';

% Optional fields ---------------------------------------------------------
dss.parallel         = 'always';
dss.display          = 'iter';
dss.solver           = 'ps';

% Run the solver ----------------------------------------------------------
dss = dss_solve(dss);
dss = dss_resimulate(dss);

%%
function J = obj_fn(U, X, dt)

pipe = create_pipe();

gamma   = 2;
T       = 10;
p_ref   = 2e5;
p_toll  = 1e4;

p = X(:,1:pipe.M);

% Intial 
delta0 = ((p(:, 1) - p_ref) ./ p_toll) .^ (2*gamma);
Sum0 = 1/T*sum(delta0)*dt;

% Terminal
delta1 = ((p(:, end) - p_ref) ./ p_toll) .^ (2*gamma);
Sum1 = 1/T*sum(delta1)*dt;

% Stage
delta2 = ((p - p_ref) ./ p_toll).^ (2*gamma);
Sum2 = 1/(pipe.L*T)*sum(sum(delta2))*pipe.dl*dt;

J = Sum0 + Sum1 + Sum2;

end

%% The state update funtion 
function X = state_update_fn(U, X, dt)
pipe = create_pipe();

p_dot = zeros(1, pipe.M);
v_dot = zeros(1, pipe.M);

p = X(1:pipe.M);
v = X(pipe.M+1:end);

i = 2 : pipe.M-1;
v_dot(i) = -1/pipe.rho * (p(i+1)-p(i-1)) / (2*pipe.dl) - ...
           pipe.f/(2*pipe.D) .* abs(v(i)).*v(i);
p_dot(i) = -pipe.rho*pipe.c^2 * (v(i+1)-v(i-1)) / (2*pipe.dl);

p = p_dot*dt + p;
v = v_dot*dt + v;

% Apply BC
v(end) = pipe.u_max - pipe.u_max * U;  % Apply the inputs                    
p(end) = p(end) + dt * pipe.rho*pipe.c^2/pipe.dl * ( v(end-1)-v(end) );

v(1) = v(1) + dt * ( 1/(pipe.rho*pipe.dl)*(pipe.P-p(2)) - pipe.f/(2*pipe.D)*v(1)*abs(v(1)) );
p(1) = pipe.P;

X(1:pipe.M) = p;
X(pipe.M+1:end) = v;

end

%%
function pipe = create_pipe()
% The parameters
pipe.L = 200;    % m
pipe.D = 100e-3; % m
pipe.rho = 1000; % kg/m3
pipe.c = 1200;   % m/s
pipe.f = 0.03;   % friction factor
pipe.P = 2e5;    % Pa

pipe.m = 16;            % Even integer, number of the pipeline segments
pipe.M = pipe.m + 1;    % Number of the pipeline nodes
pipe.dl = pipe.L / pipe.m;

% maximum velocity of the water inside the pipe
pipe.u_max = 2;

% Initial values for p and v
l = linspace(0, pipe.L, pipe.M);
pipe.p0 = pipe.P - 2 * pipe.rho * pipe.f / pipe.D * l;
pipe.v0 = pipe.u_max * ones(1, pipe.M);  % See Sec. 4 of the referenced paper
end

