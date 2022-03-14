function dss = dss_resimulate(dss)
% Resimulate the system dynamics by using the found optimal solutions
%
% Syntax:  dss = dss_resimulate(dss)
%
% Inputs:
%    dss - The data structure for the optimal control problem that already
%          contains the solutions
%
% Outputs: 
%    dss - The same data structure wit simulated optimal states in
%          dss.opt_states
%
% Author:
%   Auralius Manurung
%   Universitas Pertamina
%   auralius.manurung@ieee.org

if dss.error == 1
    return
end

% Dynamic simulation
tspan = 0:dss.T_dyn:dss.tf;

if strcmp(dss.odesolver, 'ode23')
    [~,x] = ode23( @rhs, tspan, dss.ic);
else
    [~,x] = ode45( @rhs, tspan, dss.ic);
end

x = x';

M = dss.n_inputs + dss.n_states;
figure
for k = 1 : dss.n_states    
    subplot(M, 1, k);
    hold on;
    plot(dss.hires_tvect, x(k,:));
    xlabel('Time (s)');
    s = ['State #' num2str(k)];
    ylabel(s);
end

for k = 1 : dss.n_inputs   
    subplot(M, 1, k+dss.n_states);
    hold on;
    plot(dss.hires_tvect, dss.hires_sol(k,:));
    xlabel('Time (s)')
end

dss.hires_states = x;
dss.lores_states = x(:, (1 : int32(dss.T_ocp/dss.T_dyn) : end)); % Downsampling

%--------------------------------------------------------------------------
function dxdt = rhs(t,x)
    u = interp1(dss.lores_tvect, dss.lores_sol, t, 'linear');
    dxdt = dss.state_update_fn(u, x);    
end
%--------------------------------------------------------------------------
end