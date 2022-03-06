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
N = length(dss.hires_tvect);
x(1,:) = dss.ic; % Apply IC

for k = 1 : N-1
    x(k+1,:) = dss.state_update_fn(dss.hires_sol(k), x(k,:), dss.T_dyn);
end

M = dss.n_inputs + dss.n_states;
figure
for k = 1 : dss.n_states    
    subplot(M, 1, k);
    hold on;
    plot(dss.hires_tvect, x(:,k));
    xlabel('Time (s)');
    s = ['State #' num2str(k)];
    ylabel(s);
end

for k = 1 : dss.n_inputs   
    subplot(M, 1, k+dss.n_states);
    hold on;
    plot(dss.hires_tvect, dss.hires_sol(:,k));
    xlabel('Time (s)')
end

dss.opt_states = x;