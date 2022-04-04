function dss = dss_solve(dss)
% Solve an optimal control problem by using direct single shooting method
%
% Syntax:  dss = dss_solve(dss)
%
% Inputs:
%    dss - The data structure for the optimal control problem
%
% Outputs: 
%    dss - The same data structure with optimal solutions stored in
%          dss.hires_sol and dss.lores_sol 
%
% Author:
%   Auralius Manurung
%   Universitas Pertamina
%   auralius.manurung@ieee.org

dss = sanity_check(dss);

if dss.error
    return
end

dss = prepare(dss);

target = @(U)(objfunc_runner(U, dss));

if strcmp(dss.optsolver, 'sqp')
    opts = optimoptions(@fmincon, ...
        'Display', dss.display,...
        'Algorithm', 'sqp', ...
        'SpecifyObjectiveGradient', true, ...
        'StepTolerance', 1e-10);
    U_opt = fmincon(target, dss.intial_guesses, [], [], [], [], ...
        dss.lb, dss.ub, [], opts);

elseif strcmp(dss.optsolver, 'ps')
    opts = optimoptions('patternsearch', ...
        'Display', dss.display, ...
        'UseParallel', dss.parallel);
    U_opt = patternsearch(target, dss.intial_guesses, [], [], [], [], ...
        dss.lb, dss.ub, [], opts);

else
    dss.error = 1;
    return;
end

%  hires_sol contains solution sampled with time interval T_dyn
if strcmp(dss.input_type, 'zoh')
    dss.hires_sol = interp1(dss.lores_tvect, U_opt, dss.hires_tvect, ...
                           'previous'); 
elseif strcmp(dss.input_type, 'foh')
    dss.hires_sol = interp1(dss.lores_tvect, U_opt, dss.hires_tvect, ...
                           'linear'); 
end

%  lores_sol contains solution sampled with time interval T_ocp
dss.lores_sol = U_opt;  
end

%%
function [J, grad_J] = objfunc_runner(U, dss)

% Dynamic simulation
tspan = 0:dss.T_ocp:dss.tf;

if strcmp(dss.odesolver, 'ode23')
    [~, X] = ode23(@(t,x)rhs(t, x, U), tspan, dss.ic);
else
    [~, X] = ode45(@(t,x)rhs(t, x, U), tspan, dss.ic);
end

if length(tspan) ~= length(X)
    J = 1e10; % ode solver fails, apply a very large cost
else
    J = dss.obj_fn(U, transpose(X), dss.T_ocp);
end

if nargout > 1 % gradient required
      
    % Using the Complex-Step Derivative Approximation method
    h = 1e-3; % a small number to perform perturbation
    ih = 1i*h;    
    grad_J = zeros(length(U), 1);

    for k = 1 : length(U)
        U_ = U;
        U_(k) = U(k) + ih;

        % Do the perturbation
        if strcmp(dss.odesolver, 'ode23')
            [~, X_] = ode23(@(t,x)rhs(t, x, U_), tspan, dss.ic);
        else
            [~, X_] = ode45(@(t,x)rhs(t, x, U_), tspan, dss.ic);
        end
              
        J_ = dss.obj_fn(U_, transpose(X_), dss.T_ocp);       
        grad_J(k) = imag(J_)/h;
    end % end for
    
end % end if

%--------------------------------------------------------------------------
function dxdt = rhs(t, x, u)
    if strcmp(dss.input_type, 'zoh')
        un = interp1(dss.lores_tvect, u, t, 'previous');
    elseif strcmp(dss.input_type, 'foh')
        un = interp1(dss.lores_tvect, u, t, 'linear');
    end

    dxdt = dss.state_update_fn(un, x, t);    
end
%--------------------------------------------------------------------------

end

%%
function dss = prepare(dss)

dss.tf          = (dss.n_horizon-1)*dss.T_ocp; % Final time
dss.lores_tvect = 0 : dss.T_ocp : dss.tf;   % Low-resolution time vector
dss.hires_tvect = 0 : dss.T_dyn : dss.tf;   % High-resolution time vector

end

%%
function dss = sanity_check(dss)

dss.error = 0;

% Mandatory fields --------------------------------------------------------
fields = {'n_horizon', ...
    'T_ocp', ...
    'n_inputs',...
    'n_states',...
    'lb',...
    'ub',...
    'T_dyn',...
    'intial_guesses',...
    'obj_fn',...
    'state_update_fn',...
    'ic'    
    };

for k=1:length(fields)
    if ~isfield(dss, fields{k})
        fprintf('Field <strong>%s</strong> is missing!\n\n',fields{k});
        dss.error = 1;
        break;
    end
end

% Sampling periods --------------------------------------------------------
if dss.T_dyn > dss.T_ocp
    fprintf(['<strong>T_dyn</strong> must be' ...
        ' SMALLER than <strong>T_ocp</strong>!\n\n']);
    dss.error = 1;
    return;
end

% Important dimensiouns ---------------------------------------------------
if (~isequal(size(dss.lb),[1, dss.n_horizon])|| ...
    ~isequal(size(dss.ub),[1, dss.n_horizon])|| ...
    ~isequal(size(dss.intial_guesses),[1, dss.n_horizon]))

    fprintf(['These fields must be: [1 x n_horizon] vectors:' ...
        ' <strong>lb, ub, intial_guesses</strong>!\n\n']);
    dss.error = 1;
    return;
end

% Non-mandatory fields ----------------------------------------------------
if ~isfield(dss, 'parallel')
    dss.parallel = false;
end

if ~isfield(dss, 'display')
    dss.display = 'iter';
end

if ~isfield(dss, 'optsolver')
    dss.solver = 'sqp';
end

if ~isfield(dss, 'odesolver')
    dss.solver = 'ode45';
end

if ~isfield(dss, 'input_type')
    dss.input_type = 'foh';
end

end