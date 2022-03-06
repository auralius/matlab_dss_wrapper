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

if strcmp(dss.solver, 'sqp')
    opts = optimoptions(@fmincon, ...
        'Display', dss.display,...
        'Algorithm', 'sqp', ...
        'UseParallel', dss.parallel);
    U_opt = fmincon(target, dss.intial_guesses, [], [], [], [], ...
        dss.lb, dss.ub, [], opts);

elseif strcmp(dss.solver, 'ps')
    opts = optimoptions('patternsearch', ...
        'Display', dss.display, ...
        'UseParallel', dss.parallel);
    U_opt = patternsearch(target, dss.intial_guesses, [], [], [], [], ...
        dss.lb, dss.ub, [], opts);

else
    dss.error = 1;
    return;
end

dss.hires_sol = upsampling(dss, U_opt); %  sampled with time interval T_dyn
dss.lores_sol = U_opt;                  %  sampled with time interval T_ocp
end

%%
function res = objfunc_runner(U, dss)
N = length(dss.hires_tvect);

X_upsampled = zeros(N, dss.n_states);
U_upsampled = upsampling(dss, U);

% Dynamic simulation
X_upsampled(1,:) = dss.ic; % Apply IC

for k = 1 : N-1
    X_upsampled(k+1,:) = dss.state_update_fn(U_upsampled(k), X_upsampled(k,:), dss.T_dyn);
end

X = X_upsampled((1 : int32(dss.T_ocp/dss.T_dyn) : end), :); % Downsampling
res = dss.obj_fn(U, X, dss.T_ocp);

end

%%
function dss = prepare(dss)

dss.tf          = (dss.n_horizon-1)*dss.T_ocp; % Final time
dss.lores_tvect = (0 : dss.T_ocp : dss.tf)';   % Low-resolution time vector
dss.hires_tvect = (0 : dss.T_dyn : dss.tf)';   % High-resolution time vector

end

%%
function hires = upsampling(dss, lores)

if strcmp(dss.parameterization, 'zoh') == 1
    hires = zoh_upsample_trim(lores, int32(dss.T_ocp/dss.T_dyn), ...
                              length(dss.hires_tvect));
elseif strcmp(dss.parameterization, 'foh') == 1
    hires = interp1(dss.lores_tvect, lores, dss.hires_tvect, 'linear');
end

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
    'ic',...
    'parameterization'
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
if (~isequal(size(dss.lb),[dss.n_horizon,1])|| ...
    ~isequal(size(dss.ub),[dss.n_horizon,1])|| ...
    ~isequal(size(dss.intial_guesses),[dss.n_horizon,1]))

    fprintf(['These fields must be: [n_horizon x 1] vectors:' ...
        ' <strong>lb, ub, intial_guesses</strong>!\n\n']);
    dss.error = 1;
    return;
end

% Mandatory fields --------------------------------------------------------
if ~isfield(dss, 'parallel')
    dss.parallel = false;
end

if ~isfield(dss, 'display')
    dss.display = 'iter';
end

if ~isfield(dss, 'solver')
    dss.solver = 'sqp';
end

end