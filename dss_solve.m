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

    phi{1} = @objf;
    phi{2} = @jac;

    max_iter = 1000;
    tol = 1e-5;

    % --------------------------------------------------------------------------
    function J = compute_cost(U)
        tspan = 0:dss.T_ocp:dss.tf;

        if strcmp(dss.odesolver, 'ode23')
            [~, X] = ode23(@(t,x)rhs_(t, x, U), tspan, dss.ic);
        elseif strcmp(dss.odesolver, 'ode23s')
            [~, X] = ode23(@(t,x)rhs_(t, x, U), tspan, dss.ic);
        elseif strcmp(dss.odesolver, 'ode45')
            [~, X] = ode45(@(t,x)rhs_(t, x, U), tspan, dss.ic);
        end

        if length(tspan) ~= size(X,1)
            disp("Failed solving the ODE!")
            J = inf;
        else
            U_ = reshape(U, [dss.n_inputs, dss.n_horizon]);
            J = dss.obj_fn(U_, transpose(X), dss.T_ocp);
        end
    end

    % --------------------------------------------------------------------------
    function J = objf(U)
       J = compute_cost(U)
    end

    % --------------------------------------------------------------------------

    function jac = jac(U)
        jac_ = jacobs(U, @compute_cost);
        %jac_ = numgradient("compute_cost", {guess});
        jac = transpose(jac_);
    end

    % --------------------------------------------------------------------------
    function dxdt = rhs_(t, x, u)
        u_ = reshape(u,[dss.n_inputs, dss.n_horizon]);
        if strcmp(dss.input_type, 'zoh')
            %un = interp1(dss.lores_tvect, u_, t, 'previous');
            un = costum_interp1(dss.lores_tvect, u_, t, 'previous');
        elseif strcmp(dss.input_type, 'foh')
            %un = interp1(dss.lores_tvect, u_, t, 'linear');
            un = costum_interp1(dss.lores_tvect, u_, t, 'linear');
        end

        dxdt = dss.state_update_fn(un, x, t);
    end

    % --------------------------------------------------------------------------
    if strcmp(dss.optsolver, 'sqp')
        tic

        [U_opt, obj, info, iter, nf, lambda] = sqp(
        reshape(dss.intial_guesses, [1, dss.n_horizon*dss.n_inputs]), ...
        phi, ...
        [], ...
        [], ...
        reshape(dss.lb, [1, dss.n_horizon*dss.n_inputs]), ...
        reshape(dss.ub, [1, dss.n_horizon*dss.n_inputs]), ...
        max_iter, tol)

        toc
    else
        dss.error = 1;
        return;
    end

    % --------------------------------------------------------------------------
    %  hires_sol contains solution sampled with time interval T_dyn
    U_opt_ = reshape(U_opt, [dss.n_inputs, dss.n_horizon]);

    if strcmp(dss.input_type, 'zoh')
        dss.hires_sol = costum_interp1(dss.lores_tvect, U_opt_, ...
        dss.hires_tvect, 'previous');
    elseif strcmp(dss.input_type, 'foh')
        dss.hires_sol = costum_interp1(dss.lores_tvect, U_opt_, ...
        dss.hires_tvect, 'linear');
    end

    %  lores_sol contains solution sampled with time interval T_ocp
    dss.lores_sol = U_opt_;
end

% ==============================================================================
function dss = prepare(dss)

    dss.tf          = (dss.n_horizon-1)*dss.T_ocp; % Final time
    dss.lores_tvect = 0 : dss.T_ocp : dss.tf;   % Low-resolution time vector
    dss.hires_tvect = 0 : dss.T_dyn : dss.tf;   % High-resolution time vector

end

% ==============================================================================
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
    if (~isequal(size(dss.lb),[dss.n_inputs, dss.n_horizon])|| ...
        ~isequal(size(dss.ub),[dss.n_inputs, dss.n_horizon])|| ...
        ~isequal(size(dss.intial_guesses),[dss.n_inputs, dss.n_horizon]))

        fprintf(['These fields must be: [n_inputs x n_horizon] vectors:' ...
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
