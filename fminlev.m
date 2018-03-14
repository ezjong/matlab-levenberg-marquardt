function x = fminlev(fenergy,flinear,fshowiter,tolerance,iterations,x0)
    %----------------------------------------------------------------------------
    % Levenberg-Marquardt algorithm with inexact line-search.
    %
    % Assume the nonlinear energy to be
    %
    %   EQ(x) = ||a(x) + b||^2 + g( a(x) )
    %
    % where a(x) is a nonlinear vector of x and g is a scalar function of a(x).
    % We can then linearize the vector a(x) around an operating point xb as
    %    a(d) ~ a(xb) + SIGMA_i da_di(xb) di
    % where d = (d_1 d_2 ...) is a step from the current operating point xb, i.e.
    % d := x - xb.
    %
    % Replacing the nonlinear a(x) by the linearized vector a(d) into the energy yields the
    % linearized energy EQLIN(xb,d), which we rewrite as
    %
    %   EQLIN(xb,d) = 0.5 d^T H(xb)^T d  + f(xb)^T d + 0.5 c
    %
    % INPUTS:
    %  fenergy:     callback that returns E = fenergy(x), the nonlinear energy value at x.
    %
    %  flinear:     callback [H,f,c] = fsystem(xb), that returns the (linear) energy around the
    %               the operating point xb, i.e. EQLIN(xb,d), implied by H(xb), f(xb) and c(xb).
    %
    %
    %  fshowiter    callback for displaying infos about the current iteration
    %               aka fshow_iter(x,mu,energy,progress,iteration)
    %
    %  tolerance    minimum progress tolerance
    %  iterations   maximum iterations
    %  x0           initial guess
    %  lb           lower bounds
    %  ub           upper bounds
    %
    %----------------------------------------------------------------------------


    %-------------------------
    % Some Constants
    %-------------------------
    MAX_MU        = 1e12;  % minimum trust region
    BETA0         = 0.3;   % wolfe condition c1
    BETA1         = 0.9;   % wolfe condition c2
    %MU0           = 1;    % initial (inverse) trust region radius
    MU0           = 1e-2;  % initial (inverse) trust region radius
    V             = 2;     % factor for adjusting the step length



    %------------------------------------------------------------------------
    %   MIN_MU ~ maximum trust region:
    %
    %   Influences the size of steps taken by the optimizer
    %   For highly nonlinear functions we do not want to take the optimizer too large
    %   steps. Therefore, we decrease the maximum trust region.
    %------------------------------------------------------------------------
    %MIN_MU = 0.25;
    MIN_MU = 1e-5;
    %iMIN_MU = 0.25;
    %MIN_MU = 0.25;

    %-------------------------
    % Current estimate
    %-------------------------
    x = x0;

    %-------------------------
    % Initial energy
    %-------------------------
    E0 = fenergy(x);

    %----------------------------------------
    % Current (inverse) trust region radius
    %----------------------------------------
    mu = MU0;

    %----------------------------------------
    % Display zero iteration
    %----------------------------------------
    if ~isempty(fshowiter)
        fshowiter(x,E0,mu,0,0);
    end

    %------------------------------------------------------
    % For k = 1, 2, ...
    %------------------------------------------------------
    for k = 1:iterations

        %-------------------------------------
        % (1) Compute r(p_k) and J(p_k)
        %-------------------------------------
        [H,f,c] = flinear(x);

        % do until d is acceptable
        while true

            %-----------------------------------------------------------------
            % Search direction
            %-----------------------------------------------------------------
            d = update_d(H,f,mu); % d step

            progress = norm(d);
            if progress < tolerance
                break;
            end

            %--------------------------------------------------------------
            % Test whether d is acceptable by
            %      e_mu =   ( E(x) - E(x + d)        )
            %             / ( E(x) - E_linear(x + d) )
            %--------------------------------------------------------------
            E1 = fenergy(x + d);
            E1_lin = 0.5*(d'*H*d) + f'*d + 0.5*c;
            eps_mu = (E0 - E1) / (E0 - E1_lin);

            %---------------------------------------------------
            % if e_mu <= beta0: Set mu = 2 mu and recompute d
            % if e_mu >= beta1: Set mu = mu/2 and keep d
            %--------------------------------------------------
            if eps_mu > BETA0
                break;
            end

            % e_mu <= beta0...

            % try to decrease trust region
            mu = min(V*mu,MAX_MU);
        end

        % e_mu > beta0...

        %-------------------------------------------
        % Finished if ||d|| <= tolerance
        %-------------------------------------------
        if progress < tolerance
            break;
        end

        if eps_mu >= BETA1
            %  try to increase trust region
           mu = max(mu/V,MIN_MU);
        end

        %-----------------------------------
        % Set x_k+1 = x_k + d
        %-----------------------------------
        x = x + d;

        %-----------------------------------
        % Remember energy
        %-----------------------------------
        E0 = E1;

        %-----------------------------------
        % Show iteration per callback
        %-----------------------------------
        if ~isempty(fshowiter)
            fshowiter(x,E1,mu,progress,k);
        end

    end

end

%----------------------------
% Update for step d = x-x0
%----------------------------
function d = update_d(H,f,mu)
    ndim = length(f);

    % closed-form solution to quadratic form
    A = 0.5*(H + H') + mu^2*eye(ndim);
    b = - f;

    % solve in scaled space (plus regularization)
    d = A \ b;
end

