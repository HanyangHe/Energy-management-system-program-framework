function [Xnext,Vs,Vinv,S_busInjVec,S_inv_Batt,S_inv_PV,fval,exitflag] = ...
    MicrogridSys_ODEAE_expIterSolver(X,Vs,Vinv,Ref,S_slack_ini,S_load_vec,S_DieselG_vec, ...
                                     V_DC,Xm_g,Zmv_g,Y_Bus,controlParams,K_Vdp,K_fdp, ...
                                     BusOfBattery,BusOfPV,dt,Niter)

    % ---------- persistent caches for speed ----------
    % persistent Jpat Nbus_cached
    Nbus = size(Y_Bus,1);

    % Build / refresh Jacobian sparsity pattern when network size changes
    % if isempty(Jpat) || isempty(Nbus_cached) || Nbus_cached ~= Nbus
    %     Jpat = buildJacPattern(Y_Bus);  % see local function below
    %     Nbus_cached = Nbus;
    % end

    % Use sparse for network matrices
    Y_Bus = sparse(Y_Bus);

    BatteryNum = numel(BusOfBattery);
    SolarNum   = numel(BusOfPV);

    % Initialize iterates
    Vs_iter   = Vs;
    Vinv_iter = Vinv;
    X_iter    = X;

    % Constant used in inverter voltage reconstruction
    M = (sqrt(3)/sqrt(2))*0.5*V_DC;

    for m = 1:Niter
        % ------ Build sparse inverter shunt admittances ------
        % allocate diagonal sparse matrices: number of nonzeros equals #devices
        Yinv_Batt = spalloc(Nbus,Nbus,BatteryNum);
        Yinv_PV   = spalloc(Nbus,Nbus,SolarNum);

        VinvBatt_bus = zeros(Nbus,1);
        VinvPV_bus   = zeros(Nbus,1);

        idx_X   = 0;
        idx_ref = 0;
        idx_Vin = 0;
        idx_Xm  = 0;

        if BatteryNum > 0
            for bb = 1:BatteryNum
                % 1) ODE step for GFM battery
                [Inc1, ~, Edq1] = ode4_singleStep( X_iter(idx_X+1:idx_X+6), dt, ...
                    @PID_Droopcontroller_DynModel, ...
                    Ref(idx_ref+1:idx_ref+4), ...
                    [Vinv_iter(idx_Vin+1); Vs_iter(BusOfBattery(bb))], ...
                    Xm_g(idx_Xm+1), Zmv_g(idx_Xm+1), controlParams(1:10), K_Vdp, K_fdp );

                Xtemp = X_iter(idx_X+1:idx_X+6) + Inc1*dt;
                X_iter(idx_X+1:idx_X+6,1) = Xtemp;

                % 2) Compute new Vinv from E_dq
                th  = Xtemp(2);
                Ex1 = Edq1(1)*cos(th) - Edq1(2)*sin(th);
                Ey1 = Edq1(1)*sin(th) + Edq1(2)*cos(th);
                Vinv_iter(idx_Vin+1) = M*(Ex1 + 1j*Ey1);

                % 3) Update shunt admittance and terminal voltage
                ii = BusOfBattery(bb);
                Yinv_Batt(ii,ii) = Yinv_Batt(ii,ii) + 1/(1j*Xm_g(idx_Xm+1)+Zmv_g(idx_Xm+1));
                VinvBatt_bus(ii) = Vinv_iter(idx_Vin+1);

                idx_X  = idx_X  + 6;
                idx_ref= idx_ref+ 4;
                idx_Xm = idx_Xm + 1;
                idx_Vin= idx_Vin+ 1;
            end
        else
            VinvBatt_bus = Vs_iter;  % keeps logic unchanged
        end

        if SolarNum > 0
            for ss = 1:SolarNum
                % 1) ODE step for PQ PV
                [Inc2, ~, Edq2] = ode4_singleStep( X_iter(idx_X+1:idx_X+6), dt, ...
                    @PID_PQcontroller_DynModel, ...
                    Ref(idx_ref+1:idx_ref+2), ...
                    [Vinv_iter(idx_Vin+1); Vs_iter(BusOfPV(ss))], ...
                    Xm_g(idx_Xm+1), Zmv_g(idx_Xm+1), controlParams(11:20) );

                Xtemp = X_iter(idx_X+1:idx_X+6) + Inc2*dt;
                X_iter(idx_X+1:idx_X+6,1) = Xtemp;

                % 2) Compute new Vinv from E_dq
                th  = Xtemp(2);
                Ex2 = Edq2(1)*cos(th) - Edq2(2)*sin(th);
                Ey2 = Edq2(1)*sin(th) + Edq2(2)*cos(th);
                Vinv_iter(idx_Vin+1) = M*(Ex2 + 1j*Ey2);

                % 3) Update shunt admittance and terminal voltage
                jj = BusOfPV(ss);
                Yinv_PV(jj,jj)   = Yinv_PV(jj,jj) + 1/(1j*Xm_g(idx_Xm+1)+Zmv_g(idx_Xm+1));
                VinvPV_bus(jj)   = Vinv_iter(idx_Vin+1);

                idx_X  = idx_X  + 6;
                idx_ref= idx_ref+ 2;
                idx_Xm = idx_Xm + 1;
                idx_Vin= idx_Vin+ 1;
            end
        else
            VinvPV_bus = Vs_iter;
        end

        % ------ Power-flow solve (FSOLVE) ------
        % Unknown vector: [Pslack; Re(V2..VN); Qslack; Im(V2..VN)]
        Vs_init = Vs_iter(2:end);
        Vinf    = 1+0j;


        z0 = [ real(S_slack_ini);
               real(Vs_init);
               imag(S_slack_ini);
               imag(Vs_init) ];

        % opts = optimoptions('fsolve', ...
        %     'Display','off');

        opts = optimoptions('fsolve', ...
            'Display','off',...
        'Algorithm','levenberg-marquardt');

        % opts = optimoptions('fsolve', ...
        %     'Display','off', ...
        %     'Algorithm','levenberg-marquardt', ...
        %     'JacobPattern', Jpat, ...                % sparsity pattern
        %     'FiniteDifferenceType','forward', ...    % faster FD
        %     'ScaleProblem','jacobian', ...
        %     'FunctionTolerance',1e-9, ...
        %     'StepTolerance',1e-9, ...
        %     'MaxIterations', 100);

        [z_sol, fval, exitflag] = fsolve(@(z) complex_power_flow_eq( ...
                z, Vinf, Y_Bus, Yinv_Batt, Yinv_PV, VinvBatt_bus, VinvPV_bus, ...
                S_load_vec, S_DieselG_vec), z0, opts);

        % Unpack: z = [Pr; Qi]
        n  = numel(z_sol)/2;
        Pr = z_sol(1:n);
        Qi = z_sol(n+1:end);
        VsR = [1; Pr(2:end)];
        VsI = [0; Qi(2:end)];
        Vs_iter = VsR + 1j*VsI;

        % keep warm start for next call
        if exitflag > 0
            z_last = z_sol;
        else
            % if failed, reset warm start to avoid poisoning next step
            z_last = [];
        end
    end

    % ------ outputs ------
    Xnext = X_iter;
    Vs    = Vs_iter;
    Vinv  = Vinv_iter;

    S_busInjVec = Vs .* conj(Y_Bus*Vs);
    S_inv_Batt  = Vs .* conj(Yinv_Batt * (VinvBatt_bus - Vs));
    S_inv_PV    = Vs .* conj(Yinv_PV   * (VinvPV_bus   - Vs));

    % ================= local functions =================
    function F = complex_power_flow_eq(decision_real_imag, Vinf, Ybus, Yinv_Batt,Yinv_PV, ...
                                       VinvBatt_bus, VinvPV_bus, S_load, S_DieselG)
        n    = numel(decision_real_imag)/2;
        Real = decision_real_imag(1:n);
        Imag = decision_real_imag(n+1:end);

        Vs_real = [real(Vinf); Real(2:end)];
        Vs_imag = [imag(Vinf); Imag(2:end)];
        Vbus    = Vs_real + 1i*Vs_imag;

        % Slack complex power
        Psl = Real(1);  Qsl = Imag(1);
        S_sys         = zeros(size(S_load));
        S_sys(1,1)    = Psl + 1j*Qsl;

        % Inverter injections
        S_inv = Vbus .* conj(Yinv_Batt * (VinvBatt_bus - Vbus)) + ...
                Vbus .* conj(Yinv_PV   * (VinvPV_bus   - Vbus));

        left  = Vbus .* conj(Ybus * Vbus);
        right = S_inv + S_sys + S_load + S_DieselG;

        F = [ real(left - right); imag(left - right) ];
    end

    function [IncRate,dotX,Edq] = ode4_singleStep(Xin, h, f_dyn, varargin)
        X1 = Xin;
        [Edq,k1] = f_dyn(varargin{:}, X1);

        X2 = Xin + 0.5*h*k1;
        [~,k2] = f_dyn(varargin{:}, X2);

        X3 = Xin + 0.5*h*k2;
        [~,k3] = f_dyn(varargin{:}, X3);

        X4 = Xin + h*k3;
        [~,k4] = f_dyn(varargin{:}, X4);

        IncRate = (k1 + 2*k2 + 2*k3 + k4)/6;
        dotX    = k1;
    end

    % function Jp = buildJacPattern(Ybus_full)
    %     % Jacobian sparsity for PF residual wrt variables z = [Pr;Qi]
    %     % We treat unknowns as [Pslack; Re(V2..VN); Qslack; Im(V2..VN)]
    %     N  = size(Ybus_full,1);
    %     nb = N-1;
    %     A  = spones((Ybus_full ~= 0) | (Ybus_full.' ~= 0)); % adjacency including self
    % 
    %     % Rows: N real mismatches + N imag mismatches
    %     % Cols: 1 + nb + 1 + nb = 2N unknowns
    %     Jp = spalloc(2*N, 2*N, 10*N);  % rough nnz guess
    % 
    %     % Real mismatches depend on Pslack (1), Re(V2..VN), Im(V2..VN)
    %     Jp(1:N, 1)              = 1;          % wrt Pslack
    %     Jp(1:N, 2:(1+nb))       = A(:,2:end); % wrt Re(V2..VN)
    %     Jp(1:N, (3+nb):end)     = A(:,2:end); % wrt Im(V2..VN)
    % 
    %     % Imag mismatches depend on Qslack, Re/Im of neighbors
    %     Jp((N+1):end, (2+nb))   = 1;          % wrt Qslack
    %     Jp((N+1):end, 2:(1+nb)) = A(:,2:end);
    %     Jp((N+1):end, (3+nb):end) = A(:,2:end);
    % end
end
