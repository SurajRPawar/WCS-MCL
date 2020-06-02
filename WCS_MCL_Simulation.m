%% Simulation of WCS MCL
%{
--------------------------- Description -----------------------------------
Simulation of the WCS MCL model. Both nonlinear and linearized simulations
are performed. 
About the Model : 
    - States : x = [ilv vlv xlv iao vao xao Vlv Vao]^T
    - Inputs : u = [urc ulv uao]^T
    - Outputs : y = [Plv Pao xlv-xao]^T
    - Disturbance : Qvad

Assumptions : 
- Assume that the recirculation pump follows steady state dynamics
  because we don't have information of the inductances and internal
  resistances of the motor and rotor. 
- Assume that the diaphragm compliance is linear

------------------------------ Versions -----------------------------------
v1 : Suraj R Pawar, 5-15-2020
    - Initialize
v2 : Suraj R Pawar, 6-1-2020
    - Add compliance of diaphragm
%}

close all; clear; clc;
include_us;

%% User Inputs (Only change values here)

    % ---------------------------- LV -------------------------------------
    % Voice Coil Actuator
    Lvclv = 2.4e-3;         % Voice coil Inductance (H)
    Rvclv = 3;              % Voice coil Resistance (Ohms)
    rvc = 14.5;             % Voice coil gyration constant (N/A)

    % Piston assembly. Include mass of fluid in piston mass
    m_piston_head = 0.7;    % (kg)
    dia_piston_inch = 3;    % Piston diameter (inches)
    fluid_section_cm = 5;   % length of blood plasma cylindrical...
                            % portion, area is same as piston head.  
    Cd = 0.001*31959.14283;       % Comliance of the diaphragm [N/m]
    
    % At Qloss m^3/s, we should get Ploss pressure loss due to friction
    Ploss_piston_mmHg = 50; % Frictional loss in pressure at piston assembly (mmHg)
    Qloss_piston_LPM = 10;  % (LPM)

    % Air pocket
    V0_mL = 10;             % Initial volume of air in vertical tubes (mL)

    % ----------------- Recirculation Pump --------------------------------
    %{
    % Q : Flow (m^3/s)
    % H : Pressure head (Pa)
    % u : Input voltage to Marco Pump (V)
    % Equation : Q = k2*u - k3*H
    %}
    
    k2 = 1.6824e-5;     % (m^3/s/V)
    k3 = 2.1456e-10;    % (m^3/s/Pa)

    % ---------------- Simulation Parameters ------------------------------
    t0 = 0;
    dt = 0.001; 
    tf = 2;

    % ------------------ Actuation commands -------------------------------
    krc = 3000; % Proportional gain for RC
    set_rc = 0; % Setpoint. (xlv - xa0)
    
    klv = 0.0009;
    set_lv_mmHg = 20;    % (mmHg)
    
    kao = 0.0009;
    set_ao_mmHg = 0; % (mmHg)
    
    % ------------------- Initial Conditions ------------------------------
    %{
    States : [ilv; vlv; xlv; iao; vao; xao; Vlv; Vlv]
    ilv : Current in LV VC (A)
    vlv : Velocity of LV Piston (m/s)
    xlv : Position of LV piston (m); Positive dir is away from face plate
    iao : Current in AO VC (A)
    vao : Velocity of AO Piston (m/s)
    xao : Position of LV piston (m)
    Vlv : Volume in LV capacitor (m^3)
    Vao : Volume in AO capacitor (m^3)
    %}
    
    ilv0 = 0;
    vlv0 = 0;
    xlv0 = 0;
    iao0 = 0;
    vao0 = 0;
    xao0 = 0;
    Vlv0 = 0;
    Vao0 = 0;
    
    x0 = [ilv0; vlv0; xlv0; iao0; vao0; xao0; Vlv0; Vao0];

%% Process parameters - conversions etc.
    conversions;    

%% Prepare for simulation

    t = [t0 : dt : tf];
    steps = length(t);
    num_states = numel(x0);

    % Nonlinear model
    x = zeros(num_states, steps); x(:,1) = x0;
    Qvad = zeros(1,steps);
    Qrc = zeros(1,steps);
    ulv = zeros(1,steps);   % Control input given to LV VC
    uao = zeros(1,steps);   % Control input given to AO VC
    urc = zeros(1,steps);   % Control input given to RC Pump
    Plv = zeros(1,steps);   Plv(1) = Patm*(V0/(V0 - Vlv0) - 1);
    Pao = zeros(1,steps);   Pao(1) = Patm*(V0/(V0 - Vao0) - 1);

    % Linearized model
    xlin = zeros(num_states, steps); xlin(:,1) = x0;
    Qvadlin = zeros(1,steps);
    Qrclin = zeros(1,steps);
    ulvlin = zeros(1,steps);   % Control input given to LV VC
    uaolin = zeros(1,steps);   % Control input given to AO VC
    urclin = zeros(1,steps);   % Control input given to RC Pump
    Plvlin = zeros(1,steps);   Plvlin(1) = Patm*(V0/(V0 - Vlv0) - 1);
    Paolin = zeros(1,steps);   Paolin(1) = Patm*(V0/(V0 - Vao0) - 1);
    
    parameters = [Lvclv; Lvcao; Rvclv; Rvcao; rvc; mplv; mpao; Rplv; Rpao; Ap;...
                  k2; k3; Patm; V0; Cd];   
    
%% ODE Simulation

    fprintf('Beginning ODE Simulation \n'); tic;
    for i = 2:steps
        % Get information from previous time step
        xp = x(:,i-1);  % Previous state
        xplin = xlin(:,i-1);  % Previous state
        tp = t(i-1);    % Previous time
        tc = t(i);      % Current time                

        % Calculate inputs        
        xlvp = xp(3);                   
        xaop = xp(6);
        xlvplin = xplin(3);                   
        xaoplin = xplin(6);
        erc = (xlvp - xaop) - set_rc;
        erclin = (xlvplin - xaoplin) - set_rc;
        urc(i-1) = krc*erc;               % RC
        urclin(i-1) = krc*erclin;
        
        elv = set_lv - Plv(i-1);
        elvlin = set_lv - Plvlin(i-1);
        ulv(i-1) = klv*elv;               % LV
        ulvlin(i-1) = klv*elvlin;
        
        eao = set_ao - Pao(i-1);
        eaolin = set_ao - Paolin(i-1);
        uao(i-1) = kao*eao;               % AO
        uaolin(i-1) = kao*eaolin;
        
        inputs = [ulv(i-1); uao(i-1); urc(i-1)];
        inputslin = [ulvlin(i-1); uaolin(i-1); urclin(i-1)];
        
        % ODE solver and current state
        [tout, xout] = ode45(@(t,x)model_wcs_mcl(tp,xp,parameters,inputs), [tp tc], xp);    % Nonlinear model         
        [toutlin, xoutlin] = ode45(@(t,x)model_wcs_mcl_linearized(tp,xplin,parameters,inputslin,x0), [tp tc], xplin);   % Linearized model
        xc = xout(end,:).'; % Current state
        xclin = xoutlin(end,:).';
        x(:,i) = xc;
        xlin(:,i) = xclin;

        % Get current outputs
        [xdot, y] = model_wcs_mcl(tc,xc,parameters,inputs);
        [xdotlin, ylin, sys] = model_wcs_mcl_linearized(tc,xclin,parameters,inputslin,x0);

        % Separate variables
        Qvad(i) = y(1);
        Qrc(i) = y(2);        
        Plv(i) = y(3);
        Pao(i) = y(4);
        Qvadlin(i) = ylin(1);
        Qrclin(i) = ylin(2);        
        Plvlin(i) = ylin(3);
        Paolin(i) = ylin(4);
    end
    fprintf('Finished ODE Simulation in %.2f seconds \n', toc);

%% Analysis of linearized system
    lin_wcs_mcl = ss(sys.A, sys.B, sys.C, sys.D);
    ctrb_mat = ctrb(sys.A, sys.B);
    obsv_mat = obsv(sys.A, sys.C);
    rank_ctrb = rank(ctrb_mat);
    rank_obsv = rank(obsv_mat);
    fprintf('Number of states: %d \n', num_states);
    fprintf('Rank of controllability matrix: %d \n', rank_ctrb);
    fprintf('Rank of observability matrix: %d \n', rank_obsv);
    
%% Figures
    lin_decimation = 60;    % Decimation to display points from lineraized system
    tlin_decimated = t(1:lin_decimation:end);
    xlin_decimated = xlin(:,[1:lin_decimation:end]);
    Plvlin_decimated = Plvlin(1:lin_decimation:end);
    Paolin_decimated = Paolin(1:lin_decimation:end);
    urclin_decimated = urclin(1:lin_decimation:end);
    ulvlin_decimated = ulvlin(1:lin_decimation:end);
    uaolin_decimated = uaolin(1:lin_decimation:end);
    Qrclin_decimated = Qrclin(1:lin_decimation:end);
    Qvadlin_decimated = Qvadlin(1:lin_decimation:end);
    
    %---------------------- States & Outputs ------------------------------
    figure(1);
    subplot(5,2,1); % VC Current
    hold on;
    plot(t,x(1,:));    
    plot(tlin_decimated,xlin_decimated(1,:),'ok','MarkerSize',1.5);    
    hold off;
    legend('Nonlinear','Linear');
    title('LV VC Current (A)');

    subplot(5,2,3); % LV Piston position
    hold on;
    plot(t,x(3,:)*100);    
    plot(tlin_decimated,xlin_decimated(3,:)*100,'ok','MarkerSize',1.5);    
    hold off;
    legend('Nonlinear','Linear');
    title('LV Position (cm)');
    
    subplot(5,2,5); % LV Volume
    hold on;
    plot(t,x(7,:)/mL_to_m3);    
    plot(tlin_decimated,xlin_decimated(7,:)/mL_to_m3,'ok','MarkerSize',1.5);
    hold off;
    title('LV Vol (mL)');
    
    subplot(5,2,7); % PLV
    hold on;
    plot(t,Plv/mmHg_to_Pa);
    plot(tlin_decimated,Plvlin_decimated/mmHg_to_Pa,'ok','MarkerSize',1.5);
    plot(t,(set_lv/mmHg_to_Pa)*ones(1,steps),':k');
    hold off;
    legend('Nonlinear','Linear','Setpoint');    
    title('PLV (mmHg)');
    
    subplot(5,2,2); % AO VC Current
    hold on;
    plot(t,x(4,:));    
    plot(tlin_decimated,xlin_decimated(4,:),'ok','MarkerSize',1.5);    
    hold off;
    legend('Nonlinear','Linear');
    title('AO Current (A)');

    subplot(5,2,4); % AO Piston position
    hold on;
    plot(t,x(6,:)*100);
    plot(tlin_decimated,xlin_decimated(6,:)*100,'ok','MarkerSize',1.5);
    hold off;
    legend('Nonlinear','Linear');
    title('AO Position (cm)');
    
    subplot(5,2,6); % AO Volume
    hold on;
    plot(t,x(8,:)/mL_to_m3);    
    plot(tlin_decimated,xlin_decimated(8,:)/mL_to_m3, 'ok','MarkerSize',1.5);    
    hold off;
    legend('Nonlinear','Linear');
    title('AO Vol (mL)');
    
    subplot(5,2,8); % PAO
    hold on;
    plot(t,Pao/mmHg_to_Pa);
    plot(tlin_decimated,Paolin_decimated/mmHg_to_Pa,'ok','MarkerSize',1.5);
    plot(t,(set_ao/mmHg_to_Pa)*ones(1,steps),':k');
    hold off;
    legend('Nonlinear','Linear','Setpoint');    
    title('PAO (mmHg)');      
    
    subplot(5,2,9); % LV VC Current
    hold on;
    plot(t,ulv);    
    plot(tlin_decimated,ulvlin_decimated,'ok','MarkerSize',1.5);    
    hold off;
    legend('Nonlinear','Linear');
    title('U LV (V)');
    
    subplot(5,2,10);    % AO VC Current
    hold on;
    plot(t,uao);    
    plot(tlin_decimated,uaolin_decimated,'ok','MarkerSize',1.5);
    hold off;
    legend('Nonlinear','Linear');
    title('U AO (V)');
    
    %----------------------- Flows and RC Control -------------------------
    figure;
    subplot(4,1,1); % Flows
    hold on;
    plot(t,Qrc/LPM_to_m3s,t,Qvad/LPM_to_m3s);
    plot(tlin_decimated,Qrclin_decimated/LPM_to_m3s,'ok','MarkerSize',1.5);
    hold off;
    legend('RC Nonlinear','VAD','RC Linear');
    title('Flows (L/min)');    
    
    subplot(4,1,2); % xlv - xao
    hold on;
    plot(t,(x(3,:) - x(6,:))*100);
    plot(tlin_decimated,(xlin_decimated(3,:) - xlin_decimated(6,:))*100,'ok','MarkerSize',1.5);
    hold off;
    legend('Nonlinear','Linear');
    title('Diff in positions (cm)');        
    
    subplot(4,1,3); % Individual piston positions
    hold on;
    plot(t,x(3,:)*100, t, x(6,:)*100);
    plot(tlin_decimated,xlin_decimated(3,:)*100,'ob',tlin_decimated, xlin_decimated(6,:)*100,'ok','MarkerSize',1.5);
    hold off;
    title('Positions (cm)');    
    legend('xlv Nonlin','xao Nonlin','xlv Lin','xao Lin');
    
    subplot(4,1,4); % RC Control voltage
    hold on;
    plot(t,urc);
    plot(tlin_decimated,urclin_decimated,'ok','MarkerSize',1.5);
    hold off;
    legend('Nonlinear','Linear');
    title('RC Control Voltage (V)');    