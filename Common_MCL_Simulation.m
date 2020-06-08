%% Simulation of Common MCL model
%{
--------------------------- Description -----------------------------------
Simulation of the common MCL model. Both nonlinear and linearized simulations
are performed. This common model incorporates all the elements, such as
series and parallel springs, rolling diagphragm compliance, that are found
in both the UT and WCS mock loops. The objective of preparing such a model
is to observe the effect of each parameter on the coupling. In the future,
I might try to design a MIMO controller on this model to see how each
element effects the coupling and MIMO control performance.

About the Model : 
    - States : x = [vlv vao xslv xsao vplv vpao xplv xpao xdlv xdao Vlv Vao]^T
    - Inputs : u = [Flv Fao Qrc]^T
    - Outputs : y = [Plv Pao xplv-xpao]^T
    - Disturbance : Qvad

------------------------------ Versions -----------------------------------
v1 : Suraj R Pawar, 6-6-2020
    - Initialize
v2 : Suraj R Pawar, 6-6-2020
    - Added friction at piston and diaphragm interface
    - Removed Qrc calculation from function file to simulation file
%}

close all; clear; clc;
include_us;

%% User Inputs (Only change values here)

    % ---------------------------- LV -------------------------------------    
    mlv = 0.5;  
    mplv = 0.7;    
    
    % --------------------------- AO --------------------------------------
    mao = 0.5;
    mpao = 0.7;    
    
    % -------------------------- Common -----------------------------------
    % Piston assembly. Include mass of fluid in piston mass    
    dia_piston_inch = 2;    % Piston diameter (inches)    
                            
    % Rolling diaphragm
    Cd = 0.001*31959.14283;   
    Ploss_mmHg = 50;          % Loss at piston diaphragm interface
    Qloss_LPM = 10;
    
    % Nonlinear orifice
    cdo = 0.6;                % orifice flow coefficient
    do_cm = 5;                % diameter of orifice (cm)
    
    % Series springs
    k = 0.001*9036.55;        % Spring constant (N/m)
    
    % Air pocket
    V0_mL = 10;               % Initial volume of air in vertical tubes (mL)

    % ----------------- Recirculation Pump --------------------------------
    %{
    % Q : Flow (m^3/s)
    % H : Pressure head (Pa)
    % u : Input voltage to Marco Pump (V)
    % Equation : Q = k2*u - k3*H
    %}
    
    k2 = 1*1.6824e-5;     % (m^3/s/V)
    k3 = 1*2.1456e-10;    % (m^3/s/Pa)

    % ---------------- Simulation Parameters ------------------------------
    t0 = 0;
    dt = 0.0001; 
    tf = 2;

    % ------------------ Actuation commands -------------------------------
    krc = 250;     % Proportional gain for RC
    set_rc = 0;     % Setpoint. (xlv - xa0)
        
    amp_lv = 0.05;    % (N)
    offset_lv = 0;   % (N)
    freq_lv_Hz = 1;   % (Hz)
        
    amp_ao = 0.05;     % (N)
    offset_ao = 0;   % (N)
    freq_ao_Hz = 1;   % (Hz)
    
    % ------------------- Initial Conditions ------------------------------
    %{
    States : x = [vlv vao xslv xsao vplv vpao xplv xpao xdlv xdao Vlv Vao]^T
    1. vlv   : Velocity of LV piston connector (m/s)
    2. vao   : Velocity of AO piston connector (m/s)
    3. xslv  : Compression in LV spring between piston connector and piston
               (m)
    4. xsao  : Compression in AO spring between piston connector and piston
               (m)
    5. vplv  : Velocity of LV piston (m/s)
    6. vpao  : Velocity of AO piston (m/s)
    7. xplv  : Position of LV piston (m)
    8. xpao  : Position of AO piston (m)
    9. xdlv  : Expansion of the LV rolling diaphragm (m)
    10. xdao : Expansion of the LV rolling diaphragm (m)
    11. Vlv  : Reduction in amount of air in LV tube air pocket (m^3)
    12. Vao  : Reduction in amount of air in AO tube air pocket (m^3)
    %}
    
    vlv0 = 0;    
    vao0 = 0;
    xslv0 = 0;
    xsao0 = 0;
    vplv0 = 0;
    vpao0 = 0;
    xplv0 = 0;
    xpao0 = 0;
    xdlv0 = 0;
    xdao0 = 0;
    Vlv0 = 0;
    Vao0 = 0;
        
    x0 = [vlv0; vao0; xslv0; xsao0; vplv0; vpao0; xplv0; xpao0; xdlv0; xdao0; Vlv0; Vao0];

%% Process parameters - conversions etc.
    conversions_all;    

%% Prepare for simulation

    t = [t0 : dt : tf];
    steps = length(t);
    num_states = numel(x0);

    % Nonlinear model
    x = zeros(num_states, steps); x(:,1) = x0;
    Qvad = zeros(1,steps);
    Qrc = zeros(1,steps);
    Flv = zeros(1,steps);   % Control input given to LV VC
    Fao = zeros(1,steps);   % Control input given to AO VC
    urc = zeros(1,steps);   % Control input given to RC Pump
    Plv = zeros(1,steps);  
    Pao = zeros(1,steps);  

    % Linearized model    
    xlin = zeros(num_states, steps); xlin(:,1) = x0;
    Qvadlin = zeros(1,steps);
    Qrclin = zeros(1,steps);
    Flvlin = zeros(1,steps);   % Control input given to LV VC
    Faolin = zeros(1,steps);   % Control input given to AO VC
    urclin = zeros(1,steps);   % Control input given to RC Pump
    Plvlin = zeros(1,steps);  
    Paolin = zeros(1,steps);   
    
    parameters = [mlv; k; mplv; Ap; Cd; Rlvo; Raoo; mpao; mao; V0; Patm; k2; k3; Rplv; Rpao];   
    
%% ODE Simulation
    
    % Get linearized system matrices
    [A, B, C, D, G, H] = func_linmats_all(parameters, x0);
    mats.A = A;
    mats.B = B;
    mats.C = C;
    mats.D = D;
    mats.G = G;
    mats.H = H;
    
    fprintf('Beginning ODE Simulation \n'); tic;
    for i = 2:steps
        % Get information from previous time step
        xp = x(:,i-1);          % Previous state
        xplin = xlin(:,i-1);    % Previous state
        tp = t(i-1);            % Previous time
        tc = t(i);              % Current time                
        Plvp = Plv(i-1);
        Paop = Pao(i-1);
        Plvlinp = Plvlin(i-1);
        Paolinp = Paolin(i-1);

        % Calculate inputs        
        xplv_previous = xp(7);                   
        xpao_previous = xp(8);
        xplvlin_previous = xplin(7);                   
        xpaolin_previous = xplin(8);
        erc = (xplv_previous - xpao_previous) - set_rc;
        erclin = (xplvlin_previous - xpaolin_previous) - set_rc;
        urc(i-1) = krc*erc;               % RC
        urclin(i-1) = krc*erclin;
        Qrc(i-1) = 1*(k2*urc(i-1) - k3*(Plvp - Paop));
        Qrclin(i-1) = 1*(k2*urclin(i-1) - k3*(Plvlinp - Paolinp));
        
        Flv(i-1) = amp_lv * sin(freq_lv * tp) + offset_lv; % LV
        Flvlin(i-1) = Flv(i-1);
        
        Fao(i-1) = amp_ao * sin(freq_ao * tp) + offset_ao; % AO
        Faolin(i-1) = Fao(i-1);
        
        inputs = [Flv(i-1); Fao(i-1); Qrc(i-1)];
        inputslin = [Flvlin(i-1); Faolin(i-1); Qrclin(i-1)];
        
        % ODE solver and current state
        [tout, xout] = ode23t(@(t,x)model_common_mcl(tp,xp,parameters,inputs), [tp tc], xp);    % Nonlinear model         
        [toutlin, xoutlin] = ode45(@(t,x)model_common_mcl_linearized(tp,xplin,parameters,inputslin,mats), [tp tc], xplin);   % Linearized model
        xc = xout(end,:).'; % Current state
        xclin = xoutlin(end,:).';
        x(:,i) = xc;
        xlin(:,i) = xclin;

        % Get current outputs
        [xdot, y] = model_common_mcl(tc,xc,parameters,inputs);
        [xdotlin, ylin, sys] = model_common_mcl_linearized(tc,xclin,parameters,inputslin,mats);

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
    lin_decimation = 1000;    % Decimation to display points from lineraized system
    tlin_decimated = t(1:lin_decimation:end);
    xlin_decimated = xlin(:,[1:lin_decimation:end]);
    Plvlin_decimated = Plvlin(1:lin_decimation:end);
    Paolin_decimated = Paolin(1:lin_decimation:end);
    urclin_decimated = urclin(1:lin_decimation:end);
    Flvlin_decimated = Flvlin(1:lin_decimation:end);
    Faolin_decimated = Faolin(1:lin_decimation:end);
    Qrclin_decimated = Qrclin(1:lin_decimation:end);
    Qvadlin_decimated = Qvadlin(1:lin_decimation:end);
       
    %---------------------- States & Outputs ------------------------------
    figure(1);
    subplot(5,2,1); % FLV
    hold on;
    plot(t,Flv);    
    plot(tlin_decimated,Flvlin_decimated,'ok','MarkerSize',1.5);    
    hold off;
    legend('Nonlinear','Linear');
    title('FLV (N)');

    subplot(5,2,3); % LV Piston position
    hold on;
    plot(t,x(7,:)*100);    
    plot(tlin_decimated,xlin_decimated(7,:)*100,'ok','MarkerSize',1.5);    
    hold off;
    legend('Nonlinear','Linear');
    title('LV Position (cm)');
    
    subplot(5,2,5); % LV Volume
    hold on;
    plot(t,x(11,:)/mL_to_m3);    
    plot(tlin_decimated,xlin_decimated(11,:)/mL_to_m3,'ok','MarkerSize',1.5);
    hold off;
    title('LV Vol (mL)');
    
    subplot(5,2,7); % PLV
    hold on;
    plot(t,Plv/mmHg_to_Pa);
    plot(tlin_decimated,Plvlin_decimated/mmHg_to_Pa,'ok','MarkerSize',1.5);    
    hold off;
    legend('Nonlinear','Linear');    
    title('PLV (mmHg)');
    
    subplot(5,2,2); % FAO
    hold on;
    plot(t,Fao);    
    plot(tlin_decimated,Faolin_decimated,'ok','MarkerSize',1.5);    
    hold off;
    legend('Nonlinear','Linear');
    title('FAO (N)');

    subplot(5,2,4); % AO Piston position
    hold on;
    plot(t,x(8,:)*100);
    plot(tlin_decimated,xlin_decimated(8,:)*100,'ok','MarkerSize',1.5);
    hold off;
    legend('Nonlinear','Linear');
    title('AO Position (cm)');
    
    subplot(5,2,6); % AO Volume
    hold on;
    plot(t,x(12,:)/mL_to_m3);    
    plot(tlin_decimated,xlin_decimated(12,:)/mL_to_m3, 'ok','MarkerSize',1.5);    
    hold off;
    legend('Nonlinear','Linear');
    title('AO Vol (mL)');
    
    subplot(5,2,8); % PAO
    hold on;
    plot(t,Pao/mmHg_to_Pa);
    plot(tlin_decimated,Paolin_decimated/mmHg_to_Pa,'ok','MarkerSize',1.5);    
    hold off;
    legend('Nonlinear','Linear');    
    title('PAO (mmHg)');      
    
    subplot(5,2,9); % LV diaphragm expansion
    hold on;
    plot(t,x(9,:)*100);    
    plot(tlin_decimated,xlin_decimated(9,:)*100,'ok','MarkerSize',1.5);    
    hold off;
    legend('Nonlinear','Linear');
    title('Diaphragm expansion (cm)');
    
    subplot(5,2,10);    % AO Diaphragm expansion
    hold on;
    plot(t,x(10,:)*100);    
    plot(tlin_decimated,xlin_decimated(10,:)*100,'ok','MarkerSize',1.5);
    hold off;
    legend('Nonlinear','Linear');
    title('Diaphragm expansion (cm)');
    
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
    plot(t,(x(7,:) - x(8,:))*100);
    plot(tlin_decimated,(xlin_decimated(7,:) - xlin_decimated(8,:))*100,'ok','MarkerSize',1.5);
    hold off;
    legend('Nonlinear','Linear');
    title('Diff in positions (cm)');        
    
    subplot(4,1,3); % Individual piston positions
    hold on;
    plot(t,x(7,:)*100, t, x(8,:)*100);
    plot(tlin_decimated,xlin_decimated(7,:)*100,'ob',tlin_decimated, xlin_decimated(8,:)*100,'ok','MarkerSize',1.5);
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