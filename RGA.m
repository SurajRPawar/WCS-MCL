%% RGA Analysis on WCS hMCL design
%{
-------------------------- Description ------------------------------------
RGA Analysis on WCS MCL Linearized model. To learn more about RGA analysis,
refer to Ch 3, section 3.6 and Ch 10, section 10.8 from Skogestad and
Postlethwaite.

States : x0 = [ilv0; vlv0; xlv0; iao0; vao0; xao0; Vlv0; Vao0]
Outputs : y = [Plv; Pao; xlv - xao]
Inputs : u = [urc; ulv; uao]

-------------------------- Versions ---------------------------------------
v1 : Suraj R Pawar, 5-20-2020
    - Initialize
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
    Cd = 1*31959.14283;       % Comliance of the diaphragm [N/m]
    
    % At Qloss m^3/s, we should get Ploss pressure loss due to friction
    Ploss_piston_mmHg = 1*50; % Frictional loss in pressure at piston assembly (mmHg)
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
    krc = 2500; % Proportional gain for RC
    set_rc = 0; % Setpoint. (xlv - xa0)
    
    klv = 0.0008;
    set_lv_mmHg = 50;    % (mmHg)
    
    kao = 0.0009;
    set_ao_mmHg = 5; % (mmHg)
    
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
    
    ilv = 0;
    vlv = 0;
    xlv = 0;
    iao = 0;
    vao = 0;
    xao = 0;
    Vlv = 0;
    Vao = 0;
    
%% Conversions
    conversions;
    
    parameters = [Lvclv; Lvcao; Rvclv; Rvcao; rvc; mplv; mpao; Rplv; Rpao; Ap;...
                  k2; k3; Patm; V0; Cd];   
              
    IC = [ilv, vlv, xlv, iao, vao, xao, Vlv, Vao];
    
%% Form linearized system
    [A, B, C, D, ~, ~] = func_linmats(parameters, IC);            
    syslin = ss(A,B,C,D);
    
%% RGA Matrix    
    num_meas = size(C,1);
    num_ins = size(B,2);

    w = logspace(-1,2);  % Frequency in rad/s
    steps = numel(w);        
    RGA_mat = zeros(num_meas, num_ins, steps);
    RGA_mag = zeros(num_meas, num_ins, steps);
    RGA_no = zeros(1,steps);
    
    for i = 1:steps
        G = freqresp(syslin,w(i));
        RGA_mat(:,:,i) = G.*pinv(G).'; % RGA matrix
        RGA_mag(:,:,i) = abs(RGA_mat(:,:,i));
        RGA_no(i) = sum(sum(abs(RGA_mat(:,:,i) - eye(num_meas))));
    end  
   RGA_dB = 20*log10(RGA_mag); 
   
   figure;
   set(groot,'defaultLineLineWidth',0.8);
   
   subplot(3,1,1);
   hold on;
   plot(w,squeeze(RGA_dB(1,1,:)));
   plot(w,squeeze(RGA_dB(1,2,:)),'rx-','MarkerSize',2);
   plot(w,squeeze(RGA_dB(1,3,:)),'ko-','MarkerSize',2);
   hold off;
   ylabel('Mag (dB)');
   title('PLV');
   legend('ulv','uao', 'urc');
   ylim([-120 20]);
      
   subplot(3,1,2);
   hold on;
   plot(w,squeeze(RGA_dB(2,1,:)));
   plot(w,squeeze(RGA_dB(2,2,:)),'rx-','MarkerSize',2);
   plot(w,squeeze(RGA_dB(2,3,:)),'ko-','MarkerSize',2);
   hold off;
   ylabel('Mag (dB)');
   title('PAO');
   legend('ulv','uao','urc');
   ylim([-120 20]);
   
   subplot(3,1,3);
   hold on;
   plot(w,squeeze(RGA_dB(3,1,:)));
   plot(w,squeeze(RGA_dB(3,2,:)),'rx-','MarkerSize',2);
   plot(w,squeeze(RGA_dB(3,3,:)),'ko-','MarkerSize',2);
   hold off;
   ylabel('Mag (dB)');
   title('dx = xlv - xao');
   xlabel('Frequency (Hz)');
   legend('ulv','uao','urc');
%    ylim([-50 20]);
