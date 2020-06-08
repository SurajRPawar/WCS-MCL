%% RGA Analysis on UT hMCL design
%{
-------------------------- Description ------------------------------------
RGA Analysis on common MCL Linearized model. To learn more about RGA analysis,
refer to Ch 3, section 3.6 and Ch 10, section 10.8 from Skogestad and
Postlethwaite.

- States : x = [vlv vao xslv xsao vplv vpao xplv xpao xdlv xdao Vlv Vao]^T
- Inputs : u = [Flv Fao Qrc]^T
- Outputs : y = [Plv Pao xplv-xpao]^T
- Disturbance : Qvad

-------------------------- Versions ---------------------------------------
v1 : Suraj R Pawar, 6-7-2020
    - Initialize
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
    do_cm = 1;                % diameter of orifice (cm)
    
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
    krc = 25;     % Proportional gain for RC
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
              
    parameters = [mlv; k; mplv; Ap; Cd; Rlvo; Raoo; mpao; mao; V0; Patm; k2; k3; Rplv; Rpao];   
                 
    
%% Form linearized system
    [A, B, C, D, ~, ~] = func_linmats_all(parameters, x0);            
    syslin = ss(A,B,C,D);
    
%% RGA Matrix    
    num_meas = size(C,1);
    num_ins = size(B,2);

    w = logspace(1,2);  % Frequency in rad/s
    steps = numel(w);        
    RGA_mat = zeros(num_meas, num_ins, steps);
    RGA_mag = zeros(num_meas, num_ins, steps);
    RGA_no = zeros(1,steps);
    
    for i = 1:steps
        G = freqresp(syslin,w(i));
        %RGA_mat(:,:,i) = ucrga(G);
        RGA_mat(:,:,i) = G.*inv(G).'; % RGA matrix
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
   hold off; grid on;
   ylabel('Mag (dB)');
   title('PLV');
   legend('Flv', 'Fao','urc');
   ylim([-50 20]);
   set(gca,'XScale','log');
      
   subplot(3,1,2);
   hold on;
   plot(w,squeeze(RGA_dB(2,1,:)));
   plot(w,squeeze(RGA_dB(2,2,:)),'rx-','MarkerSize',2);
   plot(w,squeeze(RGA_dB(2,3,:)),'ko-','MarkerSize',2);
   hold off; grid on;
   ylabel('Mag (dB)');
   title('PAO');
   legend('Flv', 'Fao','urc');
   ylim([-50 20]);
   set(gca,'XScale','log');
   
   subplot(3,1,3);
   hold on;
   plot(w,squeeze(RGA_dB(3,1,:)));
   plot(w,squeeze(RGA_dB(3,2,:)),'rx-','MarkerSize',2);
   plot(w,squeeze(RGA_dB(3,3,:)),'ko-','MarkerSize',2);
   hold off; grid on;
   ylabel('Mag (dB)');
   title('xplv - xpao');
   xlabel('Frequency (Hz)');
   legend('Flv', 'Fao','urc');
   ylim([-50 20]);
   set(gca,'XScale','log');
   
   figure;
   plot(w,RGA_no,'.-'); grid on;
   xlabel('Frequency (Hz)');
   ylabel('||RGA - I|| sum');
   title('RGA Number');
   set(gca, 'XScale','log');   