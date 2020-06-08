%% RGA Analysis on UT hMCL design
%{
-------------------------- Description ------------------------------------
RGA Analysis on UT MCL Linearized model. To learn more about RGA analysis,
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

    % Aortic Side
    rao   = 9.6;         % voice coil gyration constant [N/A]
    Ap    = 0.00456;     % Area of the piston [m^2]
    Rvcao = 3.9;         % Voice coil resistance [Ohms]
    Lvcao = 1.6e-3;      % Voice coil inductance [H]
    ktao  = 1*9036.55;     % Spring constant [N/m]
    kmao  = 1*9036.55;     % Spring constant [N/m]    
    mao   = 0.15;        % Mass of the piston [kg]
    Rpao  = 250;         % Bearing losses in piston [N.s/m] 
    mpao  = 0.1;         % Mass of the final piston head [kg]
    Cd    = 1*31959.14283; % Comliance of the diaphragm [N/m]
    Cao   = 4.658e-7;    % Tank capacitance [m^4.s^2/kg]
    
    % LV side
    rlv   = 14.5;        % voice coil gyration constant [N/A]    
    Rvclv = 3;           % Voice coil resistance [Ohms]
    Lvclv = 2.4e-3;      % Voice coil inductance [H]
    kmlv  = 1*9036.55;     % Spring constant [N/m]    
    mlvt  = 0.774;       % Mass of the piston [kg]
    Rplv  = 350;         % Bearing losses in piston [N.s/m] 
    mplv  = 0.1;         % Mass of the final piston head [kg]    
    Clv   = 4.658e-7;    % Tank capacitance [m^4.s^2/kg]
    
    % Pumps
    Rrc  = 4.94e8;%1.597e9;   % Recirculation pump resistance [Pa^5/m^3]
    Rvad = 1e6 ;%1.5409e5;    % VAD pump resistance [Pa^5/m^3]
    
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
    
    iao0 = 0;
    ilv0 = 0;
    xtao0 = 0;
    vao0 = 0;
    vlvt0 = 0;
    xlvt0 = 0;
    xmao0 = 0;
    xmlv0 = 0;
    vpao0 = 0;
    xpao0 = 0;
    vplv0 = 0;
    xplv0 = 0;
    xdao0 = 0;
    xdlv0 = 0;
    Vao0 = 0;
    Vlv0 = 0;   % We won't need these initial conditions because the UT mock loop is already linear
    
    IC = [iao0; ilv0; xtao0; vao0; vlvt0; xlvt0; xmao0; xmlv0; vpao0; xpao0;...
          vplv0; xplv0; xdao0; xdlv0; Vao0; Vlv0];
              
    parameters = [Lvcao; Rvcao; rao; ktao; mao; kmao; mpao; Rpao; Cd; Ap; Cao; Lvclv; Rvclv;...
                  rlv; mlvt; kmlv; mplv; Rplv; Clv; Rvad; Rrc];   
                 
    
%% Form linearized system
    [A, B, C, D, ~, ~] = func_linmats_ut_mcl(parameters, IC);            
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
   hold off;
   ylabel('Mag (dB)');
   title('PLV');
   legend('ulv', 'Qrc','uao');
   ylim([-50 20]);
   set(gca,'XScale','log');
      
   subplot(3,1,2);
   hold on;
   plot(w,squeeze(RGA_dB(2,1,:)));
   plot(w,squeeze(RGA_dB(2,2,:)),'rx-','MarkerSize',2);
   plot(w,squeeze(RGA_dB(2,3,:)),'ko-','MarkerSize',2);
   hold off;
   ylabel('Mag (dB)');
   title('PAO');
   legend('ulv', 'Qrc','uao');
   ylim([-50 20]);
   set(gca,'XScale','log');
   
   subplot(3,1,3);
   hold on;
   plot(w,squeeze(RGA_dB(3,1,:)));
   plot(w,squeeze(RGA_dB(3,2,:)),'rx-','MarkerSize',2);
   plot(w,squeeze(RGA_dB(3,3,:)),'ko-','MarkerSize',2);
   hold off;
   ylabel('Mag (dB)');
   title('xplv - xpao');
   xlabel('Frequency (Hz)');
   legend('ulv', 'Qrc','uao');
   ylim([-50 20]);
   set(gca,'XScale','log');
  
   figure;
   plot(w,RGA_no,'.-');
   xlabel('Frequency (Hz)');
   ylabel('||RGA - I|| sum');
   title('RGA Number');
   set(gca, 'XScale','log');   