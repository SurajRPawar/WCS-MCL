%% Symbolic RGA Analysis on common MCL design
%{
-------------------------- Description ------------------------------------
Symbolic RGA Analysis on common MCL Linearized model. To learn more about RGA analysis,
refer to Ch 3, section 3.6 and Ch 10, section 10.8 from Skogestad and
Postlethwaite.

Parameters  
    1. mlv   : Mass of the LV side piston connector (kg)
    2. k     : Spring constant for the spring between piston connector and
               piston head (N/m)
    3. mplv  : Mass of LV piston (kg)
    4. Ap    : Area of piston head (m^2)
    5. Ac    : Area of crevice betwee piston head and rolling diaphragm
                                                                   (m^2)
    6. Cd    : Compliance of the rolling diaphragm (N/m)
    7. Rlv   : Orifice resistance on LV side (Pa/m^3/s)
    8. Clv   : Capacitance of the LV vertical tank with air pocket (m^3/Pa)
    9. Rao   : Orifice resistance on AO side (Pa/m^3/s)
    10. Cao  : Capacitance of the AO vertical tank with air pocket (m^3/Pa)
    11. mpao : Mass of AO piston (kg)
    12. mao  : Mass of AO side piston connector (kg)

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

Outputs : y = [Plv; Pao; xplv - xpao]
Inputs : u = [Flv; Fao; Qrc]

-------------------------- Versions ---------------------------------------
v1 : Suraj R Pawar, 6-5-2020
    - Initialize
%}

close all; clear; clc;
include_us;

%% Declare symbolic parameters

    % LV Side
    syms mlv mplv Rlv real

    % LV side
    syms mao mpao Rao real

    % Inputs
    syms Qrc Flv Fao real

    % VAD
    syms Qvad real

    % Common
    syms g k Ap Ac Cd Patm V0 real
    
    % ------------------- Initial Conditions ------------------------------    
    
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
    
    IC = [vlv0; vao0; xslv0; xsao0; vplv0; vpao0; xplv0; xpao0; xdlv0; xdao0; Vlv0; Vao0];
              
    parameters = [mlv; k; mplv; Ap; Ac; Cd; Rlv; Rao; mpao; mao; V0; Patm];   
                 
    
%% Form linearized system
    [A, B, C, D, ~, ~] = func_linmats_all(parameters, IC);                
    
%% RGA Matrix    
    num_meas = size(C,1);
    num_ins = size(B,2);
    num_states = size(A,1);
    
    syms w  % Frequency in rad/s  
    syms s       
    
    G = C*(inv(s*eye(num_states) - A))*B + D;
    %RGA_mat(:,:,i) = ucrga(G);
    RGA_mat = G.*inv(G).'; % RGA matrix
    RGA_mag = abs(RGA_mat);               
   
%    figure;
%    set(groot,'defaultLineLineWidth',0.8);   
%    
%    subplot(3,1,1);
%    hold on;
%    plot(w,squeeze(RGA_dB(1,1,:)));
%    plot(w,squeeze(RGA_dB(1,2,:)),'rx-','MarkerSize',2);
%    plot(w,squeeze(RGA_dB(1,3,:)),'ko-','MarkerSize',2);
%    hold off;
%    ylabel('Mag (dB)');
%    title('PLV');
%    legend('ulv', 'Qrc','uao');
% %    ylim([-120 20]);
%    set(gca,'XScale','log');
%       
%    subplot(3,1,2);
%    hold on;
%    plot(w,squeeze(RGA_dB(2,1,:)));
%    plot(w,squeeze(RGA_dB(2,2,:)),'rx-','MarkerSize',2);
%    plot(w,squeeze(RGA_dB(2,3,:)),'ko-','MarkerSize',2);
%    hold off;
%    ylabel('Mag (dB)');
%    title('PAO');
%    legend('ulv', 'Qrc','uao');
% %    ylim([-120 20]);
%    set(gca,'XScale','log');
%    
%    subplot(3,1,3);
%    hold on;
%    plot(w,squeeze(RGA_dB(3,1,:)));
%    plot(w,squeeze(RGA_dB(3,2,:)),'rx-','MarkerSize',2);
%    plot(w,squeeze(RGA_dB(3,3,:)),'ko-','MarkerSize',2);
%    hold off;
%    ylabel('Mag (dB)');
%    title('xlvt');
%    xlabel('Frequency (Hz)');
%    legend('ulv', 'Qrc','uao');
% %    ylim([-50 20]);
%    set(gca,'XScale','log');
% 
%    figure;
%    plot(w,RGA_no,'.-');
%    xlabel('Frequency (Hz)');
%    ylabel('||RGA - I|| sum');
%    title('RGA Number');
%    set(gca, 'XScale','log');