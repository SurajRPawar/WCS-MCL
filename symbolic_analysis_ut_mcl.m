%% Symbolic analysis for UT hMCL model
%{
------------------------- Description -------------------------------------
Symbolic Analysis for UT MCL Model. Provide nonlinear equations and get
linearized matrices. Also performs symbolic controllability and
observability analysis (with linearized system matrices)

Parameters : 
    1. Lvcao : Inductance for the AO Voice coil (H)
    2. Rvcao : Resistance for AO voice coil (Ohms)
    3. rao   : Gyration constant for AO voice coil (N/A)
    4. ktao    : Spring constant for AO side top spring (N/m)
    5. mao   : Mass of AO piston connector(kg)
    6. kmao    : Spring constant for middle spring, the one between the
               voice coil piston and the end effector (N/m)
    7. mpao  : Mass of the end effector on the AO side (kg)
    8. Rpao  : Friction loss on the piston on AO side (N/m/s)
    9. Cd    : Spring compliance of the rolling diaphragm (N/m)
    10. Ap   : Area of the end effector piston for both sides (m^2)
    11. Cao  : Tank capacitance on AO side (m^3/Pa)
    12. Lvclv : Inductance for the LV voice coil (H)
    13. Rvclv : Resistance for LV voice coil (Ohm)
    14. rlv   : Gyration constant for LV voice coil (N/A)
    15. mlvt : Mass of LVT transducer on LV side (kg)
    16. kmlv : Spring compliance of top spring on LV side (N/m)
    17. mplv : Mass of LV piston (kg)
    18. Rplv : Friction loss on LV piston (N/m/s)
    19. Clv : Tank capacitance on LV side (m^3/Pa)
    20. Rvad : Loss in VAD (Pa/m^3/s)
    21. Rrc : Loss in RC (Pa/m^3/s)

States : x = [iao ilv xtao vao xao vlvt xlvt xmao xmlv vpao xpao vplv xplv
              xdao xdlv Vao Vlv]^T
    1. iao : Current in AO VC (A)
    2. ilv : Current in LV VC (A)
    3. xtao : Compression in top spring on AO side (m)
    4. vao : Velocity of AO piston connector (m/s)  
    5. vlvt : Velocity of LVT transducer (m/s)
    6. xlvt : Position of LVT transducer (m)
    7. xmao : Compression in middle spring on AO side (m)
    8. xmlv : Compression in middle spring on LV side (m)
    9. vpao : Velocity of AO piston (m/s)
    10. xpao : Position of AO piston (m)
    11. vplv : Velocity of LV piston (m/s)
    12. xplv : Position of LV piston (m)
    13. xdao : Compression in rolling diaphragm on AO side (m)
    14. xdlv : Compression in rolling diaphragm on LV side (m)
    15. Vao : Volume in AO tank (m^3)
    16. Vlv : Volume in LV tank (m^3)

------------------------- Versions ----------------------------------------
v1 : Suraj R Pawar, 6-3-2020
    - Initialize
%}

clear all; close all; clc;
savefile = 1;
filename = 'Data Files/ut_mcl_linsys.mat';

%% Define symbolics
    fprintf('Defining Symbols\n');
    % Parameters
        % LV Side
        syms Rvcao Lvcao rao ktao mao kmao Rpao mpao Cao real

        % LV side
        syms mlvt rlv Rvclv Lvclv kmlv mplv Rplv Clv real

        % RC
        syms Rrc real

        % VAD
        syms Rvad real
        
        % Common
        syms g Cd Ap real
    
    % States
        syms iao ilv xtao vao vlvt xlvt xmao xmlv vpao xpao vplv xplv real
        syms xdao xdlv Vao Vlv real
    
    % Inputs
        syms Qrc uao ulv real
        
    % Disturbance 
        syms Qvad real                    

%% State equations
    fprintf('Forming state equations \n');
    Fmao = kmao*xmao;
    Fmlv = kmlv*xmlv;
    Fdao = Cd*xdao;
    Fdlv = Cd*xdlv;   
    Pao = Vao/Cao;
    Plv = Vlv/Clv;
    P17 = Pao + Fdao/Ap;
    P28 = Plv + Fdlv/Ap;
    Q26 = (1/Rrc)*(P28 - P17);
    Q21 = (1/Rvad)*(P17 - P28);
    
    %{
        x = [iao ilv xtao vao vlvt xlvt xmao xmlv vpao xpao vplv xplv
              xdao xdlv Vao Vlv]^T
    %}
    iaodot = (1/Lvcao)*(uao - Rvcao*iao - rao*vao);
    ilvdot = (1/Lvclv)*(ulv - Rvclv*ilv - rlv*vlvt);
    xtaodot = vao;
    vaodot = (1/mao)*(rao*iao - ktao*xtao - Fmao);
    xaodot = vao;
    vlvtdot = (1/mlvt)*(rlv*ilv - Fmlv);
    xlvtdot = vlvt;    
    xmaodot = vao - vpao;
    xmlvdot = vlvt - vplv;
    vpaodot = (1/mpao)*(Fmao - Rpao*vpao - Fdao);    
    xpaodot = vpao;
    vplvdot = (1/mplv)*(Fmlv - Rplv*vplv - Fdlv);
    xplvdot = vplv;
    Vaodot = Qvad - Q21 - (Qrc - Q26);
    Vlvdot = Qrc - Q26 - Qvad + Q21;
    xdaodot = vpao + (1/Ap)*(Qvad - Q21 - Qrc + Q26);
    xdlvdot = vplv + (1/Ap)*(Qrc - Q26 - Qvad + Q21);
       
    f = [iaodot;
         ilvdot;
         xtaodot;
         vaodot;
         vlvtdot;
         xlvtdot;
         xmaodot;
         xmlvdot;
         vpaodot;
         xpaodot;
         vplvdot;
         xplvdot;
         xdaodot;
         xdlvdot;
         Vaodot;
         Vlvdot]; 

%% Outputs
    fprintf('Creating output variables \n');
    y = [P28;   % LV Side
         P17   % AO Side
         xplv - xpao];
 
%% Define states, inputs and disturbance
    fprintf('Computing Jacobians \n');
    x = [iao; 
        ilv; 
        xtao; 
        vao; 
        vlvt;   
        xlvt;
        xmao; 
        xmlv; 
        vpao;   
        xpao;
        vplv;
        xplv;
        xdao; 
        xdlv; 
        Vao; 
        Vlv];
    u = [Qrc; uao; ulv];
    d = Qvad;

    % Linearized matrices
    A = jacobian(f,x);
    B = jacobian(f,u);
    G = jacobian(f,d);
    C = jacobian(y,x);
    D = jacobian(y,u);
    H = jacobian(y,d);

%% Data Logging
    if savefile == 1
        save(filename, 'A', 'B', 'C', 'D', 'G', 'H');
    end
