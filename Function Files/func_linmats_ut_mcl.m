function [A, B, C, D, G, H] = func_linmats_ut(parameters, IC)
% Get linearized matrices for UT hMCL
%{
----------------------------- Description ---------------------------------
Function takes all parameters of the UT hMCL, and returns the linearized
system matrices. The system matrices will make up the state space
representation of the hMCL as follows : 
    xdot = Ax + Bu + Gd
    y    = Cx + Du + Hd
Where, x represents states [iao ilv xtao vao xao vlvt xlvt xmao xmlv vpao xpao vplv xplv
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

u represents inputs : u = [ulv; uao; Qrc], where each u is the voltage to
the respective actuator (V)

d represents the disturbance to the system, which is the Qvad (m^3/s)

---------------------------- Inputs ---------------------------------------
parameters : Array of hMCL parameters
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

IC  : Initial conditions around which we linearize   
    
----------------------------- Outputs -------------------------------------
System matrices A, B, C, D, G, H

----------------------------- Versions ------------------------------------
v1 : Suraj R Pawar, 6-3-2020
    - Initialize
%}

% Extract Parameters
    Lvcao = parameters(1);
    Rvcao = parameters(2);
    rao = parameters(3);
    ktao = parameters(4);
    mao = parameters(5);
    kmao = parameters(6);
    mpao = parameters(7);
    Rpao = parameters(8);
    Cd = parameters(9);
    Ap = parameters(10);
    Cao = parameters(11);
    Lvclv = parameters(12);
    Rvclc = parameters(13);
    rlv = parameters(14);
    mlvt = parameters(15);
    kmlv = parameters(16);
    mplv = parameters(17);
    Rplv = parameters(18);
    Clv = parameters(19);
    Rvad = parameters(20);
    Rrc = parameters(21);
    
% Initial Conditions
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
    Vlv0 = 0;
    
% Load file
    filename = 'ut_mcl_linsys.mat';
    ut_mcl_linsys = load(filename);
    
% Linear System

    A = double(subs(ut_mcl_linsys.A));
    B = double(subs(ut_mcl_linsys.B));
    C = double(subs(ut_mcl_linsys.C));
    D = double(subs(ut_mcl_linsys.D));
    G = double(subs(ut_mcl_linsys.G));
    H = double(subs(ut_mcl_linsys.H));   
end

