function [A, B, C, D, G, H] = func_linmats_all(parameters, IC)
% Get linearized matrices for UT hMCL
%{
----------------------------- Description ---------------------------------
Function takes all parameters of the common MCL, and returns the linearized
system matrices. The system matrices will make up the state space
representation of the hMCL as follows : 
    xdot = Ax + Bu + Gd
    y    = Cx + Du + Hd
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

u represents inputs : u = [Flv; Fao; Qrc]

d represents the disturbance to the system, which is the Qvad (m^3/s)

---------------------------- Inputs ---------------------------------------
parameters : Array of hMCL parameters
    1. mlv   : Mass of the LV side piston connector (kg)
            2. k     : Spring constant for the spring between piston connector and
                       piston head (N/m)
            3. mplv  : Mass of LV piston (kg)
            4. Ap    : Area of piston head (m^2)    
            5. Cd    : Compliance of the rolling diaphragm (N/m)
            6. Rlvo  : Orifice resistance on LV side (Pa/m^3/s)    
            7. Raoo  : Orifice resistance on AO side (Pa/m^3/s)    
            8. mpao  : Mass of AO piston (kg)
            9. mao   : Mass of AO side piston connector (kg)
            10. V0   : Initial volume of air pocket (m^3)
            11. Patm : Atmospheric pressure (Pa)
            12. k2   : Parameter for Qrc
            13. k3   : Parameter for Qrc
            14. Rplv : Friction resistance at LV piston and diaphragm
                       interface (Pa/m^3/s)
            15. Rpao : Friction resistance at AO piston and diaphragm
                       interface (Pa/m^3/s)

IC  : Initial conditions around which we linearize   
    
----------------------------- Outputs -------------------------------------
System matrices A, B, C, D, G, H

----------------------------- Versions ------------------------------------
v1 : Suraj R Pawar, 6-5-2020
    - Initialize
%}

% Extract Parameters
    mlv = parameters(1);
    k = parameters(2);
    mplv = parameters(3);
    Ap = parameters(4);    
    Cd = parameters(5);
    Rlv0 = parameters(6);    
    Rao0 = parameters(7);    
    mpao = parameters(8);
    mao = parameters(9);
    V0 = parameters(10);
    Patm = parameters(11);
    Rplv = parameters(14);
    Rpao = parameters(15);
    
% Initial Conditions
    vlv = IC(1);
    vao  = IC(2);
    xslv = IC(3);
    xsao = IC(4);
    vplv = IC(5);
    vpao = IC(6);
    xplv = IC(7);
    xpao = IC(8);
    xdlv = IC(9);
    xdao = IC(10);
    Vlv = IC(11);
    Vao = IC(12);
    Qrc = 2e-5;
    Qvad = 2e-5;
    
% Load file
    filename = 'common_mcl_linsys.mat';
    common_mcl_linsys = load(filename);
    
% Linear System

    A = double(subs(common_mcl_linsys.A));
    B = double(subs(common_mcl_linsys.B));
    C = double(subs(common_mcl_linsys.C));
    D = double(subs(common_mcl_linsys.D));
    G = double(subs(common_mcl_linsys.G));
    H = double(subs(common_mcl_linsys.H));   
end

