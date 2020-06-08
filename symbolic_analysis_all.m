%% Symbolic analysis for all MCL models
%{
------------------------- Description -------------------------------------
Symbolic Analysis for a common MCL model that incorporates all the elements
found in both UT and WCS MCL. The objective is to do a parametric study to 
analyze which elements contribute to the decoupling. 
Provide nonlinear equations and get linearized matrices.

Parameters  
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
    12. Rplv : Friction resistance at LV piston and diaphragm
                       interface (Pa/m^3/s)
    13. Rpao : Friction resistance at AO piston and diaphragm
               interface (Pa/m^3/s)

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

------------------------- Versions ----------------------------------------
v1 : Suraj R Pawar, 6-5-2020
    - Initialize
v2 : Suraj R Pawar, 6-6-2020
    - Removed the crevice type model of the rolling diaphragm
    - The crevice approach complicates things unnecessarily 
v3 : Suraj R Pawar, 6-6-2020
    - Added friction loss at the piston diaphragm interface    
%}

clear all; close all; clc;
savefile = 1;
filename = 'Data Files/common_mcl_linsys.mat';

%% Define symbolics
    fprintf('Defining Symbols\n');
    % Parameters
        % LV Side
        syms mlv mplv Rlvo Rplv real

        % AO side
        syms mao mpao Raoo Rpao real

        % Inputs
        syms Qrc Flv Fao real

        % VAD
        syms Qvad real
        
        % Common
        syms g k Ap Cd Patm V0 real
    
    % States
        syms vlv vao xslv xsao vplv vpao xplv xpao xdlv xdao Vlv Vao real              

%% State equations
    fprintf('Forming state equations \n');
    
    % Variables used in state equations
        % LV Side
        Fslv = k*xslv;
        Fdlv = Cd*xdlv;        
        Ptlv = Patm*(V0/(V0 - Vlv) - 1); % LV vertical tank pressure, due to air pocket
        Qtlv = vplv*Ap + Qrc - Qvad;     % Flow into the vertical tube
        Prlv = Rlvo * (Qtlv)^2;          % Pressure admitted by nonlinear orifice
        Plv = Ptlv + Prlv;
        
        % AO side
        Fsao = k*xsao;
        Fdao = Cd*xdao;
        Qtao = vpao*Ap + Qvad - Qrc;        
        Ptao = Patm*(V0/(V0 - Vao) - 1); % AO vertical tank pressure, due to air pocket
        Prao = Raoo * (Qtao)^2;
        Pao = Ptao + Prao;
        
        
    % LV side
    vlvdot = (1/mlv)*(Flv - Fslv);
    xslvdot = vlv - vplv;
    vplvdot = (1/mplv)*(Fslv - Plv*Ap - Cd*xdlv - Rplv*vplv);
    xplvdot = vplv;
    xdlvdot = vplv;
    Vlvdot = vplv*Ap + Qrc - Qvad;    
    
    % AO side
    vaodot = (1/mao)*(Fao - Fsao);
    xsaodot = vao - vpao;
    vpaodot = (1/mpao)*(Fsao - Pao*Ap - Cd*xdao - Rpao*vpao);
    xpaodot = vpao;
    xdaodot = vpao;
    Vaodot = vpao*Ap + Qvad - Qrc;    
    
    f = [vlvdot;
         vaodot;
         xslvdot;
         xsaodot;
         vplvdot;
         vpaodot;
         xplvdot;
         xpaodot;
         xdlvdot;
         xdaodot;
         Vlvdot;
         Vaodot];        

%% Outputs
    fprintf('Creating output variables \n');
    y = [Plv; % LV Side
         Pao;   % AO Side
         xplv - xpao];
 
%% Define states, inputs and disturbance
    fprintf('Computing Jacobians \n');
    x = [vlv;
         vao;
         xslv;
         xsao;
         vplv;
         vpao;
         xplv;
         xpao;
         xdlv;
         xdao;
         Vlv;
         Vao];
     
    u = [Flv; Fao; Qrc];
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
        fprintf('Measurement filed saved with filename: %s \n', filename);
    end
