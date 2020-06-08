%% Conversion for variables
%{
!!! Only call this script from within the Common MCL Simulation file !!!
------------------------- Description -------------------------------------
Convert all user input values to SI units. Also declare necessary constants

----------------------------- Versions ------------------------------------
v1 : Suraj R Pawar, 6-6-2020
    - Initialize
v2 : Suraj R Pawar, 6-6-2020
    - Added friction resistance at piston and diaphragm interface
v3 : Suraj R Pawar, 6-7-2020
    - Fixed unit conversion of actuation force
%}

% ------------------------ Constants --------------------------------------
    inch_to_m = 0.0254; 
    mmHg_to_Pa = 133.322;
    cm_to_m = 0.01;
    LPM_to_m3s = 1.6667e-5;
    mL_to_m3 = 1e-6;
    density_blood_plasma = 1025;    % (kg/m^3)        
    Patm = 101325;                  % (Pa)

% ----------------- Actuators ---------------------------------------------
    freq_lv = 2 * pi * freq_lv_Hz;
    freq_ao = 2 * pi * freq_ao_Hz;
    
% -------- Piston assembly parameter calculations (both sides) ------------
    % Geometry
    dia_piston = dia_piston_inch * inch_to_m;
    radius_piston = dia_piston/2;
    area_piston = pi * radius_piston^2;        
    Ap = area_piston;                               % Use this parameter for both LV and AO side            
    
    % Friction resistance
    Qloss_piston = Qloss_LPM * LPM_to_m3s;
    Ploss_piston = Ploss_mmHg * mmHg_to_Pa;
    Rplv = (Ploss_piston / Qloss_piston)*Ap^2;             % Resistance of piston assembly (Pa-s/m^3)
    
% -------------------------- Nonlinear orifice ----------------------------
    do = do_cm * cm_to_m;
    ro = do/2;
    area_orifice = pi * ro^2;
    R_orifice = density_blood_plasma / (2 * cdo^2 * area_orifice^2); % Orifice flow = sqrt(P/R)
    
% --------------------------- Air pocket ----------------------------------
    V0 = V0_mL * mL_to_m3;

% ----------------------------- Set variables -----------------------------
    Rlvo = R_orifice;
    Raoo = R_orifice;
    Rpao = Rplv;