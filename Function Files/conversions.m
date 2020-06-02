%% Conversion for variables
%{
!!! Only call this script from within the WCS Simulation file !!!
------------------------- Description -------------------------------------
Convert all user input values to SI units. Also declare necessary constants

----------------------------- Versions ------------------------------------
v1 : Suraj R Pawar, 6-1-2020
    - Initialize
%}

% ------------------------ Constants --------------------------------------
    inch_to_m = 0.0254; 
    mmHg_to_Pa = 133.322;
    cm_to_m = 0.01;
    LPM_to_m3s = 1.6667e-5;
    mL_to_m3 = 1e-6;
    density_blood_plasma = 1025;    % (kg/m^3)        
    Patm = 101325;                  % (Pa)

% ----------------- Convert actuator control gains ------------------------
    set_lv = set_lv_mmHg * mmHg_to_Pa;
    set_ao = set_ao_mmHg * mmHg_to_Pa;

% -------- Piston assembly parameter calculations (both sides) ------------
    % Geometry
    dia_piston = dia_piston_inch * inch_to_m;
    radius_piston = dia_piston/2;
    area_piston = pi * radius_piston^2;        
    Ap = area_piston;                               % Use this parameter for both LV and AO side
    fluid_section = fluid_section_cm * cm_to_m;     % This is the fluid that always sits in the glass chambers        
    vol_fluid_section = area_piston * fluid_section;

    % Mass
    m_fluid = vol_fluid_section * density_blood_plasma;
    mplv = m_piston_head + m_fluid;                 % Mass of piston head and fluid section (kg)

    % Friction resistance
    Qloss_piston = Qloss_piston_LPM * LPM_to_m3s;
    Ploss_piston = Ploss_piston_mmHg * mmHg_to_Pa;
    Rplv = (Ploss_piston / Qloss_piston)*Ap^2;             % Resistance of piston assembly (Pa-s/m^3)

% --------------------------- Air pocket ----------------------------------
    V0 = V0_mL * mL_to_m3;

% -------------------- Set variables for AO side --------------------------

    % AO Voice Coil actuator
    Lvcao = Lvclv;
    Rvcao = Rvclv;

    % Piston assembly
    mpao = mplv;
    Rpao = Rplv;