function [xdot,y, sys] = model_common_mcl_linearized(t,x,parameters,inputs,mats)
% Funciton file for Linearized MCL (Windmill Cardiovascular Systems)
%{
------------------------------ Description --------------------------------
State equuations for the common MCL. The model file must be called using a
solver in a main simulation file. This model contains linearized state
equations.

--------------------------------- Inputs ----------------------------------
t         :  Time (s)
x         :  x = [vlv vao xslv xsao vplv vpao xplv xpao xdlv xdao Vlv Vao]^T
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
mats        : Linearized system matrices
parameters : 
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

------------------------ Outputs ------------------------------------------
Qvad : VAD flow rate (m^3/s)
Qrc : RC flow rate (m^3/s)
Plv : LV pressure (Pa)
Pao : AO pressure (Pa)

---------------------------- Versions -------------------------------------
v1 : Suraj R Pawar, 6-6-2020
    - Initialize
v2 : Suraj R Pawar, 6-6-2020
    - Removed the crevice type model of the rolling diaphragm
    - The crevice approach complicates things unnecessarily 
v3 : Suraj R Pawar, 6-6-2020
    - Added friction loss at the piston diaphragm interface
    - Shifted Qrc from function file to outside, in the simulation file
    - Had to shift Qrc from here because Qrc needs Plv, while Plv needs
    Qrc. Therefore, shifted Qrc calcultion to simulation file so that it
    can use previous values of Plv and Pao
v4 : Suraj R Pawar, 6-7=2020
    - Linearized matrices are fed into function file as input
    - This saves processing by not computing the linearized matrices every
    time
%}

%% Extract Parameters
    mlv = parameters(1);
    k = parameters(2);
    mplv = parameters(3);
    Ap = parameters(4);    
    Cd = parameters(5);
    Rlvo = parameters(6);    
    Raoo = parameters(7);    
    mpao = parameters(8);
    mao = parameters(9);
    V0 = parameters(10);
    Patm = parameters(11);
    k2 = parameters(12);
    k3 = parameters(13);
    Rplv = parameters(14);
    Rpao = parameters(15);
    
%% Extract states
    vlv = x(1);
    vao  = x(2);
    xslv = x(3);
    xsao = x(4);
    vplv = x(5);
    vpao = x(6);
    xplv = x(7);
    xpao = x(8);
    xdlv = x(9);
    xdao = x(10);
    Vlv = x(11);
    Vao = x(12);

%% Inputs   
    Flv = inputs(1);
    Fao = inputs(2);
    Qrc = inputs(3);
    u = inputs;
    
%% Compute other signals         
    % Qvad
    SV = 10;    % Stroke volume (mL)
    T = 0.5;    % Ejection time (sec)
    Qvad = SV*2*60/(1000*T)*sin(2*pi*t/(2*T))^2/60000;  % (m^3/s)
    
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
    
%% State equations
    % State Space Representation
    %{
        xdot = A*x + B*u + G*Qvad
        y    = C*x + D*u + H*Qvad
    %}
        
    xdot = mats.A*x + mats.B*u + mats.G*Qvad;
    
%% Outputs
    y = [Qvad; Qrc; Plv; Pao];
    % Form the structure sys containing all linearized system matrices
    sys.A = mats.A;
    sys.B = mats.B;
    sys.G = mats.G;
    sys.C = mats.C;
    sys.D = mats.D;
    sys.H = mats.H;
end