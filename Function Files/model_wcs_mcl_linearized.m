function [xdot,y, sys] = model_wcs_mcl_linearized(t,x,parameters,inputs,IC)
% Funciton file for Linearized MCL (Windmill Cardiovascular Systems)
%{
------------------------------ Description --------------------------------
State equuations for the WCS MCL. The model file must be called using a
solver in a main simulation file. This model contains linearized state
equations.

--------------------------------- Inputs ----------------------------------
t         :  Time (s)
x         :  States
             ilv : Current in LV Voice Coil actuator (A)
             vlv : Velocity of LV side piston head (m/s)
             xlv : Position of the LV piston head (m)
             iao : Current in AO Voice Coil actuator (A)
             vao : Velocity of AO side piston head (m/s)
             xao : Position of the AO piston head (m)
             Vlv : Volume in LV capacitor tank (m^3)
             Vao : Volume in AO capacitor tank (m^3)
IC        : Initial conditions used for linearizing system matrices
parameters : Lvclv, Lvcao : LV and AO VC Inductance (H)
             Rvclv, Rvcao : LV and AO VC Resistances (Ohm)
             rvc : LV and AO VC Solenoid gyration constant (N/A)
             mplv, mpao : Mass of LV and AO Piston heads including 
                          fluid section (kg)
             Rplv, Rpao : Friction resistance on piston heads (N/m/s)
             Ap : Area of piston head (m^2)
             k2, k3 : Constants for RC pump Q-H curve. 
                      Q = k2*urc - k3*H (m^3/s), where
                      urc = RC pump input voltage (V)
                      H = Pressure head across pump (Pa)
             Patm : Atmospheric pressure (Pa)
             V0 : Initial volume of air in LV and AO capacitor tanks (m^3)
             Cd : Compliance of diaphragm (N/m)

------------------------ Outputs ------------------------------------------
Qvad : VAD flow rate (m^3/s)
Qrc : RC flow rate (m^3/s)
Plv : LV pressure (Pa)
Pao : AO pressure (Pa)

---------------------------- Versions -------------------------------------
v1 : Suraj R Pawar, 5-15-2020
    - Initialize
v2 : Suraj R Pawar, 6-1-2020
    - Include state equations due to diaphragm
%}

%% Extract Parameters
    Lvclv = parameters(1);
    Lvcao = parameters(2);
    Rvclv = parameters(3);
    Rvcao = parameters(4);
    rvc = parameters(5);
    mplv = parameters(6);
    mpao = parameters(7);
    Rplv = parameters(8);
    Rpao = parameters(9);
    Ap = parameters(10);
    k2 = parameters(11);
    k3 = parameters(12);
    Patm = parameters(13);
    V0 = parameters(14);
    Cd = parameters(15);
    
%% Extract states
    ilv = x(1);
    vlv = x(2);
    xlv = x(3);
    iao = x(4);
    vao = x(5);
    xao = x(6);
    Vlv = x(7);
    Vao = x(8);

%% Inputs   
    ulv = inputs(1);
    uao = inputs(2);
    urc = inputs(3);
    u = inputs;
    
%% Compute other signals
    % Pressures in LV and AO
    Plv = Patm*(V0/(V0 - Vlv) - 1);    % Gauge LV pressure (Pa)
    Pao = Patm*(V0/(V0 - Vao) - 1);    % Gauge AO pressure (Pa)
    
    % Qrc
    Qrc = k2*urc - k3*(Plv - Pao);      
    
    % Qvad
    SV = 30;    % Stroke volume (mL)
    T = 0.5;    % Ejection time (sec)
    Qvad = SV*2*60/(1000*T)*sin(2*pi*t/(2*T))^2/60000;  % (m^3/s)
    
    % State Space Representation
    %{
        xdot = A*x + B*u + G*Qvad
        y    = C*x + D*u + H*Qvad
    %}
    
    [A, B, C, D, G, H] = func_linmats(parameters,IC);
    xdot = A*x + B*u + G*Qvad;
    
%% Outputs
    y = [Qvad; Qrc; Plv; Pao];
    % Form the structure sys containing all linearized system matrices
    sys.A = A;
    sys.B = B;
    sys.G = G;
    sys.C = C;
    sys.D = D;
    sys.H = H;
end