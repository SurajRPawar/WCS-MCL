function [A, B, C, D, G, H] = func_linmats(parameters, IC)
% Get linearized matrices for WCS hMCL
%{
----------------------------- Description ---------------------------------
Function takes all parameters of the WCS hMCL, and returns the linearized
system matrices. The system matrices will make up the state space
representation of the hMCL as follows : 
    xdot = Ax + Bu + Gd
    y    = Cx + Du + Hd
Where, x represents states [ilv; vlv; iao; vao; Vlv; Vao; dx]
    ilv : Current in LV Voice Coil actuator (A)
    vlv : Velocity of LV side piston head (m/s)
    iao : Current in AO Voice Coil actuator (A)
    vao : Velocity of AO side piston head (m/s)
    Vlv : Volume in LV capacitor tank (m^3)
    Vao : Volume in AO capacitor tank (m^3)
    dx  : Difference in position of pistons = xlv - xao (m)

u represents inputs : u = [ulv; uao; urc], where each u is the voltage to
the respective actuator (V)

d represents the disturbance to the system, which is the Qvad (m^3/s)

Notes : 
    - Gear pump equation used is 
      Qrc = k2*urc - k3*(Plv - Pao)

---------------------------- Inputs ---------------------------------------
parameters : Array of hMCL parameters
            1. Lvclv  : LV voice coil inductance (H)
            2. Lvcao  : LV voice coil inductance (H)
            3. Rvclv  : LV voice coil resistance (Ohm)
            4. Rvcao  : AO voice coil resistance (Ohm)
            5. rvc    : LV and AO voice coil gyration constant (N/A)
            6. mplv   : Mass of LV piston (kg)
            7. mpao   : Mass of AO piston (kg)
            8. Rplv   : Friction resistance on LV piston (N/m/s)
            9. Rpao   : Friction resistance on AO piston (N/m/s)
            10. Ap    : Area of pistons on LV and AO side (m^2) 
            11. k2    : Constant for gear pump equation
            12. k3    : Constant for gear pump equation
            13. Patm  : Atmospheric pressure constant (Pa)
            14. V0    : Initial volume of air in capacitor tanks (m^3)
            15. Cd    : Diaphragm compliance (N/m)

IC          : Initial conditions around which we linearize
            1. ilv
            2. vlv
            3. xlv
            4. iao
            5. vao
            6. xao
            7. Vlv
            8. Vao

----------------------------- Outputs -------------------------------------
System matrices A, B, C, D, G, H

----------------------------- Versions ------------------------------------
v1 : Suraj R Pawar, 6-1-2020
    - Initialize
v2 : Suraj R Pawar, 6-1-2020
    - Add diaphragm compliance
%}

% A Matrix
    Ap = parameters(10);
    Lvcao = parameters(2);
    Lvclv = parameters(1);
    Patm = parameters(13);
    Rpao = parameters(9);
    Rplv = parameters(8);
    Rvcao = parameters(4);
    Rvclv = parameters(3);
    V0 = parameters(14);
    mpao = parameters(7);
    mplv = parameters(6);
    rvc = parameters(5);
    k3 = parameters(12);
    Cd = parameters(15);
    Vao = IC(8);
    Vlv = IC(7);


    A = [-Rvclv/Lvclv, -rvc/Lvclv, 0, 0, 0, 0, 0, 0;
              rvc/mplv, -Rplv/mplv, -Cd/mplv, 0, 0, 0, -(Ap*Patm*V0)/(mplv*(V0 - Vlv)^2), 0;
              0, 1, 0, 0, 0, 0, 0, 0;
              0, 0, 0, -Rvcao/Lvcao, -rvc/Lvcao, 0, 0, 0;
              0, 0, 0, rvc/mpao, -Rpao/mpao, -Cd/mpao, 0, -(Ap*Patm*V0)/(mpao*(V0 - Vao)^2);
              0, 0, 0, 0, 1, 0, 0, 0;
              0, Ap, 0, 0, 0, 0, -(Patm*V0*k3)/(V0 - Vlv)^2, (Patm*V0*k3)/(V0 - Vao)^2;
              0, 0, 0, 0, Ap, 0, (Patm*V0*k3)/(V0 - Vlv)^2, -(Patm*V0*k3)/(V0 - Vao)^2];

% B Matrix
    k2 = parameters(11);

    B = zeros(8,3);
    B(1,1) = 1/Lvclv;
    B(4,2) = 1/Lvcao;
    B(7,3) = k2;
    B(8,3) = -k2;

% C Matrix
    C = zeros(3,8);
    C(1,7) = (Patm*V0)/(V0 - Vlv)^2;
    C(2,8) = (Patm*V0)/(V0 - Vao)^2;
    C(3,3) = 1;
    C(3,6) = -1;

% D Matrix
    D = zeros(3,3);

% G Matrix
    G = zeros(8,1);
    G(7,1) = -1;
    G(8,1) = 1;

% H Matrix
    H = zeros(3,1);
end

