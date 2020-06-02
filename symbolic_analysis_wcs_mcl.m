%% Symbolic analysis for WCS hMCL model
%{
------------------------- Description -------------------------------------
Symbolic Analysis for WCS MCL Model. Provide nonlinear equations and get
linearized matrices. Also performs symbolic controllability and
observability analysis (with linearized system matrices)

------------------------- Versions ----------------------------------------
v1 : Suraj R Pawar, 5-18-2020
    - Initialize
v2 : Suraj R Pawar, 6-1-2020
    - Add diaphragm compliance
%}

clear all; close all; clc;

% Define symbolics
syms Lvclv Lvcao Rvclv Rvcao rvc mplv mpao Rplv Rpao Ap k2 k3 Patm V0 Cd real
syms ilv vlv iao vao Vlv Vao xlv xao real
syms urc ulv uao real
syms Qvad real

% Calculate quantities needed for state dynamics
Plv = Patm*(V0/(V0 - Vlv) - 1);    
Pao = Patm*(V0/(V0 - Vao) - 1);    
Qrc = k2*urc - k3*(Plv - Pao);      

% State Equations and outputs
ilvdot = (1/Lvclv)*(ulv - Rvclv*ilv - rvc*vlv);
vlvdot = (1/mplv)*(rvc*ilv - Rplv*vlv - Ap*Plv - Cd*xlv);
xlvdot = vlv;
Vlvdot = vlv*Ap + Qrc - Qvad;
iaodot = (1/Lvcao)*(uao - Rvcao*iao - rvc*vao);
vaodot = (1/mpao)*(rvc*iao - Rpao*vao - Ap*Pao - Cd*xao);
xaodot = vao;
Vaodot = vao*Ap - Qrc + Qvad;

f = [ilvdot;
     vlvdot; 
     xlvdot
     iaodot;
     vaodot;  
     xaodot
     Vlvdot;
     Vaodot]; 

y = [Plv; 
     Pao;
     (xlv - xao)];
 
% Define states, inputs and disturbance
x = [ilv; vlv; xlv; iao; vao; xao;  Vlv; Vao];
u = [ulv; uao; urc];
d = Qvad;

% Linearized matrices
A = jacobian(f,x);
B = jacobian(f,u);
G = jacobian(f,d);
C = jacobian(y,x);
D = jacobian(y,u);
H = jacobian(y,d);

% Analysis
num_states = numel(x);
num_ins = numel(u);
num_meas = numel(y);

WC = sym(zeros(num_states, num_states*num_ins));
for i = 1:num_states
    start_index = num_ins*(i-1) + 1;
    stop_index = num_ins*i;    
    WC(:,[start_index : stop_index]) = (A^(i-1))*B;
end
rank_controllability = rank(WC)

WO = sym(zeros(num_states*num_meas, num_states));
for i = 1:num_states
    start_index = num_meas*(i-1) + 1;
    stop_index = num_meas*i;    
    WO([start_index : stop_index],:) = C*A^(i-1);
end

rank_observability = rank(WO)

