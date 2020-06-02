function [xdot] = model_vc(t,x,parameters,inputs)
% Function file for voice coil model
% Test this model before integrating equations into bigger model

%% Extract parameters
    Lvclv = parameters(1);
    Rvclv = parameters(2);
    rvc = parameters(3);
    mplv = parameters(4);
    Rplv = parameters(5);

%% Inputs
    u = inputs(1);
    F = inputs(2);
    
%% Extract states
    ilv = x(1);
    vlv = x(2);
    xlv = x(3);
    
%% State equations
    ilvdot = (1/Lvclv)*(u - ilv*Rvclv - rvc*vlv);
    vlvdot = (1/mplv)*(rvc*ilv - Rplv*vlv - F);
    xlvdot = vlv;
    
    xdot = [ilvdot; vlvdot; xlvdot];
end

