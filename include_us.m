% Asks MATLAB to add listed folders to its path
%{
v1 : Suraj Pawar, 6-1-2020
    - Initialize
%}

folders = {
        'Function Files';   
        'Data Files';
        'Results';
        };
    
for i = 1 : length(folders)
    addpath(genpath(folders{i}));
end