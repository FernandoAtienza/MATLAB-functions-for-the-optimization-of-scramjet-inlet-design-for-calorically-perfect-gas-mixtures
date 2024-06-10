% Ideal/required temperature at the scramjet combustion chamber entrance
% for calorically perfect gases
%
% @author: Fernando Atienza González
%          Senior Undergraduate Aerospace Engineer Student
%          Universidad Carlos III de Madrid
%
% Last updated: 7th June 2024
% -------------------------------------------------------------------------

function [T_id, g, g_f, R_f] = Ideal_Temp(f,self)

% Air properties
R_a = 287.05;           % Air Gas Constant
g = 1.4;                % Air Specific Heat Ratio
T_inj = 300;            % Injection Temperature [K]

f_Rg = @(species) self.C.R0 / (self.DB.(species).mm * 1e-3);

% Ratio of specific heats and specific gas constant data obtained from
% Combustion turbo toolbox (by Alberto Cuadra Lara)
    switch f
        case 1  % Hydrogen (H2)
            g_f = species_gamma('H2', T_inj, self.DB);              % H2 Specific heat ratio
            R_f = f_Rg('H2');                                       % H2 Gas Constant
            f_st = (3*2)/(103*2);                                   % Fuel(H2)/air mass flow ratio
            T_ig = 845.15;                                          % H2 Ignition Temperature [K]

        case 2  % Methane (CH4)
            g_f = species_gamma('CH4', T_inj, self.DB);             % Methane Specific heat ratio
            R_f = f_Rg('CH4');                                      % Methane Gas Constant
            f_st = (36 + 3*4)/(103*(4+4));                          % Fuel (Methane)/air mass flow ratio
            T_ig = 810.15;                                          % Methane Ignition Temperature [K]
    
        case 3  % Ethane (C2H6)
            g_f = species_gamma('C2H6', T_inj, self.DB);            % Ethane Specific heat ratio
            R_f = f_Rg('C2H6');                                     % Ethane Gas Constant
            f_st = (36*2 + 3*6)/(103*(4*2 + 6));                    % Fuel (Ethane)/air mass flow ratio 
            T_ig = 745.15;                                          % Ethane Ignition Temperature [K]
    
        case 4  % Octane (C8H18)
            g_f = species_gamma('C8H18_isooctane', T_inj, self.DB); % Octane Specific heat ratio
            R_f = f_Rg('C8H18_isooctane');                          % Octane Gas Constant
            f_st = (36*8 + 3*18)/(103*(4*8+18));                    % Fuel (Octane)/air mass flow ratio 
            T_ig = 479.15;                                          % Octane Ignition Temperature [K]
    end

    Cp_f = ( g_f/(g_f-1) ) * R_f;           % Fuel Specific Heat
    Cp_a = (g/(g-1)) * R_a;                 % Air Specific Heat

    T_id = f_st * (Cp_f/Cp_a) * (T_ig - T_inj) + T_ig;      % Ideal Temperature at the Combustion Chamber

end