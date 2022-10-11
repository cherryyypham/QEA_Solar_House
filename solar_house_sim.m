function [Time,T] = solar_house_sim()
    % HOUSE CONSTRUCTION PARAMETERS
    B_HEIGHT = 3.2; % Wall Height (m)
    b_depth = 2.2; % Building Depth (m) - calculated
    B_LENGTH = 5.1; % Building Length (m)
    TOP_PANE_HEIGHT = 0.4;
    BOTTOM_PANE_HEIGHT = 0.2;
    SA_PANES = TOP_PANE_HEIGHT * b_depth + BOTTOM_PANE_HEIGHT * b_depth;
    HEIGHT_WINDOW = 2.6; % Window Height (m)
    b_thickness = 100; % Wall thickness (m) - sweeping through
    esu_thickness = 100; % Storage unit thickness (m) - sweeping through
    
    % THERMAL-ODE-RELATED PARAMETERS
    k_walls = 0.04; % Conductivity of fiberglass insulation (W/m-K)
    p_esu = 3000; % Energy Storage Unit (ESU) Density (kg/m^3)
    c_esu = 800;  % ESU 
    h_indoors = 15;
    h_outdoors = 30;
    T_air = -3;
    T_storage_i = 20; % Randomly picking from given range of ()
    h_window = 1.4;
    
    % CALCULATING AREAS
    SA_esu = (B_LENGTH - 2 * b_thickness) * b_depth;
    SA_walls_in = 2*((B_LENGTH - 2 * b_thickness) * b_depth) + B_HEIGHT * b_depth + SA_PANES;
    SA_walls_out = 2 * B_LENGTH * b_depth + B_HEIGHT * b_depth + SA_PANES + 4 * b_thickness;
    SA_window = HEIGHT_WINDOW * b_depth;
    
    % ODE SUBEQUATIONS
    % Resistance
    R_air_in = 1./(h_indoors * SA_esu);
    R_floor_cond = b_thickness./(k_walls * SA_walls_in);
    R_floor_conv = 1./(h_indoors * SA_walls_in);
    R_window = 1./(h_window * SA_window);
    R_air_out = 1./(h_outdoors * SA_walls_out);
    R = R_air_in + 1./(1./(R_floor_cond + R_floor_conv) + 1./(R_window)) + R_air_out;
    
    % Time values in hours
    t_i = 0;
    t_f = 24 * 3600;
    
    % Initial Temperature
    T_i = [T_air T_storage_i];
    
    % ODE CALCULATION
    [Time, T] = ode45(@rate_func, [t_i t_f], T_i);
    
    
    function res = rate_func(t, T_i)  
        % Solar flux
        q = -1 * cos((pi * t) / (12 * 3600)) + 224 * cos((pi * t) / (6 * 3600)) + 210;
        q_in = q * SA_window;
        q_out = (T_i(0) - T_i(1)) / R;
        dTdt = q_in/C - q_out/C;
        res = dTdt;
    end
end
