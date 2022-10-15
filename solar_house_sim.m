function [Time,T] = solar_house_sim(w_thickness, esu_thickness)
    % HOUSE CONSTRUCTION PARAMETERS
    W_HEIGHT = 3.2; % Wall Height (m)
    W_DEPTH = 5; % Building Depth (m) - calculated
    W_LENGTH = 5.1; % Building Length (m)
    TOP_PANE_HEIGHT = 0.4;
    BOTTOM_PANE_HEIGHT = 0.2;
    SA_PANES = TOP_PANE_HEIGHT * W_DEPTH + BOTTOM_PANE_HEIGHT * W_DEPTH;
    HEIGHT_WINDOW = 2.6; % Window Height (m)
%     w_thickness = 0.02; % Wall thickness (m) - sweeping through
%     esu_thickness = 0.5; % Storage unit thickness (m) - sweeping through
    
    % THERMAL-ODE-RELATED PARAMETERS
    K_WALLS = 0.04; % Conductivity of fiberglass insulation (W/m-K)
    P_ESU = 3000; % Energy Storage Unit (ESU) Density (kg/m^3)
    C_ESU = 800;  % ESU 
    H_INDOORS = 15;
    H_OUTDOORS = 30;
    T_AIR = -3;
    T_storage_i = 20; % Randomly picking from given range of ()
    H_WINDOW = 1.4;

    v_esu = esu_thickness * W_DEPTH * W_LENGTH;
    C = P_ESU * v_esu * C_ESU;

    % CALCULATING AREAS
    SA_esu = (W_LENGTH - 2 * w_thickness) * W_DEPTH;
    SA_walls_in = 2*((W_LENGTH - 2 * w_thickness) * W_DEPTH) + W_HEIGHT * W_DEPTH + SA_PANES + 2*((W_LENGTH - 2 * w_thickness) * W_HEIGHT);
    SA_walls_out = 2 * W_LENGTH * W_DEPTH + W_HEIGHT * W_DEPTH + SA_PANES + 2 * W_LENGTH * W_HEIGHT;
    SA_window = HEIGHT_WINDOW * W_DEPTH;
    
    % ODE SUBEQUATIONS
    % Resistance
    R_air_in = 1./(H_INDOORS * SA_esu);
    R_floor_cond = w_thickness./(K_WALLS * SA_walls_in);
    R_floor_conv = 1./(H_INDOORS * SA_walls_in);
    R_window = 1./(H_WINDOW * SA_window);
    R_air_out = 1./(H_OUTDOORS * SA_walls_out);
    R = R_air_in + 1./(1./(R_floor_cond + R_floor_conv) + 1./(R_window)) + R_air_out;
    
    % Time values in hours
    t_i = 0;
    t_f = 30 * 24 * 3600;
    
    % Initial Temperature
    
    % ODE CALCULATION
    [Time, T] = ode45(@rate_func, [t_i t_f], T_storage_i);
    
    
    function res = rate_func(t, T_storage_i)  
        % Solar flux
        q = -1 * cos((pi * t) / (12 * 3600)) + 224 * cos((pi * t) / (6 * 3600)) + 210;
        q_in = q * SA_window;
        q_out = (T_storage_i - T_AIR) / R;
        dTdt = q_in/C- q_out/C;
        res = dTdt;
    end
end
