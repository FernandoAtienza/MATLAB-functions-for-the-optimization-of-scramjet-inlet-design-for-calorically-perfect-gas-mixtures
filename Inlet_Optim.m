% Compute Air properties at each scramjet inlet stage given Free-stream Mach number,
% altitude and number of ramps for calorically perfect gaseous mixtures.
%
% @author: Fernando Atienza González
%          Senior Undergraduate Aerospace Engineer Student
%          Universidad Carlos III de Madrid
%
% Last updated: 7th June 2024
% -------------------------------------------------------------------------
function [beta, theta, M, M3_id, Mn, p_ratio, rho_ratio, T_ratio, pt_ratio, T, Tt, T_id, p, rho, TPR, g_f, R_f] = Inlet_Optim(M1, n, h, f, self)
    
    % Error and step definition
    error = 0.0001;                       % Error of 0.01%
    step = error * 0.1;                     % Step of PR_guess update

    % Variable Preallocation:
    M = zeros(length(M1), n+2);             % Mach number
    M3_id = zeros(length(M1), 1);           % Ideal Mach number at the combsution chamber
    Mn = zeros(length(M1), n+1);            % Normal Mach number
    T = zeros(length(M1), n+2);             % Temperature
    Tt = zeros(length(M1), n+1);            % Total Temperature
    p = zeros(length(M1), n+2);             % Pressure
    beta = zeros(length(M1), n+1);          % Shock wave angle
    theta = zeros(length(M1), n+1);         % Flow defelction angle
    p_ratio = zeros(length(M1), n+1);       % Pressure ratio across shock wave
    pt_ratio = zeros(length(M1), n+1);      % Stagnation Pressure ratio across shock wave
    rho_ratio = zeros(length(M1), n+1);     % Density ratio across shock wave
    rho = zeros(length(M1), n+2);           % Density
    T_ratio = zeros(length(M1), n+1);       % Temperature ratio across shock wave
    TPR = zeros(length(M1), 1);             % Total Pressure Recovery (Ptf/Pt0)

    M(:, 1) = M1;                                   % Initial Mach number

    [rho(:,1),~,T(:,1),p(:,1),~,~,~] = atmos(h);    % Ambient Temperature [K] & Pressure [Pa] @ altitude h [m]

    [T_id, g, g_f, R_f] = Ideal_Temp (f, self);                     % Ideal Temperature required at the Combustion Chamber
                                                    
    PR_guess = 1;                                % External Stagnation Pressure Ratio Guess
    
    % Calculation of air properties at each stage:
    for j = 1:length(M1)
        
        % Ideal Mach number at the Combustion Chamber
        M3_id(j,1) = sqrt( (2/(g-1)) * (( (1 + (g-1)*0.5*((M(j,1))^2)) * (T(j,1)./T_id) - 1)) );

        a = 0;            % Variable to ensure that the Optimum External Pressure Ratio is found
        b = 0;            % Variable to avoid entering in an infite loop

        % Initialization: shock wave angle guess
        beta (j,1) = 15;

        while a==0
                
            % EXTERNAL Compressions:
            for i = 1:n 
                
                % Shock Wave Angle (Guessing the Optimal Pressure Ratio) w.r.t the flow
                beta (j,i) = newton_raphson(beta(j,i), beta_fun2(g, M(j,i), PR_guess), dbeta_fun2(g, M(j,i)));
                
                % Flow turning angle
                theta(j,i) = atand(2 * cotd(beta(j,i)) * ((M(j,i)^2 * (sind(beta(j,i)))^2) - 1) ./ ...
                    (((g + cosd(2 * beta(j,i))) * M(j,i)^2) + 2));

                % Mach number
                M(j,i+1) = sqrt((1 + (g-1)*0.5*((M(j,i) * sind(beta(j,i)))^2))./ ...
                    (g*((M(j,i) * sind(beta(j,i)))^2) - (g-1)*0.5))./(sind(beta(j,i) - theta (j,i)));
    
                % Normal Mach number
                Mn(j, i) = M(j, i) * sind(beta(j,i));
        
                % Pressure ratio
                p_ratio(j,i) = (2 * g * (M(j,i) * sind(beta(j,i)))^2 - (g-1) ) / (g+1);
    
                % Pressure
                p(j,i+1) = p (j,i) * p_ratio(j,i);
                
                % Density ratio 
                rho_ratio(j,i) = ( (g+1) * ((M(j,i) * sind(beta(j,i)))^2) ) / ( (g-1)* (M(j,i) * sind(beta(j,i)))^2 + 2);
                
                % Density
                rho(j,i+1) = rho_ratio(j,i) * rho(j,i);

                % Temperature ratio
                T_ratio(j,i) = p_ratio(j,i) * rho_ratio(j,i)^(-1);
    
                % Temperature
                T(j,i+1) = T (j,i) * T_ratio(j,i);
    
                % Total Temperature
                Tt(j,i) = T (j,i) * (1 + (g-1)*0.5 * (M(j,i)^2));
            
                % Total Pressure ratio: Double check with the PR given
                pt_ratio(j,i) = rho_ratio(j,i)^(g/(g-1)) * p_ratio(j,i)^(-1/(g-1));
                
                % Update guess
                beta (j,i + 1) = beta (j,i);
            end
            
            % INTERNAL compression Stage

            % Flow turning angle
            theta(j,end) = sum (theta(j,1:end-1));

            % Shock wave angle (Using the M-θ-β relationship)
            beta (j, end) = newton_raphson( 1.5 * beta (j, end), Mtb_fun2(g, M(j,end-1), theta(j, end)), dMtb_fun2(g, M(j,end-1)));

            % Mach number at the Combustion Chamber
            M(j,end) = sqrt((1 + (g-1)*0.5*((M(j,end-1) * sind(beta(j,end)))^2))./ ...
                    (g*((M(j,end-1) * sind(beta(j,end)))^2) - (g-1)*0.5))./(sind(beta(j,end) - theta (j,end)));

            % Normal Mach number
            Mn(j, end) = M(j, end-1) * sind(beta(j,end));
    
            % Pressure ratio
            p_ratio(j,end) = (2*g*((M(j,end-1) * sind(beta(j,end)))^2) - (g-1))/(g+1);
    
            % Pressure
            p(j,end) = p (j,end-1) * p_ratio(j,end);
            
            % Density ratio
            rho_ratio(j,end) = ((g+1)*((M(j,end-1) * sind(beta(j,end)))^2))/((g-1)*((M(j,end-1) * sind(beta(j,end)))^2) +2);
            
            % Density
            rho(j, end) = rho_ratio(j, end) * rho(j, end - 1);

            % Temperature ratio
            T_ratio(j,end) = p_ratio(j,end) * rho_ratio(j,end)^(-1);
    
            % Temperature
            T (j,end) = T_ratio(j,end) * T (j, end-1);
    
            % Total Temperature
            Tt(j,end) = T (j,end) * (1 + (g-1)*0.5 * (M(j,end)^2));
    
            % Total Pressure ratio:
            pt_ratio(j,end) = rho_ratio(j,end)^(g/(g-1)) * p_ratio(j,end)^(-1/(g-1));
          
            % Total Pressure Recovery (Ptf/Pt0) after ALL compressions:
            TPR (j) = prod(pt_ratio(j,:));              
                
            
            % Check results and PR_guess update
            if (abs((T(j,end)-T_id)/T_id) <= error)   % The Temperature obtained is within the error criteria (exits the loop)
                b = 0;
                a = 1;

            elseif T(j,end) < T_id                    % Pressure ratio needs to decrease
                PR_guess = PR_guess - step;
                b = b + 1;
            else                                      % Pressure ratio needs to increase
                PR_guess = PR_guess + step; 
                b = b + 1;
            end
            
            if b == 1/error && j>1
                    a = 1;
                    b = 0;
                    fprintf ('The Optimum Pressure ratio wasn´t found for M∞ = %.2f\n', M(j,1))
                    M(j,:) = NaN;             % Mach number
                    Mn(j,:) = NaN;            % Normal Mach number
                    T(j,:) = NaN;             % Temperature
                    Tt(j,:) = NaN;            % Total Temperature
                    p(j,:) = NaN;             % Pressure
                    beta(j,:) = NaN;          % Shock wave angle
                    theta(j,:) = NaN;         % Flow defelction angle
                    p_ratio(j,:) = NaN;       % Pressure ratio across shock wave
                    pt_ratio(j,:) = NaN;      % Stagnation Pressure ratio across shock wave
                    rho_ratio(j,:) = NaN;       % Density ratio across shock wave
                    T_ratio(j,:) = NaN;       % Temperature ratio across shock wave
                    TPR(j,:) = NaN;           % Total Pressure Recovery (Ptf/Pt0)
            end
            
             
        end % Optimum Condition Found
    end % Next M nº

end