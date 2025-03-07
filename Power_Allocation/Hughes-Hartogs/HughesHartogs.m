clear all;
close all;
clc;

%% Setting Data

vtChannel = 10*rand(1,10); %Channel Modeling
dPmax = 200; %Budget Power
dNoise_Power = 1; %Noise Power per Sub-Carrier
dBand_Carrier = 1; %Sub-Carrier Bandwidth
iN = 10; %discrete number of Sub-Carriers

%% 1) Water-Pouring Result

    [vtPower, vtChannel_new] = WP_result(vtChannel, dNoise_Power, dPmax); %Generating the Power Vector by Water-Pouring Result
    vtSNR_WP = SNR_func(vtPower, vtChannel_new, dNoise_Power); %Generating the SNR Vector
    vtSub_Carrier_Capacity_WP = Sub_Carrier_Capacity_func(dBand_Carrier, vtSNR_WP); %Generating the Sub Carrier Capacity Vector
    dEntire_System_Capacity_WP = sum(vtSub_Carrier_Capacity_WP); %Calculation of Entire System Capacity

%% 2) Equal Power Allocation Result

    dPower_Sub_Carrier = dPmax/iN; %Calculation of Power Allocation
    vtSNR_EPA = SNR_func_2(dPower_Sub_Carrier, vtChannel, dNoise_Power); %Generating the SNR Vector
    vtSub_Carrier_Capacity_EPA = Sub_Carrier_Capacity_func(dBand_Carrier, vtSNR_EPA); %Generating the Sub Carrier Capacity Vector
    dEntire_System_Capacity_EPA = sum(vtSub_Carrier_Capacity_EPA); %Calculation of Entire System Capacity

%% 3) Optimal Solution by fmincon

    %Setting fmincon parameters
        % Initial solution required by the solver.
        x0 = ones(iN,1);

        % Models in matrix form the total power constraint.
        A = ones(1,iN);
        b = dPmax;

        % Define the upper/lower bounds
        lb = zeros(iN,1);
        ub = inf(iN,1);

        f = @(x)objective(x, vtChannel);

    % Calls the solver
    [vtPower_Optimal, fval, exitflag] = fmincon(f,x0,A,b,[],[],lb,ub);
    
    vtPower_Optimal = vtPower_Optimal'; %Transposing the Power Vector.
        
    vtSNR_OS = SNR_func(vtPower_Optimal, vtChannel, dNoise_Power); %Generating the SNR Vector
    vtSub_Carrier_Capacity_OS = Sub_Carrier_Capacity_func(dBand_Carrier, vtSNR_OS); %Generating the Sub Carrier Capacity Vector
    dEntire_System_Capacity_OS = sum(vtSub_Carrier_Capacity_OS); %Calculation of Entire System Capacity

 %% 4) Comparison
    for k = 1:50
        
        vtChannel = 10*rand(1,10); %Channel Modeling
        
        % 1) Water-Pouring Result -----------------------------------------------------------------------------------------------------------
            [vtPower, vtChannel_new] = WP_result(vtChannel, dNoise_Power, dPmax); %Generating the Power Vector by Water-Pouring Result
            vtSNR_WP = SNR_func(vtPower, vtChannel_new, dNoise_Power); %Generating the SNR Vector
            vtSub_Carrier_Capacity_WP = Sub_Carrier_Capacity_func(dBand_Carrier, vtSNR_WP); %Generating the Sub Carrier Capacity Vector
            dEntire_System_Capacity_WP = sum(vtSub_Carrier_Capacity_WP); %Calculation of Entire System Capacity
            
        % 2) Equal Power Allocation Result --------------------------------------------------------------------------------------------------
            dPower_Sub_Carrier = dPmax/iN; %Calculation of Power Allocation
            vtSNR_EPA = SNR_func_2(dPower_Sub_Carrier, vtChannel, dNoise_Power); %Generating the SNR Vector
            vtSub_Carrier_Capacity_EPA = Sub_Carrier_Capacity_func(dBand_Carrier, vtSNR_EPA); %Generating the Sub Carrier Capacity Vector
            dEntire_System_Capacity_EPA = sum(vtSub_Carrier_Capacity_EPA); %Calculation of Entire System Capacity
            
        % 3) Optimal Solution by fmincon ---------------------------------------------------------------------------------------------------
            %Setting fmincon parameters
                % Initial solution required by the solver.
                x0 = ones(iN,1);

                % Models in matrix form the total power constraint.
                A = ones(1,iN);
                b = dPmax;

                % Define the upper/lower bounds
                lb = zeros(iN,1);
                ub = inf(iN,1);

            % Calls the solver
            [vtPower_Optimal, fval, exitflag] = fmincon(@(x)objective(x, vtChannel),x0,A,b,[],[],lb,ub);

            vtPower_Optimal = vtPower_Optimal'; %Transposing the Power Vector.

            vtSNR_OS = SNR_func(vtPower_Optimal, vtChannel, dNoise_Power); %Generating the SNR Vector
            vtSub_Carrier_Capacity_OS = Sub_Carrier_Capacity_func(dBand_Carrier, vtSNR_OS); %Generating the Sub Carrier Capacity Vector
            dEntire_System_Capacity_OS = sum(vtSub_Carrier_Capacity_OS); %Calculation of Entire System Capacity
            
            %Generating the Capacity Vectors -----------------------------------------------------------------------------------------------
            vtWaterPouring(1,k) = dEntire_System_Capacity_WP;
            vtEqualAllocation(1,k) = dEntire_System_Capacity_EPA;
            vtOptimalSolution(1,k) = dEntire_System_Capacity_OS;
    end
    
    figure(1)
    hold on;
%     axis equal;
    cdfplot(lin2db(vtWaterPouring(1,:)));
    cdfplot(lin2db(vtEqualAllocation(1,:)));
    cdfplot(lin2db(vtOptimalSolution(1,:)));
    legend('Water-Pouring','Equal Allocation','Optimal Solution');

%% Functions
% 1) Water-Pouring Result
function [vtPn_opt, vtChannel_new] = WP_result(vtChannel, dNoise_Power, dPmax)
    vtPn_opt = [];
    vtH = [];
    
    %Calculo do termo constante da solução analítica
    for j = 1:length(vtChannel)
        dNoise_Gain = dNoise_Power./(vtChannel(1,j)^2);
        vtH = [vtH dNoise_Gain];
    end
    dSum_Noise_Gain = sum(vtH);
    
    %Calculo do vetor de potências para as condições iniciais do WP
    vtSub_Carrier_negative = [];
    
    for i = 1:length(vtChannel)
        dPn_opt_calc = (1/length(vtChannel))*(dSum_Noise_Gain + dPmax) - dNoise_Power./(vtChannel(1,i)^2);
        
        %se o valor da potência for negativo o índice da portadora é guardado
        if dPn_opt_calc < 0
            vtSub_Carrier_negative = [vtSub_Carrier_negative i];
        end
        
        vtPn_opt = [vtPn_opt dPn_opt_calc]; %Vetor para armazenamento das potências
    end
    
    %Verifica a existência de portadora com valor de potência negativo
    if length(vtSub_Carrier_negative) > 0
        while length(vtSub_Carrier_negative) > 0
            
            vtChannel_aux = [];
            
            %Remover os ganhos das portadoras descartadas
            for j = 1:length(vtChannel)
                if vtPn_opt(j) > 0
                    vtChannel_aux = [vtChannel_aux vtChannel(j)];
                end
            end
            
            vtChannel = vtChannel_aux;
            
            %O Water-Pouring é novamente aplicado
                vtPn_opt = [];
                vtH = [];

                %Cálculo do termo constante da solução analítica
                for j = 1:length(vtChannel)
                    dNoise_Gain = dNoise_Power./(vtChannel(1,j)^2);
                    vtH = [vtH dNoise_Gain];
                end
                dSum_Noise_Gain = sum(vtH);

                %Cálculo das Potências
                vtSub_Carrier_negative = [];

                for i = 1:length(vtChannel)
                    dPn_opt_calc = (1/length(vtChannel))*(dSum_Noise_Gain + dPmax) - dNoise_Power./(vtChannel(1,i)^2);
                    if dPn_opt_calc < 0
                        vtSub_Carrier_negative = [vtSub_Carrier_negative i];
                    end
                    vtPn_opt = [vtPn_opt dPn_opt_calc];
                end
        end
        
        vtChannel_new = vtChannel; %Update the Gain Vector
        
    else
        vtChannel_new = vtChannel; %Update the Gain Vector
    end
end

% 2) SNR Function
function [vtSNR] = SNR_func(vtPower, vtChannel, dNoise_Power)
    vtSNR = [];
    for i = 1:length(vtChannel)
        dAux = vtPower(1,i).*(vtChannel(1,i)^2)./dNoise_Power;
        vtSNR = [vtSNR dAux];
    end
end

% 3) SNR Function - Question 2
function [vtSNR_2] = SNR_func_2(dPower_Sub_Carrier, vtChannel, dNoise_Power)
    vtSNR_2 = [];
    for i = 1:length(vtChannel)
        dAux = dPower_Sub_Carrier.*(vtChannel(1,i)^2)./dNoise_Power;
        vtSNR_2 = [vtSNR_2 dAux];
    end
end

% 4) Capacity Function
function [vtSub_Carrier_Capacity] = Sub_Carrier_Capacity_func(dBand_Carrier, vtSNR)
    vtSub_Carrier_Capacity = [];
    for i = 1:length(vtSNR)
        vtAux = dBand_Carrier*log2(1 + vtSNR(1,i));
        vtSub_Carrier_Capacity = [vtSub_Carrier_Capacity vtAux];
    end
end

% 5) Converting from linear scale to dB
function g = lin2db(x)
    g = 10*log10(x);
end

% 6) Optimal Solution Expression
function [dEntire_System_Capacity] = optimal_solution_Capacity(vtPower, vtChannel)
    dNoise_Power = 1;
    dBand_Carrier = 1;
    
    vtCap = [];
    for i = 1:length(vtChannel)
        x = log2(1 + (vtPower(i).*(vtChannel(i).^2))./dNoise_Power);
        vtCap = [vtCap x];
    end
    dEntire_System_Capacity = dBand_Carrier*sum(vtCap);
end

% 7) Objective Function
function [obj] = objective(vtPower, vtChannel)
    obj = -optimal_solution_Capacity(vtPower, vtChannel);
end
