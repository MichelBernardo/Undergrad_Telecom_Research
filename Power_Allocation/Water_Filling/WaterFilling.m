clear all;
close all;
clc;

%% Setting Data

vtChannel = 10*rand(1,10); %Channel Modeling
dPmax = 200; %Budget Power
dNoise_Power = 1; %Noise Power per Sub-Carrier
dBand_Carrier = 1; %Sub-Carrier Bandwidth
iS = 10; %discrete number of Sub-Carriers
iM = 11; %number of MCS avaliable

%% 1) Algorithm 1: The Hughes-Hartogs Bit-Loading Algorithm

    %Calculation the necessary power for each MCS per carrier
    vtSNR_min = [0; 1; 3; 7; 15; 31; 63; 127; 255; 511; 1023]; %Minimum SNR for each MCS
    mtPot_mn = dNoise_Power.*vtSNR_min./(vtChannel.^2); %isolating Pn in expression (3.1)
    
    %Calculation the incremental power matrix
    mtDPot_mn = zeros(size(mtPot_mn)); %copy the mtPot_mn to the inicial incremental power matrix's status
    
    for i = 1:iS %Cycles through each 
        for j = 2:iM %Cycles through each Sub Carrier
            mtDPot_mn(j,i) = mtPot_mn(j,i) - mtPot_mn((j-1),i); %incremental power matrix
        end
    end

    dP_used_HH = 0; %Total Power used by HH
    vtMCS_n_HH = zeros(1, iS); %MCS used in each Sub Carrier by HH
    vtPn_HH = zeros(1, iS); %Power in each Sub Carrier by HH
    
    while dP_used_HH < dPmax & sum(vtMCS_n_HH) ~= iS*iM
        dP_used_HH_last = dP_used_HH; %save the last value of Total Power Used
        vtMCS_n_HH_last = vtMCS_n_HH; %save the last vtMCS_n_HH's status
        vtPn_HH_last = vtPn_HH; %save the last vtPn_HH's status

        [dDP_min, n] = min(mtDPot_mn(1,:)); %step 9
        vtMCS_n_HH(n) = vtMCS_n_HH(n) + 1; %step 10
        vtPn_HH(n) = mtPot_mn(vtMCS_n_HH(n), n); %step 11
        dP_used_HH = dP_used_HH + dDP_min; %step 12

        if dP_used_HH > dPmax %if the new dP_used_HH's value is bigger than dPmax available:
            dP_used_HH = dP_used_HH_last; %the dP_used_HH'value be back to his last value
            vtMCS_n_HH = vtMCS_n_HH_last; %the vtMCS_n_HH's value be back to his last value
            vtPn_HH = vtPn_HH_last; %the vtPn_HH's value be back to his last value
            break %And The Loop While is interrupted
        end

        %step 13
        vtAux = mtDPot_mn(:,n); %Define vtAux as the n-column vector of mtDPot_mn
        vtAux(1) = []; %delete the first position of vtAux
        vtAux = [vtAux ; inf]; %Put inf as the last position of vtAux
        mtDPot_mn(:,n) = vtAux; %Update the n-column vector to vtAux
    end

    vtRate_n_HH = Rate_function_MCS(vtMCS_n_HH);
    dRate_HH = sum(vtRate_n_HH);

%% 2) Equal Power Allocation Result
    dSub_Carrier_Power = dPmax/iS; %assigns equal power value to each sub-carrier
    vtSNR_EPA = dSub_Carrier_Power.*(vtChannel.^2)./dNoise_Power; %applying the expression (3.1)
    vtRate_n_EPA = Rate_function_SNR(vtSNR_EPA); %Sets the rate for the SNR level for each sub-carrier
    dRate_Epa = sum(vtRate_n_EPA); %Calculates the Total Rate Transmitte
    vtMCS_n_EPA = MCS_function_Rate(vtRate_n_EPA); %Sets the vector of MCS used for each Sub-Carrier

%% 3) Optimal Solution by intlinprog
    %Setting data
    mtPot_nm = mtPot_mn';
    vtR_m = [0 1 2 3 4 5 6 7 8 9 10];
    
    %Setting intlinprog parameters        
        % Define the index vector (integer constraints)
            intcon = [1:(iS*iM)]; %Sets all decision variables as integer
        
        % Linear inequality constraints
            % first Inequality constraints matrix
            mtA_1 = zeros(iS,(iS*iM));
            for i = 0:(iS - 1)
                mtAaux = zeros(1,(iS*iM));
                for j = 1:iM
                    mtAaux(1,:) = [zeros(1,iM*i) ones(1,iM) zeros(1,iM*(iS - 1 - i))];
                end
                mtA_1((i + 1),:) = mtAaux;
            end
    
            %Second inequality constraints matrix
            mtA_2 = [];
            for i = 1:iS
                mtA_2 = [mtA_2 mtPot_nm(i,:)];
            end
            
            %Inequality constraints matrix A and b
            mtA = [mtA_1; mtA_2];
            vtb = [ones(iS,1); dPmax];
        
        % Equality constraints
            Aeq = [];
            beq = [];
        
        % Low and Upper bounds
            vtlb = zeros(1,iS*iM);
            vtub = ones(1,iS*iM);

        %objective function for the solver
            vtf = objective(vtR_m, iS);

        % Calls the solver
            [vtAll_Decision_Variables, dRate_Optimal] = intlinprog(vtf,intcon,mtA,vtb,[],[],vtlb,vtub);
            dRate_Optimal = -dRate_Optimal;

    %Sets the MCS vector by intlinprog
        %Sets All_Decision_Variables as a matrix M-by-S
        mtMCS_Optimal = zeros(iM,iS);
        vtMCS_Optimal = zeros(1,iS);
        for i = 1:iS
                mtMCS_Optimal(:,i) = vtAll_Decision_Variables(((i*iM)-iS):(i*iM)); %Sets All_Decision_Variables as a matrix M-by-S
                [M,I] = max(mtMCS_Optimal(:,i)); %assigns the index of the bit "1" to I
                vtMCS_Optimal(i) = I; %assigns I to the i-position of vtMCS_Optimal
        end

%% 4) Comparison

    % Setting Data
    dPmax = 200; %Budget Power
    dNoise_Power = 1; %Noise Power per Sub-Carrier
    dBand_Carrier = 1; %Sub-Carrier Bandwidth
    iS = 10; %discrete number of Sub-Carriers
    iM = 11; %number of MCS avaliable

    T = 1000; %Total de repetições da comparação

    vtHH = zeros(1,T);
    vtEPA = zeros(1,T);
    vtOptimal = zeros(1,T);

    for k = 1:T
        vtChannel = 10*rand(1,10); %Channel Modeling

        % 1) Algorithm 1: The Hughes-Hartogs Bit-Loading Algorithm ---------------------------------------------------------------------------------------------------
            %Calculation the necessary power for each MCS per carrier
            vtSNR_min = [0; 1; 3; 7; 15; 31; 63; 127; 255; 511; 1023]; %Minimum SNR for each MCS
            mtPot_mn = dNoise_Power.*vtSNR_min./(vtChannel.^2); %isolating Pn in expression (3.1)
            
            %Calculation the incremental power matrix
            mtDPot_mn = mtPot_mn; %copy the mtPot_mn to the inicial incremental power matrix's status
            
            for i = 1:iS %Cycles through each MCS
                for j = 2:iM %Cycles through each Sub Carrier
                    mtDPot_mn(j,i) = mtPot_mn(j,i) - mtPot_mn((j-1),i); %incremental power matrix
                end
            end
            
            dP_used_HH = 0; %Total Power used by HH
            vtMCS_n_HH = zeros(1, iS); %MCS used in each Sub Carrier by HH
            vtPn_HH = zeros(1, iS); %Power in each Sub Carrier by HH
            
            while dP_used_HH < dPmax & sum(vtMCS_n_HH) ~= iS*iM
                dP_used_HH_last = dP_used_HH; %save the last value of Total Power Used
                vtMCS_n_HH_last = vtMCS_n_HH; %save the last vtMCS_n_HH's status
                vtPn_HH_last = vtPn_HH; %save the last vtPn_HH's status
        
                [dDP_min, n] = min(mtDPot_mn(1,:)); %step 9
                vtMCS_n_HH(n) = vtMCS_n_HH(n) + 1; %step 10
                vtPn_HH(n) = mtPot_mn(vtMCS_n_HH(n), n); %step 11
                dP_used_HH = dP_used_HH + dDP_min; %step 12
        
                if dP_used_HH > dPmax %if the new dP_used_HH's value is bigger than dPmax available:
                    dP_used_HH = dP_used_HH_last; %the dP_used_HH'value be back to his last value
                    vtMCS_n_HH = vtMCS_n_HH_last; %the vtMCS_n_HH's value be back to his last value
                    vtPn_HH = vtPn_HH_last; %the vtPn_HH's value be back to his last value
                    break %And The Loop While is interrupted
                end
        
                %step 13
                vtAux = mtDPot_mn(:,n); %Define vtAux as the n-column vector of mtDPot_mn
                vtAux(1) = []; %delete the first position of vtAux
                vtAux = [vtAux ; inf]; %Put inf as the last position of vtAux
                mtDPot_mn(:,n) = vtAux; %Update the n-column vector to vtAux
            end
        
            vtRate_n_HH = Rate_function_MCS(vtMCS_n_HH);
            dRate_HH = sum(vtRate_n_HH);

        % 2) Equal Power Allocation Result ----------------------------------------------------------------------------------------------------------------------------
            dSub_Carrier_Power = dPmax/iS; %assigns equal power value to each sub-carrier
            vtSNR_EPA = dSub_Carrier_Power.*(vtChannel.^2)./dNoise_Power; %applying the expression (3.1)
            vtRate_n_EPA = Rate_function_SNR(vtSNR_EPA); %Sets the rate for the SNR level for each sub-carrier
            dRate_Epa = sum(vtRate_n_EPA); %Calculates the Total Rate Transmitte
            vtMCS_n_EPA = MCS_function_Rate(vtRate_n_EPA); %Sets the vector of MCS used for each Sub-Carrier

        % 3) Optimal Solution by intlinprog --------------------------------------------------------------------------------------------------------------------------
            %Setting data
            mtPot_nm = mtPot_mn';
            vtR_m = [0 1 2 3 4 5 6 7 8 9 10];
            
            %Setting intlinprog parameters        
                % Define the index vector (integer constraints)
                    intcon = [1:(iS*iM)]; %Sets all decision variables as integer
                
                % Linear inequality constraints
                    % first Inequality constraints matrix
                    mtA_1 = zeros(iS,(iS*iM));
                    for i = 0:(iS - 1)
                        mtAaux = zeros(1,(iS*iM));
                        for j = 1:iM
                            mtAaux(1,:) = [zeros(1,iM*i) ones(1,iM) zeros(1,iM*(iS - 1 - i))];
                        end
                        mtA_1((i + 1),:) = mtAaux;
                    end
            
                    %Second inequality constraints matrix
                    mtA_2 = [];
                    for i = 1:iS
                        mtA_2 = [mtA_2 mtPot_nm(i,:)];
                    end
                    
                    %Inequality constraints matrix A and b
                    mtA = [mtA_1; mtA_2];
                    vtb = [ones(iS,1); dPmax];
                
                % Equality constraints
                    Aeq = [];
                    beq = [];
                
                % Low and Upper bounds
                    vtlb = zeros(1,iS*iM);
                    vtub = ones(1,iS*iM);
        
                %objective function for the solver
                    vtf = objective(vtR_m, iS);
        
                % Calls the solver
                    [vtAll_Decision_Variables, dRate_Optimal] = intlinprog(vtf,intcon,mtA,vtb,[],[],vtlb,vtub);
                    dRate_Optimal = -dRate_Optimal;

        % Generating the Rate Vectors -------------------------------------------------------------------------------------------------------------------------------- 
            vtHH(1,k) = dRate_HH;
            vtEPA(1,k) = dRate_Epa;
            vtOptimal(1,k) = dRate_Optimal;

    end

    figure(1)
    hold on;
    cdfplot(lin2db(vtHH(1,:)));
    cdfplot(lin2db(vtEPA(1,:)));
    cdfplot(lin2db(vtOptimal(1,:)));
    legend('Hughes-Hartogs','Equal Power Allocation','Optimal Solution');

%% Functions
    % 1) Rate by SNR
    function [vtRate] = Rate_function_SNR(vtSNR)
        vtRate = zeros(1,length(vtSNR));
        for k = 1:length(vtSNR)
            if vtSNR(k)>=0 & vtSNR(k)<1
                vtRate(k) = 0;
            elseif vtSNR(k)>=1 & vtSNR(k)<3
                vtRate(k) = 1;
            elseif vtSNR(k)>=3 & vtSNR(k)<7
                vtRate(k) = 2;
            elseif vtSNR(k)>=7 & vtSNR(k)<15
                vtRate(k) = 3;
            elseif vtSNR(k)>=15 & vtSNR(k)<31
                vtRate(k) = 4;
            elseif vtSNR(k)>=31 & vtSNR(k)<63
                vtRate(k) = 5;
            elseif vtSNR(k)>=63 & vtSNR(k)<127
                vtRate(k) = 6;
            elseif vtSNR(k)>=127 & vtSNR(k)<255
                vtRate(k) = 7;
            elseif vtSNR(k)>=255 & vtSNR(k)<511
                vtRate(k) = 8;
            elseif vtSNR(k)>=511 & vtSNR(k)<1023
                vtRate(k) = 9;
            else
                vtRate(k) = 10;
            end
        end
    end
    
    % 2) MCS by Rate
    function [vtMCS] = MCS_function_Rate(vtRate)
        vtMCS = zeros(1,length(vtRate));
        for k = 1:length(vtRate)
            if vtRate(k) == 0
                vtMCS(k) = 1;
            elseif vtRate(k) == 1
                vtMCS(k) = 2;
            elseif vtRate(k) == 2
                vtMCS(k) = 3;
            elseif vtRate(k) == 3
                vtMCS(k) = 4;
            elseif vtRate(k) == 4
                vtMCS(k) = 5;
            elseif vtRate(k) == 5
                vtMCS(k) = 6;
            elseif vtRate(k) == 6
                vtMCS(k) = 7;
            elseif vtRate(k) == 7
                vtMCS(k) = 8;
            elseif vtRate(k) == 8
                vtMCS(k) = 9;
            elseif vtRate(k) == 9
                vtMCS(k) = 10;
            else
                vtMCS(k) = 11;
            end
        end
    end

    % 3) Rate by MCS
    function [vtRate] = Rate_function_MCS(vtMCS)
        vtRate = zeros(1,length(vtMCS));
        for k = 1:length(vtMCS)
            if vtMCS(k) == 1
                vtRate(k) = 0;
            elseif vtMCS(k) == 2
                vtRate(k) = 1;
            elseif vtMCS(k) == 3
                vtRate(k) = 2;
            elseif vtMCS(k) == 4
                vtRate(k) = 3;
            elseif vtMCS(k) == 5
                vtRate(k) = 4;
            elseif vtMCS(k) == 6
                vtRate(k) = 5;
            elseif vtMCS(k) == 7
                vtRate(k) = 6;
            elseif vtMCS(k) == 8
                vtRate(k) = 7;
            elseif vtMCS(k) == 9
                vtRate(k) = 8;
            elseif vtMCS(k) == 10
                vtRate(k) = 9;
            else
                vtRate(k) = 10;
            end
        end
    end

    % 4) Optimal Solution Expression
    function [All_Decision_Variables] = optimal_solution_function(r_m, S)
    r_m_aux = [];
    for i = 1:S
        r_m_aux = [r_m_aux r_m];
    end
        All_Decision_Variables = r_m_aux;
    end
    
    % 5) Objective Function
    function [obj] = objective(r_m, S)
        obj = -optimal_solution_function(r_m, S);
    end

    % 6) Converting from linear scale to dB
    function g = lin2db(x)
        g = 10*log10(x);
    end
