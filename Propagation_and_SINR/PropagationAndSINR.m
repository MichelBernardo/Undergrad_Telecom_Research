clear;
clc;

% Activite 01
%% Setting Data
iNumBS = 4; %Number of Base Stations
iNumTM = 4; %Number of Mobile Terminals per Base Station
dRadius = 500; %BS's Coverage Radius
dTheta = 2*pi/3; %Angle of positions

vtBSPos = [0 dRadius*exp(-j*[0 dTheta -dTheta])]; %Vector of BS's positions
mtBSPos = [0 0
           dRadius 0
           dRadius*cos(dTheta) dRadius*sin(dTheta)
           dRadius*cos(-dTheta) dRadius*sin(-dTheta)]; %Matriz of BS's positions
%% Question 1

%Setting Circles
    n = 500; %Number of points
    t = linspace(0,2*pi,n); %Points of plotting
    seno = dRadius*sin(t); %Calculating the sine.
    coss = dRadius*cos(t); %Calculating the cosine.
    
    %Parameterization of each coverage area
    x0 = mtBSPos(1,1) + seno;
    y0 = mtBSPos(1,2) + coss;
    x1 = mtBSPos(2,1) + seno;
    y1 = mtBSPos(2,2) + coss;
    x2 = mtBSPos(3,1) + seno;
    y2 = mtBSPos(3,2) + coss;
    x3 = mtBSPos(4,1) + seno;
    y3 = mtBSPos(4,2) + coss;
    
%Plotting The BS' positions and their coverage area
plot(mtBSPos(1,1), mtBSPos(1,2), 'b^',...
     mtBSPos(2,1), mtBSPos(2,2), 'g^',...
     mtBSPos(3,1), mtBSPos(3,2), 'r^',...
     mtBSPos(4,1), mtBSPos(4,2), 'k^',...
     x0, y0,'b.',...
     x1, y1,'g.',...
     x2, y2,'r.',...
     x3, y3,'k.');

 axis equal
 
 %% Question 2
 
 %Generation of all TM coordinates with reference to origin
 for i = 1:iNumBS
     mtTM = []; %Coordinates of TM j
     for j = 1:iNumTM
         vtRef_00 = dRadius*(2*rand(1,2)-1); %Random coordinates
         dDist = norm(vtRef_00 - [0 0]); %Distance between the generated point to the origin
         while dDist > dRadius %checks if the generated point is inside the coverage area
             vtRef_00 = dRadius*(2*rand(1,2)-1); %generates new coordinates
             dDist = norm(vtRef_00 - [0 0]); %checks the new distance
         end
         mtTM = [mtTM ; vtRef_00]; %TM's coordinates for BS i
     end
     mtCooTM_00(:,:,i) = mtTM; %Associates the coordinates with BS.
     mtCooTM(:,:,i) = mtCooTM_00(:,:,i) + mtBSPos(i,:); %TM's coordinates with respect to the respective BS.
 end
 
 %plotting
    %BS Positions
    figure(1)
    plot(mtBSPos(1,1), mtBSPos(1,2), 'b^',...   %BS positions
         mtBSPos(2,1), mtBSPos(2,2), 'g^',...
         mtBSPos(3,1), mtBSPos(3,2), 'r^',...
         mtBSPos(4,1), mtBSPos(4,2), 'k^',...
         x0, y0,'b.',...                        %BS Coverage area
         x1, y1,'g.',...
         x2, y2,'r.',...
         x3, y3,'k.',...
         mtCooTM(1:iNumTM,1,1),mtCooTM(1:iNumTM,2,1),'b*',...   %TM positions
         mtCooTM(1:iNumTM,1,2),mtCooTM(1:iNumTM,2,2),'g*',...
         mtCooTM(1:iNumTM,1,3),mtCooTM(1:iNumTM,2,3),'r*',...
         mtCooTM(1:iNumTM,1,4),mtCooTM(1:iNumTM,2,4),'k*');
 axis equal
 
 %% Question 3
 % (a)Mean Path Loss ------------------------------------------------------
 %Transposing the TM Coordinate Matrix
 for i = 1:iNumBS
     mtCooTM_t(:,:,i) = mtCooTM(:,:,i)';
 end
 
 iCount = 1; %Counter of TM
 
 %Distance Matrix
 for i = 1:iNumBS %cycles through each BS
     for j = mtCooTM_t(:,:,i) %provides TM's coordinates (in colunm vector) as arguments
         mtDist = [];
         for ii = 1:iNumBS
             jj = j'; %Transposing the colunm vector to row vector
             dDist = norm(jj-mtBSPos(ii,:)); %Distance between the TM j and the BS i
             mtDist = [mtDist dDist]; %Matrix of distance between the TM j to the all BS.
         end
         mtDistTM(iCount,:) = mtDist; %Associates the mtDist to the respective TM
         iCount = iCount + 1; %Next TM
     end
 end
 
 iCount = iCount - 1; % Total Number of TM
 
 for j = 1:iCount %cycles through each TM in order
     for i = 1:iNumBS
        mtPLoss_dB(j,i) = PathLoss(mtDistTM(j,i)); %The Main Path Loss for the TM j (in dB)
        mtPloss_li(j,i) = db2lin(mtPLoss_dB(j,i)); %The Main Path Loss for the TM j (in linear scale)
        mtPGain_li(j,i) = (mtPloss_li(j,i)).^(-1); %The Travel Gain for the TM j (in linear scale)
        mtPGain_dB(j,i) = lin2db(mtPGain_li(j,i)); %The Travel Gain for the TM j (in dB)
     end
     
 end 
 
 % (b)Shadowing -----------------------------------------------------------
 
 iSigma_b = 8; %standard deviation in dB
 
 for j = 1:iCount %cycles through each TM in order
    mtS_dB = [];
    mtS_li = [];
    for i = 1:iNumBS %cycles through each BS in order
        vtSGain_dB = iSigma_b*randn; %random variable: sigma = 8dB e mu = 0
        vtSGain_li = 10^(vtSGain_dB/10); %Gain associated with shading for TM j
        mtS_dB = [mtS_dB vtSGain_dB]; %Matrix of shading (in dB) to the TM j due all BS
        mtS_li = [mtS_li vtSGain_li]; %Matrix of shading (in linear scale) to the TM j due all BS
    end

    mtSGain_dB(j,:) = mtS_dB; %Associates the mtS_dB to the TM j
    mtSGain_li(j,:) = mtS_li; %Associates the mtS_li to the TM j
 end
 
 % (c)Fast Fading --------------------------------------------------------------
 
 iSigma_c = 1/sqrt(2); %standard deviation
 
 for j = 1:iCount %cycles through each TM
     mtRay_li = [];
     mtRay_dB = [];
     for i = 1:iNumBS %cycles through each BS
         cH_ij = iSigma_c*randn + iSigma_c*randn*i; %Random variable with Rayleigh distribution for TM j for the BS i
         vtRayleigh_li = (abs(cH_ij))^2; %Gain associated with Fast Fading fo TM j due to BS i
         vtRayleigh_dB = 10*log10(vtRayleigh_li); %Gain associated with Fast Fading fo TM j due to BS i (in dB)
         mtRay_li = [mtRay_li vtRayleigh_li]; %Matrix of fast fading for the TM j with respect all the BS
         mtRay_dB = [mtRay_dB vtRayleigh_dB]; %Matrix of fast fading for the TM j with respect all the BS (in dB)
     end
     mtRayleigh_li(j,:) = mtRay_li; %Associates the mtRay_li to the TM j
     mtRayleigh_dB(j,:) = mtRay_dB; %Associates the mtRay_dB to the TM j
 end
 
 %% Question 4
 dPT = 43; %Transmitted Power from BS to TM [dBm]
 dpT = dbm2lin(dPT); %Transmitted Power from BS to TM [linear scale]
 
 for j = 1:iCount %cycles through each TM
     mtpR_aux = [];
     mtPR_aux = [];
     for i = 1:iNumBS %cycles through each BS
         vtpR = RPow(dpT, mtPGain_li(j,i), mtSGain_li(j,i), mtRayleigh_li(j,i)); %Received power from the BS i to the TM j
         vtPR = lin2dbm(vtpR); %Received Power in dBm
         mtpR_aux = [mtpR_aux vtpR]; %Received Power for the TM j from all BS
         mtPR_aux = [mtPR_aux vtPR]; %Received Power for the TM j from all BS (in dBm)
     end
     mtpR(j,:) = mtpR_aux; %Associates mtpR_aux to the TM j
     mtPR(j,:) = mtPR_aux; %Associates mtPR_aux to the TM j
 end
 
 %% Question 5
 
 dPN = -116; %Average Noise Power in dBm
 dpN = dbm2lin(dPN); %Average Noise Power in linear scale
 
 %Rearranging the received power data
 mtpR_t = mtpR';
 
 for i = 1:iNumBS %cycles through each BS
     vtAux = [];
     mtAux = [];
     for j = (1 + (i-1)*iNumTM):(i*iNumTM) %cycles through each TM which belongs to the BS i
         vtAux = mtpR_t(1:iNumBS,j); % vtAux <---- the received power from all BS to the TM j
         mtAux = [mtAux vtAux]; %Matrix of received power data for the TMs which belongs to the BS i
     end
     mtpR_5(:,:,i) = mtAux; %Associates mtAux to the BS i
 end
 
 %Calculation of SINR
 for i = 1:iNumBS %Cycles through each BS
     sum_pR = 0;
     mtAux = [];
     vtAux = [];
     for j = 1:iNumTM  %cycles through each TM
         vtAux = mtpR_5(1:iNumBS,j,i);
         dpR_ii = vtAux(i);
         vtAux(i) = [];
         sum_pR = sum(vtAux); %Sum of received power for the TM j for the BS i when i != j
         sinr_tm = sinr(dpR_ii, sum_pR, dpN); %The SINR
         mtAux = [mtAux sinr_tm]; %Matrix of SINR for each TM which belongs to the BS i
     end
     mtSINR(i,:) = mtAux; %Associates mtAux to the BS i
 end
 
 %% Question 6
 iI = 100; %Repetitions
 for k = 1:iI
        
     % Question 2 ---------------------------------------------------------
     %Generation of all TM coordinates with reference to origin
     for i = 1:iNumBS
         mtTM = []; %Coordinates of TM j
         for j = 1:iNumTM
             vtRef_00 = dRadius*(2*rand(1,2)-1); %Random coordinates
             dDist = norm(vtRef_00 - [0 0]); %Distance between the generated point to the origin
             while dDist > dRadius %checks if the generated point is inside the coverage area
                 vtRef_00 = dRadius*(2*rand(1,2)-1); %generates new coordinates
                 dDist = norm(vtRef_00 - [0 0]); %checks the new distance
             end
             mtTM = [mtTM ; vtRef_00]; %TM's coordinates for BS i
         end
         mtCooTM_00(:,:,i) = mtTM; %Associates the coordinates with BS.
         mtCooTM(:,:,i) = mtCooTM_00(:,:,i) + mtBSPos(i,:); %TM's coordinates with respect to the respective BS.
     end

     %plotting
        %BS Positions
%         plot(mtBSPos(1,1), mtBSPos(1,2), 'b^',...   %BS positions
%              mtBSPos(2,1), mtBSPos(2,2), 'g^',...
%              mtBSPos(3,1), mtBSPos(3,2), 'r^',...
%              mtBSPos(4,1), mtBSPos(4,2), 'k^',...
%              x0, y0,'b.',...                        %BS Coverage area
%              x1, y1,'g.',...
%              x2, y2,'r.',...
%              x3, y3,'k.',...
%              mtCooTM(1:iNumTM,1,1),mtCooTM(1:iNumTM,2,1),'b*',...   %TM positions
%              mtCooTM(1:iNumTM,1,2),mtCooTM(1:iNumTM,2,2),'g*',...
%              mtCooTM(1:iNumTM,1,3),mtCooTM(1:iNumTM,2,3),'r*',...
%              mtCooTM(1:iNumTM,1,4),mtCooTM(1:iNumTM,2,4),'k*');
%      axis equal

     % Question 3 ---------------------------------------------------------
     % (a)Mean Path Loss --------------------------------------------------
     %Transposing the TM Coordinate Matrix
     for i = 1:iNumBS
         mtCooTM_t(:,:,i) = mtCooTM(:,:,i)';
     end

     iCount = 1; %Counter of TM

     %Distance Matrix
     for i = 1:iNumBS %cycles through each BS
         for j = mtCooTM_t(:,:,i) %provides TM's coordinates (in colunm vector) as arguments
             mtDist = [];
             for ii = 1:iNumBS
                 jj = j'; %Transposing the colunm vector to row vector
                 dDist = norm(jj-mtBSPos(ii,:)); %Distance between the TM j and the BS i
                 mtDist = [mtDist dDist]; %Matrix of distance between the TM j to the all BS.
             end
             mtDistTM(iCount,:) = mtDist; %Associates the mtDist to the respective TM
             iCount = iCount + 1; %Next TM
         end
     end

     iCount = iCount - 1; % Total Number of TM

     for j = 1:iCount %cycles through each TM in order
         for i = 1:iNumBS
            mtPLoss_dB(j,i) = PathLoss(mtDistTM(j,i)); %The Main Path Loss for the TM j (in dB)
            mtPloss_li(j,i) = db2lin(mtPLoss_dB(j,i)); %The Main Path Loss for the TM j (in linear scale)
            mtPGain_li(j,i) = (mtPloss_li(j,i)).^(-1); %The Travel Gain for the TM j (in linear scale)
            mtPGain_dB(j,i) = lin2db(mtPGain_li(j,i)); %The Travel Gain for the TM j (in dB)
         end

     end 

     % (b)Shadowing -------------------------------------------------------

     iSigma_b = 8; %standard deviation in dB

     for j = 1:iCount %cycles through each TM in order
        mtS_dB = [];
        mtS_li = [];
        for i = 1:iNumBS %cycles through each BS in order
            vtSGain_dB = iSigma_b*randn; %random variable: sigma = 8dB e mu = 0
            vtSGain_li = 10^(vtSGain_dB/10); %Gain associated with shading for TM j
            mtS_dB = [mtS_dB vtSGain_dB]; %Matrix of shading (in dB) to the TM j due all BS
            mtS_li = [mtS_li vtSGain_li]; %Matrix of shading (in linear scale) to the TM j due all BS
        end

        mtSGain_dB(j,:) = mtS_dB; %Associates the mtS_dB to the TM j
        mtSGain_li(j,:) = mtS_li; %Associates the mtS_li to the TM j
     end

     % (c)Fast Fading -----------------------------------------------------

     iSigma_c = 1/sqrt(2); %standard deviation

     for j = 1:iCount %cycles through each TM
         mtRay_li = [];
         mtRay_dB = [];
         for i = 1:iNumBS %cycles through each BS
             cH_ij = iSigma_c*randn + iSigma_c*randn*i; %Random variable with Rayleigh distribution for TM j for the BS i
             vtRayleigh_li = (abs(cH_ij))^2; %Gain associated with Fast Fading fo TM j due to BS i
             vtRayleigh_dB = 10*log10(vtRayleigh_li); %Gain associated with Fast Fading fo TM j due to BS i (in dB)
             mtRay_li = [mtRay_li vtRayleigh_li]; %Matrix of fast fading for the TM j with respect all the BS
             mtRay_dB = [mtRay_dB vtRayleigh_dB]; %Matrix of fast fading for the TM j with respect all the BS (in dB)
         end
         mtRayleigh_li(j,:) = mtRay_li; %Associates the mtRay_li to the TM j
         mtRayleigh_dB(j,:) = mtRay_dB; %Associates the mtRay_dB to the TM j
     end

     % Question 4 ---------------------------------------------------------
     dPT = 43; %Transmitted Power from BS to TM [dBm]
     dpT = dbm2lin(dPT); %Transmitted Power from BS to TM [linear scale]

     for j = 1:iCount %cycles through each TM
         mtpR_aux = [];
         mtPR_aux = [];
         for i = 1:iNumBS %cycles through each BS
             vtpR = RPow(dpT, mtPGain_li(j,i), mtSGain_li(j,i), mtRayleigh_li(j,i)); %Received power from the BS i to the TM j
             vtPR = lin2dbm(vtpR); %Received Power in dBm
             mtpR_aux = [mtpR_aux vtpR]; %Received Power for the TM j from all BS
             mtPR_aux = [mtPR_aux vtPR]; %Received Power for the TM j from all BS (in dBm)
         end
         mtpR(j,:) = mtpR_aux; %Associates mtpR_aux to the TM j
         mtPR(j,:) = mtPR_aux; %Associates mtPR_aux to the TM j
     end

     % Question 5 ---------------------------------------------------------

     dPN = -116; %Average Noise Power in dBm
     dpN = dbm2lin(dPN); %Average Noise Power in linear scale

     %Rearranging the received power data
     mtpR_t = mtpR';

     for i = 1:iNumBS %cycles through each BS
         vtAux = [];
         mtAux = [];
         for j = (1 + (i-1)*iNumTM):(i*iNumTM) %cycles through each TM which belongs to the BS i
             vtAux = mtpR_t(1:iNumBS,j); % vtAux <---- the received power from all BS to the TM j
             mtAux = [mtAux vtAux]; %Matrix of received power data for the TMs which belongs to the BS i
         end
         mtpR_5(:,:,i) = mtAux; %Associates mtAux to the BS i
     end

     %Calculation of SINR
     for i = 1:iNumBS %Cycles through each BS
         sum_pR = 0;
         mtAux = [];
         vtAux = [];
         for j = 1:iNumTM  %cycles through each TM
             vtAux = mtpR_5(1:iNumBS,j,i);
             dpR_ii = vtAux(i);
             vtAux(i) = [];
             sum_pR = sum(vtAux); %Sum of received power for the TM j for the BS i when i != j
             sinr_tm = sinr(dpR_ii, sum_pR, dpN); %The SINR
             mtAux = [mtAux sinr_tm]; %Matrix of SINR for each TM which belongs to the BS i
         end
         mtSINR(i,:) = mtAux; %Associates mtAux to the BS i
     end

    %----------------------------------------------------------------------
    
    mtSINR_6(:,:,k) = mtSINR; %Associates mtSINR with the repetition k
 end
 
%Setting the SINR ploting

 for k = 1:iI
     for i = 1:iNumBS
        mtSINR_plot(i,k) = mtSINR_6(i,1,k);
     end
 end
 
  figure(2)
  hold on;
 cdfplot(lin2db(mtSINR_plot(1,:)));
 cdfplot(lin2db(mtSINR_plot(2,:)));
 cdfplot(lin2db(mtSINR_plot(3,:)));
 cdfplot(lin2db(mtSINR_plot(4,:)));
 
  %axis equal;
 
 %% Functions
% 1) The Mean Path Loss Function
function PL_ij = PathLoss(d)
    PL_ij = 128.1 + 36.7*log10(d);
end

% 2) Converting from linear scale to dB
function g = lin2db(x)
    g = 10*log10(x);
end

% 3) Convertion from dB to linaer scale
function y = db2lin(x)
    y = 10.^(x./10);
end

% 4) Convertion from linear scale to dBm
function [y] = lin2dbm(x)
    y = 10*log10(x./1e-3);
end

% 5) Convertion from dBm to linear scale
function [y] = dbm2lin(x)
    y = 10.^(x./10-3);
end

% 6) Received Power from BS to TM
function [dpR] = RPow(dpT, PGain_li, SGain_li, Rayleigh_li)
    dpR = dpT*PGain_li*SGain_li*Rayleigh_li;
end

% 7) SINR
function sinr_tm = sinr(dpR, sum_pR, dpN)
    sinr_tm = dpR/(sum_pR + dpN);
end