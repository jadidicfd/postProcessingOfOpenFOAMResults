%% Header
% @author: M.Jadidi
% @program name: temporal auto correlation
% @dependency: OpenFoam data in the format of "probes" MUST be availabe in the working Dir. 
% @dependency: "mjReadData.m" MUST be MUST be run before this m-file.

% @task: Generating an temporal auto correlation from a time signal
%        the time signal can be diveded into "NT_segment" part for smooth
%        plot


% @see: [1]	M. Jadidi, H. K. Param, A. Revell, and Y. Mahmoudi, "Flow leakage and Kelvin–Helmholtz instability of turbulent flow over porous media," Physics of Fluids, vol. 34, no. 10, p. 105114, 2022/10/01 2022, doi: 10.1063/5.0111195.
% @see: [2]	M. Jadidi, H. K. Param, A. Revell, and Y. Mahmoudi, "Large eddy simulations of turbulent heat transfer in packed bed energy storage systems," Journal of Energy Storage, vol. 59, p. 106449, 2023/03/01/ 2023, doi: https://doi.org/10.1016/j.est.2022.106449.
% @see: [3]	M. Jadidi, A. Revell, and Y. Mahmoudi, "Pore-scale large eddy simulation of turbulent flow and heat transfer over porous media," Applied Thermal Engineering, vol. 215, p. 118916, 2022/10/01/ 2022, doi: https://doi.org/10.1016/j.applthermaleng.2022.118916.
% @see: https://uk.mathworks.com/matlabcentral/fileexchange/17204-kernel-density-estimation


% @created: May 2021,
% @version 01 (May 2021)
% @version 02 (June 2021)


%*******************IMPORTANT NOTE*****************************************
% The variable "Tmax" should be declared by the user. it should be less than
% "T" calculated by the "mjReadData.m"
%**************************************************************************

close all;

%% Given by user: should be compeleted based on your data:
Tmax = 0.2;        

%% Calculated values:

nTMax = int32(Tmax/del_t);
NT_segment = floor(r/nTMax);
dataTotal_T= nTMax*NT_segment;       %T stands for time not Temp.
ttt=[0:del_t:Tmax];            %sec

%% segmentation of original data for figure dispaly

uPrime_ijk_T = zeros(nTMax,c,NT_segment);       % c: Number of probes + 1 (columns of data_U), 
                                                 % Be careful: the first column is time
vPrime_ijk_T = zeros(nTMax,c,NT_segment);                                              
wPrime_ijk_T = zeros(nTMax,c,NT_segment);
TPrime_ijk_T = zeros(nTMax,c,NT_segment);

 k=1;
 z=1;
     for i=1:dataTotal_T
         uPrime_ijk_T(z,:,k)= uSigma(i,:);       %T stands for time not Temp.
         vPrime_ijk_T(z,:,k)= vPrime(i,:);       %T stands for time not Temp.
         wPrime_ijk_T(z,:,k)= wPrime(i,:);       %T stands for time not Temp.
         TPrime_ijk_T(z,:,k)= TPrime(i,:);       %T stands for time not Temp.
         
        
         z=z+1;
         if isequal (mod(i,nTMax),0)
             k=k+1;
             z=1;             
         end
     end


%%

uPrime_T_Rijk = zeros(nTMax,c,NT_segment);
vPrime_T_Rijk = zeros(nTMax,c,NT_segment);
wPrime_T_Rijk = zeros(nTMax,c,NT_segment);
TPrime_T_Rijk = zeros(nTMax,c,NT_segment);

%% this way of auto correlation calculation did not lead to accurate results! i dont know why!!!

% for k=1:nTMax
%     for j = 3:c
%         for i=1:nTMax
%             uPrime_Rijk(i,j,k)= (uPrime_ijk(1,j,k).*uPrime_ijk(i,j,k))./((uPrime_ijk(1,j,k).*uPrime_ijk(1,j,k)));
%             vPrime_Rijk(i,j,k)= (vPrime_ijk(1,j,k).*vPrime_ijk(i,j,k))./((vPrime_ijk(1,j,k).*vPrime_ijk(1,j,k)));
%             wPrime_Rijk(i,j,k)= (wPrime_ijk(1,j,k).*wPrime_ijk(i,j,k))./((wPrime_ijk(1,j,k).*wPrime_ijk(1,j,k)));
%             TPrime_Rijk(i,j,k)= (TPrime_ijk(1,j,k).*TPrime_ijk(i,j,k))./((TPrime_ijk(1,j,k).*TPrime_ijk(1,j,k)));
%         end
%     end
% end

%% the following method for auto correlation calculation is baed on "Lars Davidson" metho in his book
for k=1:NT_segment
    for j=2:c
            uPrime_T_Rijk(:,j,k)= autocorr(uPrime_ijk_T(:,j,k),nTMax-1);
            vPrime_T_Rijk(:,j,k)= autocorr(vPrime_ijk_T(:,j,k),nTMax-1);
            wPrime_T_Rijk(:,j,k)= autocorr(wPrime_ijk_T(:,j,k),nTMax-1);
            TPrime_T_Rijk(:,j,k)= autocorr(TPrime_ijk_T(:,j,k),nTMax-1);
    end
end

uPrime_T_RijMean = mean(uPrime_T_Rijk,3);     %T stands for time not Temp.
vPrime_T_RijMean = mean(vPrime_T_Rijk,3);     %T stands for time not Temp.
wPrime_T_RijMean = mean(wPrime_T_Rijk,3);   %T stands for time not Temp.
TPrime_T_RijMean = mean(TPrime_T_Rijk,3);   %the second T stands for time not Temp.


%% Figures for auto corrlation 
for j=2:c
    figure('Name',['--> probe: ',num2str(j-2)]);
    
    subplot(4,2,2)
    hold on;
    kisi=[0:del_t:(size(uPrime_T_RijMean, 1)-1)*del_t]';
    plot(kisi,uPrime_T_RijMean(:,j),'-r','LineWidth',2)
    grid on;
    xlabel('\bfkisi (sec)')
    ylabel('\bfTemporal R_u_u')
    title(['Temporal auto-correlations: uPrime, probe ',num2str(j-2)],'FontSize',12)
    hold off;
    
    subplot(4,2,4)
    hold on;
    plot(kisi,vPrime_T_RijMean(:,j),'-b','LineWidth',2)
    grid on;
    xlabel('\bfkisi (sec)')
    ylabel('\bfTemporal R_v_v')
    title(['Temporal auto-correlations of: vPrime, probe ',num2str(j-2)],'FontSize',12)
    hold off;
    
    subplot(4,2,6)
    hold on;
    plot(kisi,wPrime_T_RijMean(:,j),'-g','LineWidth',2)
    grid on;
    xlabel('\bfkisi (sec)')
    ylabel('\bfTemporal R_w_w')
    title(['Temporal auto-correlations: wPrime, probe ',num2str(j-2)],'FontSize',12)
    hold off;
    
    subplot(4,2,8)
    hold on;
    plot(kisi,TPrime_T_RijMean(:,j),'-y','LineWidth',2)
    grid on;
    xlabel('\bfkisi (sec)')
    ylabel('\bfTemporal R_T_T')
    title(['Temporal auto-correlations: TPrime, probe ',num2str(j-2)],'FontSize',12)
    hold off;
    
    %fluctuation figures
    k=1; %change k for fluctution of different segment, k=1 is first segment    
    
    subplot(4,2,1)
    hold on;
    tttt=[0:del_t:(size(uPrime_ijk_T(:,j,1), 1)-1)*del_t]';
    plot(tttt,uPrime_ijk_T(:,j,1),'-r','LineWidth',2)
    grid on;
    xlabel('\bft (sec)')
    ylabel('\bfuPrime')
    title(['fluctuation for probe: ',num2str(j-2)],'FontSize',12)
    hold off;
    
    subplot(4,2,3)
    hold on;
    tttt=[0:del_t:(size(vPrime_ijk_T(:,j,k), 1)-1)*del_t]';
    plot(tttt,vPrime_ijk_T(:,j,k),'-b','LineWidth',2)
    grid on;
    xlabel('\bft (sec)')
    ylabel('\bfvPrime')
    title(['fluctuation for probe: ',num2str(j-2)],'FontSize',12)
    hold off;
    
    subplot(4,2,5)
    hold on;
    tttt=[0:del_t:(size(wPrime_ijk_T(:,j,k), 1)-1)*del_t]';
    plot(tttt,wPrime_ijk_T(:,j,k),'-g','LineWidth',2)
    grid on;
    xlabel('\bft (sec)')
    ylabel('\bfwPrime')
    title(['fluctuation for probe: ',num2str(j-2)],'FontSize',12)
    hold off;
    
    subplot(4,2,7)
    hold on;
    tttt=[0:del_t:(size(TPrime_ijk_T(:,j,k), 1)-1)*del_t]';
    plot(tttt,TPrime_ijk_T(:,j,k),'-y','LineWidth',2)
    grid on;
    xlabel('\bft (sec)')
    ylabel('\bfTPrime')
    title(['fluctuation for probe: ',num2str(j-2)],'FontSize',12)
    hold off;
    
    set(gcf, 'Position',  [300, 50, 1000, 700])
end

%% Integral time scale, integralT,
   %The integral time scale represents the time scale of the large energy-containing eddies
   % You should carry out the integration only up to the ...
   % point when the autocorrelation goes negative (Lars Davidson)
 