%% Header
% @author: M.Jadidi

% @program name: OpenFoamToMatlabImport
% @dependency: OpenFoam data in the format of "probes in OpenFOAM" MUST be... 
% availabe in the working Dir. 
% @task: import, sort, calculate and dispay time signal from OpenFOAM


% @see: [1]	M. Jadidi, H. K. Param, A. Revell, and Y. Mahmoudi, "Flow leakage and Kelvin–Helmholtz instability of turbulent flow over porous media," Physics of Fluids, vol. 34, no. 10, p. 105114, 2022/10/01 2022, doi: 10.1063/5.0111195.
% @see: [2]	M. Jadidi, H. K. Param, A. Revell, and Y. Mahmoudi, "Large eddy simulations of turbulent heat transfer in packed bed energy storage systems," Journal of Energy Storage, vol. 59, p. 106449, 2023/03/01/ 2023, doi: https://doi.org/10.1016/j.est.2022.106449.
% @see: [3]	M. Jadidi, A. Revell, and Y. Mahmoudi, "Pore-scale large eddy simulation of turbulent flow and heat transfer over porous media," Applied Thermal Engineering, vol. 215, p. 118916, 2022/10/01/ 2022, doi: https://doi.org/10.1016/j.applthermaleng.2022.118916.
% @see: https://uk.mathworks.com/matlabcentral/fileexchange/17204-kernel-density-estimation


% @created: May 2021,
% @version 01 (May 2021)
% @version 02 (June 2021)

%*******************IMPORTANT NOTE*****************************************

% The variable "del_t" should be declared by the user
% The variable "delimiterIn" should be declared by the user
% The variable "headerlinesIn" should be declared by the user

%**************************************************************************

clear; close all;

%% Parameters given by user: should be compeleted based on your data:
del_t = 0.00005;  % time step size in your CFD (s)

%% reading data of probes from Openfoam to matlab

delimiterIn = ' ';   % data MUST be separeted form each other by " space"
headerlinesIn = 30;  % variable start to record after line 99. the lines 
                     %before 99 are location of probes
%%%For multiple file inputs
                 
myStructure = importdata('p1',delimiterIn,headerlinesIn);
data_p1 = myStructure.data;
myStructure = importdata('T1',delimiterIn,headerlinesIn);
data_T1 = myStructure.data;
myStructure = importdata('U1',delimiterIn,headerlinesIn);
data_U1 = myStructure.data;

myStructure = importdata('p2',delimiterIn,headerlinesIn);
data_p2 = myStructure.data;
myStructure = importdata('T2',delimiterIn,headerlinesIn);
data_T2 = myStructure.data;
myStructure = importdata('U2',delimiterIn,headerlinesIn);
data_U2 = myStructure.data;

data_p = [data_p1;data_p2];
data_T = [data_T1;data_T2];
data_U = [data_U1;data_U2];  

%%%For single file input

% myStructure = importdata('p',delimiterIn,headerlinesIn);
% data_p = myStructure.data;
% myStructure = importdata('T',delimiterIn,headerlinesIn);
% data_T = myStructure.data;
% myStructure = importdata('U',delimiterIn,headerlinesIn);
% data_U = myStructure.data;

% myStructure = importdata('nut',delimiterIn,headerlinesIn);
% data_nut = myStructure.data;

%% Parameters given by user: should be compeleted based on your data:
c = size(data_T,2);            %number of probes + 1 (columns of data_U),                             %Be careful: the first column is time
r = size (data_T, 1);        %number of time step (rows of data_U)
%% Calculated values:
T = del_t*r;    % total actual sample time (s)
f_s =1/del_t;   % sampling frequency (Hz)
t = [0:del_t:T-del_t]'; % time, t (s)
%% sperating data_U which is a vector to its components: (U_x, U_y, U_z)
Ux = zeros(r,c);
Uy = zeros(r,c);
Uz = zeros(r,c);

% first column of each variable (U_x, U_y, U_z, T, P, ...) should be time
Ux(:,1) = data_U(:,1);
Uy(:,1) = data_U(:,1);
Uz(:,1) = data_U(:,1);

m=2;            % the first column (m=1) is time
for j = 2:c
    Ux(:,j) = data_U(:,m);      %U_x components statrs at second column.
    Uy(:,j) = data_U(:,m+1);    %U_y components statrs at third column.
    Uz(:,j) = data_U(:,m+2);    %U_y components statrs at forth column.
    m=m+3;
end
%% Example: for testing with random number
% Ux=[randn(200,2); randn(200,1),randn(200,1);];
% Uy=[randn(200,2); randn(200,1),randn(200,1);];
% Uz=[randn(200,2); randn(200,1),randn(200,1);];
% data_T=[randn(200,2); randn(200,1),randn(200,1);];
% data_p=[randn(200,2); randn(200,1),randn(200,1);];
% 
% r=size (data_T,1);
% c=size (data_T,2);

%% calculating Mean value and fluctuations u' = u - <u>; T'= T-<T>

uMean = mean(Ux);
vMean = mean(Uy);
wMean = mean(Uz);
TMean = mean(data_T);
pMean = mean(data_p);

uPrime = zeros(r,c);
vPrime = zeros(r,c);
wPrime = zeros(r,c);
TPrime = zeros(r,c);
pPrime = zeros(r,c);

for j= 2:c       %j=1 is time column
    uPrime(:, j) = Ux(:, j) - uMean(j);
    vPrime(:, j) = Uy(:, j) - vMean(j);
    wPrime(:, j) = Uz(:, j) - wMean(j);
    TPrime(:, j) = data_T(:, j) - TMean(j);
    pPrime(:, j) = data_p(:, j) - pMean(j);
end

%% Calculating RMS and nondimensionalizing the fluctuations
   %(uSigma=u'/uRMS; vSigmav=v'/vRMS; TSigmaT=T'/TRMS; pSigmap=p'/pRMS)

uRMS = rms(uPrime);
vRMS = rms(vPrime);
wRMS = rms(wPrime);
TRMS = rms(TPrime);
pRMS = rms(pPrime);

uSigma = zeros(r,c);
vSigma = zeros(r,c);
wSigma = zeros(r,c);
TSigma = zeros(r,c);
pSigma = zeros(r,c);

for j= 2:c       %j=1 is time column
    uSigma(:, j) = uPrime(:, j) / uRMS(j);
    vSigma(:, j) = vPrime(:, j) / vRMS(j);
    wSigma(:, j) = wPrime(:, j) / wRMS(j);
    TSigma(:, j) = TPrime(:, j) / TRMS(j);
    pSigma(:, j) = pPrime(:, j) / pRMS(j);
end

%% figures 

% Velocity, Temp. and pressure
for j= 2:c       %j=1 is time column
    
    figure('Name',['--> probe: ',num2str(j-2)]);
    set(gcf, 'Position',  [50, 50, 1400, 700])
    
    subplot(5,3,1)
    hold on;
    plot(t,Ux(:,j),'-r','LineWidth',2)
    grid on;
    xlabel('\bft (sec)')
    ylabel('\bfUx')
    title(['Streamwise velocity for probe: ',num2str(j-2)],'FontSize',12)
    hold off;
    
    subplot(5,3,4)
    hold on;
    plot(t,Uy(:,j),'-b','LineWidth',2)
    grid on;
    xlabel('\bft (sec)')
    ylabel('\bfUy')
    title(['Vertical velocity  for probe: ',num2str(j-2)],'FontSize',12)
    hold off;
    
    subplot(5,3,7)
    hold on;
    plot(t,Uz(:,j),'-g','LineWidth',2)
    grid on;
    xlabel('\bft (sec)')
    ylabel('\bfUz')
    title(['Spanwise velocity for probe: ',num2str(j-2)],'FontSize',12)
    hold off;
    
    subplot(5,3,10)
    hold on;
    plot(t,data_T(:,j),'-y','LineWidth',2)
    grid on;
    xlabel('\bft (sec)')
    ylabel('\bfT')
    title(['Temperature for probe: ',num2str(j-2)],'FontSize',12)
    hold off;
    
    subplot(5,3,13)
    hold on;
    plot(t,data_p(:,j),'-k','LineWidth',2)
    grid on;
    xlabel('\bft (sec)')
    ylabel('\bfp')
    title(['Pressure for probe: ',num2str(j-2)],'FontSize',12)
    hold off;
    
  % Fluctuations
    subplot(5,3,2)
    hold on;
    plot(t,uPrime(:,j),'-r','LineWidth',2)
    grid on;
    xlabel('\bft (sec)')
    ylabel('\bfuPrime')
    title(['fluctuation for probe: ',num2str(j-2)],'FontSize',12)
    hold off;
    
    subplot(5,3,5)
    hold on;
    plot(t,vPrime(:,j),'-b','LineWidth',2)
    grid on;
    xlabel('\bft (sec)')
    ylabel('\bfvPrime')
    title(['fluctuation for probe: ',num2str(j-2)],'FontSize',12)
    hold off;
    
    subplot(5,3,8)
    hold on;
    plot(t,wPrime(:,j),'-g','LineWidth',2)
    grid on;
    xlabel('\bft (sec)')
    ylabel('\bfwPrime')
    title(['fluctuation for probe: ',num2str(j-2)],'FontSize',12)
    hold off;
    
    subplot(5,3,11)
    hold on;
    plot(t,TPrime(:,j),'-y','LineWidth',2)
    grid on;
    xlabel('\bft (sec)')
    ylabel('\bfTPrime')
    title(['fluctuation for probe: ',num2str(j-2)],'FontSize',12)
    hold off;
    
    subplot(5,3,14)
    hold on;
    plot(t,pPrime(:,j),'-k','LineWidth',2)
    grid on;
    xlabel('\bft (sec)')
    ylabel('\bfpPrime')
    title(['fluctuation for probe: ',num2str(j-2)],'FontSize',12)
    hold off;
    
    % Nondimensional fluctuations
    subplot(5,3,3)
    hold on;
    plot(t,uSigma(:,j),'-r','LineWidth',2)
    grid on;
    xlabel('\bft (sec)')
    ylabel('\bfuSigma')
    title(['nondimentsional fluctuations for probe: ',num2str(j-2)],'FontSize',12)
    hold off;
    
    subplot(5,3,6)
    hold on;
    plot(t,vSigma(:,j),'-b','LineWidth',2)
    grid on;
    xlabel('\bft (sec)')
    ylabel('\bfvSigma')
    title([' nondimentsional fluctuations for probe: ',num2str(j-2)],'FontSize',12)
    hold off;
    
    subplot(5,3,9)
    hold on;
    plot(t,wSigma(:,j),'-g','LineWidth',2)
    grid on;
    xlabel('\bft (sec)')
    ylabel('\bfwSigma')
    title(['nondimentsional fluctuations for probe: ',num2str(j-2)],'FontSize',12)
    hold off;
    
    subplot(5,3,12)
    hold on;
    plot(t,TSigma(:,j),'-y','LineWidth',2)
    grid on;
    xlabel('\bft (sec)')
    ylabel('\bfTSigma')
    title(['nondimentsional fluctuations for probe: ',num2str(j-2)],'FontSize',12)
    hold off;
    
    subplot(5,3,15)
    hold on;
    plot(t,pSigma(:,j),'-k','LineWidth',2)
    grid on;
    xlabel('\bft (sec)')
    ylabel('\bfpSigma')
    title(['nondimentsional fluctuations for probe: ',num2str(j-2)],'FontSize',12)
    hold off;
    
end