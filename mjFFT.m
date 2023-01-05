%% Header
% *@author: M.Jadidi*

% @program name: FFT 
% @dependency: mjReadingData.m
% @dependency: OpenFoam data in the format of "probes" MUST be availabe in the working Dir.

% @task: Generating an FFT from a time signal - including Hanning windows
%        the time signal can be diveded into "N_segment" part for smooth
%        plot

% @see: [1]	M. Jadidi, H. K. Param, A. Revell, and Y. Mahmoudi, "Flow leakage and Kelvin–Helmholtz instability of turbulent flow over porous media," Physics of Fluids, vol. 34, no. 10, p. 105114, 2022/10/01 2022, doi: 10.1063/5.0111195.
% @see: [2]	M. Jadidi, H. K. Param, A. Revell, and Y. Mahmoudi, "Large eddy simulations of turbulent heat transfer in packed bed energy storage systems," Journal of Energy Storage, vol. 59, p. 106449, 2023/03/01/ 2023, doi: https://doi.org/10.1016/j.est.2022.106449.
% @see: [3]	M. Jadidi, A. Revell, and Y. Mahmoudi, "Pore-scale large eddy simulation of turbulent flow and heat transfer over porous media," Applied Thermal Engineering, vol. 215, p. 118916, 2022/10/01/ 2022, doi: https://doi.org/10.1016/j.applthermaleng.2022.118916.
% @see: https://uk.mathworks.com/matlabcentral/fileexchange/17204-kernel-density-estimation
% @see: Spectral Analysis, John M. Cimbala, Penn State University, https://www.me.psu.edu/cimbala/me345/

% @created: May 2021
% @version 01(May 2021)
% @version 02(June 2021)

%*******************IMPORTANT NOTE*****************************************

% The variable "N_segment" should be declared by the user

%**************************************************************************

close all
       
 %% Given by user: should be compeleted based on your data:
N = r;       % N=total actual number of data points;  
             % r=number of time step (rows of data_U)

N_segment = 10;      % the full length signal is devided into N_segment for smooth plot generation
n_segment = 2*floor((r/N_segment)/2);       % n_segment is number of point(or data) in each segment
                                            % n_segment Round to nearest even integer
dataTotal=int32(N_segment*n_segment);

%% Calculated values:

TT = del_t*n_segment;               % total actual sample time for one seqment (s)
f_s =1/del_t;                       % sampling frequency (Hz)
del_f = 1/TT;                       % (Hz)
f_fold = f_s/2;                     % folding frequency = max frequency of FFT (Hz)
n_freq = n_segment/2;               % number of discrete frequencies

tt = [0:del_t:TT-del_t]';           % time, t (s)
frequency = [0:del_f:f_fold]';      % frequency (Hz)      

%% Example (for validation)
% clear
% close all
% 
% N_segment = 2;      % the full length signal is devided into N_segment for smooth plot generation
% n_segment = 2000;       % n_segment is number of point(or data) in each segment
%                         % n_segment Round to nearest even integer
% dataTotal=int32(N_segment*n_segment);
% 
% del_t = 0.001;
% T = del_t*n_segment;                % total actual sample time for one seqment (s)
% f_s =1/del_t;                       % sampling frequency (Hz)
% del_f = 1/T;                        % (Hz)
% f_fold = f_s/2;                     % folding frequency = max frequency of FFT (Hz)
% n_freq = n_segment/2;               % number of discrete frequencies
% 
% t = [0:del_t:T-del_t]';             % time, t (s)
% frequency = [0:del_f:f_fold]';      % frequency (Hz)   

% % Form a signal containing a 50 Hz sinusoid of amplitude 0.7 and
% % a 120 Hz sinusoid of amplitude 1.
% S = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t);
% % Corrupt the signal with zero-mean white noise with a variance of 4.
% g = S;
% DC = mean(g);   % DC = mean value of input signal (V) (average of all the useful data)
% g_uncoupled = g-DC; % uncoupled
% u_Hann = 0.5*(1-cos(2*pi*t/T)); % u_Hanning(t)
% g_Hann = g_uncoupled.*u_Hann; % g(t)*u_Hanning(t)
% G_Hann = fft(g_Hann); % G(omega) with Hanning window
% Magnitude_Hann = abs(G_Hann)*sqrt(8/3)/(n_segment/2); % |F|*sqrt(8/3)/(N/2)
% plot(frequency, Magnitude_Hann(1:length(frequency)))
%--------------------------------------------------------------------------

%% segmentation of original data

uPrime_ijk = zeros(n_segment,c,N_segment);       % c: Number of probes + 1 (columns of data_U), 
                                                 % Be careful: the first column is time
vPrime_ijk = zeros(n_segment,c,N_segment);                                              
wPrime_ijk = zeros(n_segment,c,N_segment);
TPrime_ijk = zeros(n_segment,c,N_segment);

uSigma_ijk = zeros(n_segment,c,N_segment);       % c: Number of probes + 1 (columns of data_U), 
                                                 % Be careful: the first column is time
vSigma_ijk = zeros(n_segment,c,N_segment);                                              
wSigma_ijk = zeros(n_segment,c,N_segment);
TSigma_ijk = zeros(n_segment,c,N_segment);

 k=1;
 z=1;
     for i=1:dataTotal
         uPrime_ijk(z,:,k)= uSigma(i,:);
         vPrime_ijk(z,:,k)= vPrime(i,:);
         wPrime_ijk(z,:,k)= wPrime(i,:);
         TPrime_ijk(z,:,k)= TPrime(i,:);
         
         uSigma_ijk(z,:,k)= uSigma(i,:);
         vSigma_ijk(z,:,k)= vSigma(i,:);
         wSigma_ijk(z,:,k)= wSigma(i,:);
         TSigma_ijk(z,:,k)= TSigma(i,:);
         
         z=z+1;
         if isequal (mod(i,n_segment),0)
             k=k+1;
             z=1;             
         end
     end
     
     uPrime_ijkVAR=var(uPrime_ijk);
     vPrime_ijkVAR=var(vPrime_ijk);
     wPrime_ijkVAR=var(wPrime_ijk);
     TPrime_ijkVAR=var(TPrime_ijk);
     %% Frequency Spectrum:
u_Hann = 0.5*(1-cos(2*pi*tt/TT)); % u_Hanning(t)

g_Hann_uPrime = zeros(n_segment,c,N_segment);
G_Hann_uPrime = zeros(n_segment,c,N_segment);
Magnitude_Hann_uPrime = zeros(n_segment,c,N_segment);

g_Hann_vPrime = zeros(n_segment,c,N_segment);
G_Hann_vPrime = zeros(n_segment,c,N_segment);
Magnitude_Hann_vPrime = zeros(n_segment,c,N_segment);

g_Hann_wPrime = zeros(n_segment,c,N_segment);
G_Hann_wPrime = zeros(n_segment,c,N_segment);
Magnitude_Hann_wPrime = zeros(n_segment,c,N_segment);

g_Hann_TPrime = zeros(n_segment,c,N_segment);
G_Hann_TPrime = zeros(n_segment,c,N_segment);
Magnitude_Hann_TPrime = zeros(n_segment,c,N_segment);

for k=1:N_segment
    for j=2:c
        g_Hann_uPrime(:,j,k) = uPrime_ijk(:,j,k).*u_Hann; % g(t)*u_Hanning(t)
        G_Hann_uPrime(:,j,k) = fft(g_Hann_uPrime(:,j,k)); % G(omega) with Hanning window
        Magnitude_Hann_uPrime(:,j,k) = abs(G_Hann_uPrime(:,j,k))*sqrt(8/3)/(n_segment/2); % |F|*sqrt(8/3)/(N/2)
        
        g_Hann_vPrime(:,j,k) = vPrime_ijk(:,j,k).*u_Hann; % g(t)*u_Hanning(t)
        G_Hann_vPrime(:,j,k) = fft(g_Hann_vPrime(:,j,k)); % G(omega) with Hanning window
        Magnitude_Hann_vPrime(:,j,k) = abs(G_Hann_vPrime(:,j,k))*sqrt(8/3)/(n_segment/2); % |F|*sqrt(8/3)/(N/2)
        
        g_Hann_wPrime(:,j,k) = wPrime_ijk(:,j,k).*u_Hann; % g(t)*u_Hanning(t)
        G_Hann_wPrime(:,j,k) = fft(g_Hann_wPrime(:,j,k)); % G(omega) with Hanning window
        Magnitude_Hann_wPrime(:,j,k) = abs(G_Hann_wPrime(:,j,k))*sqrt(8/3)/(n_segment/2); % |F|*sqrt(8/3)/(N/2)
        
        g_Hann_TPrime(:,j,k) = TPrime_ijk(:,j,k).*u_Hann; % g(t)*u_Hanning(t)
        G_Hann_TPrime(:,j,k) = fft(g_Hann_TPrime(:,j,k)); % G(omega) with Hanning window
        Magnitude_Hann_TPrime(:,j,k) = abs(G_Hann_TPrime(:,j,k))*sqrt(8/3)/(n_segment/2); % |F|*sqrt(8/3)/(N/2)
    end
end

Magnitude_Hann_uPrimeMean=mean(Magnitude_Hann_uPrime,3);
Magnitude_Hann_vPrimeMean=mean(Magnitude_Hann_vPrime,3);
Magnitude_Hann_wPrimeMean=mean(Magnitude_Hann_wPrime,3);
Magnitude_Hann_TPrimeMean=mean(Magnitude_Hann_TPrime,3);


%% 

for j = 2:c
    
    figure('Name',['--> probe: ',num2str(j-2)]);
    
    %loglog plot
    subplot(4,2,1)
    loglog(frequency,Magnitude_Hann_uPrimeMean(1:length(frequency),j),'-r','LineWidth',2)
    
    hold on;
    
    negative5over3=frequency.^(-5/3);   %-5/3 law for energy cascade
    loglog(frequency,negative5over3,'k.','MarkerSize',0.5)

    grid on;
    xlabel('\bffrequency, \itf \rm\bf(Hz)')
    ylabel('\bf|F|')
    title(['uPrime PSD with Hanning window, probe ',num2str(j-2)],'FontSize',12)
    hold off;
    
    subplot(4,2,3)
    loglog(frequency,Magnitude_Hann_vPrimeMean(1:length(frequency),j),'-b','LineWidth',2)
    hold on;
    
    negative5over3=frequency.^(-5/3);   %-5/3 law for energy cascade
    loglog(frequency,negative5over3,'k.','MarkerSize',0.5)
    
    grid on;
    xlabel('\bffrequency, \itf \rm\bf(Hz)')
    ylabel('\bf|F|')
    title(['vPrime PSD with Hanning window, probe ',num2str(j-2)],'FontSize',12)
    hold off;
    
    subplot(4,2,5)
    loglog(frequency,Magnitude_Hann_wPrimeMean(1:length(frequency),j),'-g','LineWidth',2)
    hold on;
    
    negative5over3=frequency.^(-5/3);   %-5/3 law for energy cascade
    loglog(frequency,negative5over3,'k.','MarkerSize',0.5)
    
    grid on;
    xlabel('\bffrequency, \itf \rm\bf(Hz)')
    ylabel('\bf|F|')
    title(['wPrime PSD with Hanning window, probe ',num2str(j-2)],'FontSize',12)
    hold off;
    
    subplot(4,2,7)
    loglog(frequency,Magnitude_Hann_TPrimeMean(1:length(frequency),j),'-y','LineWidth',2)
    hold on;
    
    negative5over3=frequency.^(-5/3);   %-5/3 law for energy cascade
    loglog(frequency,negative5over3,'k.','MarkerSize',0.5)
    
    grid on;
    xlabel('\bffrequency, \itf \rm\bf(Hz)')
    ylabel('\bf|F|')
    title(['TPrime PSD with Hanning window, probe',num2str(j-2)],'FontSize',12)
    hold off;
    
    %semilog plot 
    subplot(4,2,2) 
    semilogx(frequency,Magnitude_Hann_uPrimeMean(1:length(frequency),j),'-r','LineWidth',2) 
    hold on;
    set(gcf, 'Position',  [400, 50, 400, 700])
    grid on;
    xlabel('\bffrequency, \itf \rm\bf(Hz)')
    ylabel('\bf|F|')
    title(['uPrime PSD with Hanning window, probe ',num2str(j-2)],'FontSize',12)
    hold off;
    
    subplot(4,2,4)
    semilogx(frequency,Magnitude_Hann_vPrimeMean(1:length(frequency),j),'-b','LineWidth',2)
    hold on;
    grid on;
    xlabel('\bffrequency, \itf \rm\bf(Hz)')
    ylabel('\bf|F|')
    title(['vPrime PSD with Hanning window, probe ',num2str(j-2)],'FontSize',12)
    hold off;
    
    subplot(4,2,6)
    semilogx(frequency,Magnitude_Hann_wPrimeMean(1:length(frequency),j),'-g','LineWidth',2)
    hold on;
    grid on;
    xlabel('\bffrequency, \itf \rm\bf(Hz)')
    ylabel('\bf|F|')
    title(['wPrime PSD with Hanning window, probe ',num2str(j-2)],'FontSize',12)
    hold off;
    
    subplot(4,2,8)
    semilogx(frequency,Magnitude_Hann_TPrimeMean(1:length(frequency),j),'-y','LineWidth',2)
    hold on;
    grid on;
    xlabel('\bffrequency, \itf \rm\bf(Hz)')
    ylabel('\bf|F|')
    title(['TPrime PSD with Hanning window, probe',num2str(j-2)],'FontSize',12)
    hold off;
    
    set(gcf, 'Position',  [200, 50, 1000, 700]) 
end
