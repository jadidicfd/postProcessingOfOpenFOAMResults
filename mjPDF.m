%% Header
% @author: M.Jadidi

% @program name: pdf
% @dependency: OpenFoam data in the format of "probes in OpenFOAM" MUST be... 
% availabe in the working Dir. 
% @task: calculate and dispay pdf, Skewness, Flatness of time signal from OpenFOAM


% @see: [1]	M. Jadidi, H. K. Param, A. Revell, and Y. Mahmoudi, "Flow leakage and Kelvin–Helmholtz instability of turbulent flow over porous media," Physics of Fluids, vol. 34, no. 10, p. 105114, 2022/10/01 2022, doi: 10.1063/5.0111195.
% @see: [2]	M. Jadidi, H. K. Param, A. Revell, and Y. Mahmoudi, "Large eddy simulations of turbulent heat transfer in packed bed energy storage systems," Journal of Energy Storage, vol. 59, p. 106449, 2023/03/01/ 2023, doi: https://doi.org/10.1016/j.est.2022.106449.
% @see: [3]	M. Jadidi, A. Revell, and Y. Mahmoudi, "Pore-scale large eddy simulation of turbulent flow and heat transfer over porous media," Applied Thermal Engineering, vol. 215, p. 118916, 2022/10/01/ 2022, doi: https://doi.org/10.1016/j.applthermaleng.2022.118916.
% @see: https://uk.mathworks.com/matlabcentral/fileexchange/17204-kernel-density-estimation


% @created: May 2021,
% @version 01 (May 2021)
% @version 02 (June 2021)

%*******************IMPORTANT NOTE*****************************************
% The variable "nBins" should be declared by the user
%**************************************************************************

close all;
%% Probability density/histogram
close all

% n1=floor(1+3.3*log(size(uPrime,1)));  %Sturgis rule
% n2=floor(2*(size(uPrime,1))^(1/3));  %Sturgis rule
% nBins=floor((n1+n2)/2);

nBins=13;

norm_uPrime=zeros(1,c);
norm_vPrime=zeros(1,c);
norm_wPrime=zeros(1,c);
norm_TPrime=zeros(1,c);

Skewness_uPrime=zeros(1,c);
Skewness_vPrime=zeros(1,c);
Skewness_wPrime=zeros(1,c);
Skewness_TPrime=zeros(1,c);

Flatness_uPrime=zeros(1,c);
Flatness_vPrime=zeros(1,c);
Flatness_wPrime=zeros(1,c);
Flatness_TPrime=zeros(1,c);

pdf_uPrime=cell(1,c);
pdf_vPrime=cell(1,c);
pdf_wPrime=cell(1,c);
pdf_TPrime=cell(1,c);

pdf_uPrime_normalized=cell(1,c);
pdf_vPrime_normalized=cell(1,c);
pdf_wPrime_normalized=cell(1,c);
pdf_TPrime_normalized=cell(1,c);

for j=2:c
    [pdf_uPrime{1,j}(:,1),pdf_uPrime{1,j}(:,2)]= hist(uPrime(:,j),nBins);
    [pdf_vPrime{1,j}(:,1),pdf_vPrime{1,j}(:,2)]= hist(vPrime(:,j),nBins);
    [pdf_wPrime{1,j}(:,1),pdf_wPrime{1,j}(:,2)]= hist(wPrime(:,j),nBins);
    [pdf_TPrime{1,j}(:,1),pdf_TPrime{1,j}(:,2)]= hist(TPrime(:,j),nBins);
    
    
    norm_uPrime(1,j)=sum(pdf_uPrime{1,j}(:,1))*(pdf_uPrime{1,j}(3,2)-pdf_uPrime{1,j}(2,2));
    norm_vPrime(1,j)=sum(pdf_vPrime{1,j}(:,1))*(pdf_vPrime{1,j}(3,2)-pdf_vPrime{1,j}(2,2));
    norm_wPrime(1,j)=sum(pdf_wPrime{1,j}(:,1))*(pdf_wPrime{1,j}(3,2)-pdf_wPrime{1,j}(2,2));
    norm_TPrime(1,j)=sum(pdf_TPrime{1,j}(:,1))*(pdf_TPrime{1,j}(3,2)-pdf_TPrime{1,j}(2,2));
    
    
    pdf_uPrime_normalized{1,j}(:,1)=pdf_uPrime{1,j}(:,1)/norm_uPrime(1,j);
    pdf_vPrime_normalized{1,j}(:,1)=pdf_vPrime{1,j}(:,1)/norm_vPrime(1,j);
    pdf_wPrime_normalized{1,j}(:,1)=pdf_wPrime{1,j}(:,1)/norm_wPrime(1,j);
    pdf_TPrime_normalized{1,j}(:,1)=pdf_TPrime{1,j}(:,1)/norm_TPrime(1,j);
    
    pdf_uPrime_normalized{1,j}(:,2)=pdf_uPrime{1,j}(:,2);
    pdf_vPrime_normalized{1,j}(:,2)=pdf_vPrime{1,j}(:,2);
    pdf_wPrime_normalized{1,j}(:,2)=pdf_wPrime{1,j}(:,2);
    pdf_TPrime_normalized{1,j}(:,2)=pdf_TPrime{1,j}(:,2);
    
    figure(j)
    set(gcf, 'Position',  [400, 50, 400, 700])
    subplot(4,1,1)
    hold on;
    plot(pdf_uPrime_normalized{1,j}(:,2),pdf_uPrime_normalized{1,j}(:,1),'-r','linewidth',2)
    grid on;
    xlabel('uPrime','FontSize',12)
    ylabel('pdf','FontSize',12)
    title(['pdf uPrime for probe: ',num2str(j-2)],'FontSize',12)
    hold off;
    
    subplot(4,1,2)
    hold on;
    plot(pdf_vPrime_normalized{1,j}(:,2),pdf_vPrime_normalized{1,j}(:,1),'-b','linewidth',2)
    grid on;
    xlabel('vPrime','FontSize',12)
    ylabel('pdf','FontSize',12)
    title(['pdf vPrime for probe: ',num2str(j-2)],'FontSize',12)
    hold off;
    
    subplot(4,1,3)
    hold on;
    plot(pdf_wPrime_normalized{1,j}(:,2),pdf_wPrime_normalized{1,j}(:,1),'-g','linewidth',2)
    grid on;
    xlabel('u: ','FontSize',12)
    ylabel('pdf','FontSize',12)
    title(['pdf wPrime for probe: ',num2str(j-2)],'FontSize',12)
    hold off;
    
    subplot(4,1,4)
    hold on;
    plot(pdf_TPrime_normalized{1,j}(:,2),pdf_TPrime_normalized{1,j}(:,1),'-y','linewidth',2)
    grid on;
    xlabel('TPrime','FontSize',12)
    ylabel('pdf','FontSize',12)
    title(['pdf TPrime for probe: ',num2str(j-2)],'FontSize',12)
    hold off;
    
    
    %Skewness and flatness
    Skewness_uPrime(1,j)=mean(uPrime(:,j).^3)/uRMS(1,j)^3;
    Skewness_vPrime(1,j)=mean(vPrime(:,j).^3)/vRMS(1,j)^3;
    Skewness_wPrime(1,j)=mean(wPrime(:,j).^3)/wRMS(1,j)^3;
    Skewness_TPrime(1,j)=mean(TPrime(:,j).^3)/TRMS(1,j)^3;
    
    Flatness_uPrime(1,j)=mean(uPrime(:,j).^4)/uRMS(1,j)^4;
    Flatness_vPrime(1,j)=mean(vPrime(:,j).^4)/vRMS(1,j)^4;
    Flatness_wPrime(1,j)=mean(wPrime(:,j).^4)/wRMS(1,j)^4;
    Flatness_TPrime(1,j)=mean(TPrime(:,j).^4)/TRMS(1,j)^4;
end

