%% Header
% @author: M.Jadidi

% @program name: 3D histogram and JPDF surface and contour
% @dependency: OpenFoam data in the format of "probes in OpenFOAM" MUST be... 
% availabe in the working Dir.
% @dependency: mjReadingData.m

% @task: Calculate 2D and 3D historan, compute JPDF and dispay in format of contour and surface

% @see: [1]	M. Jadidi, H. K. Param, A. Revell, and Y. Mahmoudi, "Flow leakage and Kelvin–Helmholtz instability of turbulent flow over porous media," Physics of Fluids, vol. 34, no. 10, p. 105114, 2022/10/01 2022, doi: 10.1063/5.0111195.
% @see: [2]	M. Jadidi, H. K. Param, A. Revell, and Y. Mahmoudi, "Large eddy simulations of turbulent heat transfer in packed bed energy storage systems," Journal of Energy Storage, vol. 59, p. 106449, 2023/03/01/ 2023, doi: https://doi.org/10.1016/j.est.2022.106449.
% @see: [3]	M. Jadidi, A. Revell, and Y. Mahmoudi, "Pore-scale large eddy simulation of turbulent flow and heat transfer over porous media," Applied Thermal Engineering, vol. 215, p. 118916, 2022/10/01/ 2022, doi: https://doi.org/10.1016/j.applthermaleng.2022.118916.
% @see: https://uk.mathworks.com/matlabcentral/fileexchange/17204-kernel-density-estimation


% @created: May 2021,
% @version 01 (May 2021)
% @version 02 (June 2021)

%*******************IMPORTANT NOTE*****************************************

% The variables "x_pdf" and "y_pdf" should be declared by the user
% "x_pdf" and "y_pdf" control the number of bins for Histogram and JPDF

%**************************************************************************

close all;

x_pdf = -3:1:3; % axis x, which you want to see for pdf (surf and contour)
y_pdf = -3:1:3; % axis y, which you want to see for pdf (surf and contour)
[X_pdf,Y_pdf] = meshgrid(x_pdf,y_pdf); %important for "surf" - makes defined grid

pdf_uv = cell(1,c);
pdf_uv_normalize = cell(1,c);
integralOverDensityPlot_uv = zeros(1,c);
for j=3:c
 
    pdf_uv{1,j} = hist3([uSigma(:,j) , vSigma(:,j)],{x_pdf y_pdf}); %standard hist3 (calculated for yours axis)
    pdf_uv_normalize{1,j} = ((pdf_uv{1,j})'./length(vSigma(:,j)));    %normalization means devide it by length of
    
    %This give the joint density plot in 3D. Which can be checked by calculating the integral over the surface with.
    %the following integral should be one:
    integralOverDensityPlot_uv(1,j) = sum(trapz(pdf_uv_normalize{1,j}));
end

pdf_uw = cell(1,c);
pdf_uw_normalize = cell(1,c);
integralOverDensityPlot_uw = zeros(1,c);
for j=3:c
    
    pdf_uw{1,j} = hist3([uSigma(:,j) , wSigma(:,j)],{x_pdf y_pdf}); %standard hist3 (calculated for yours axis)
    pdf_uw_normalize{1,j} = ((pdf_uw{1,j})'./length(uSigma(:,j)));    %normalization means devide it by length of
    
    %This give the joint density plot in 3D. Which can be checked by calculating the integral over the surface with.
    %the following integral should be one:
    integralOverDensityPlot_uw(1,j) = sum(trapz(pdf_uw_normalize{1,j}));
end

pdf_vw = cell(1,c);
pdf_vw_normalize = cell(1,c);
integralOverDensityPlot_vw = zeros(1,c);
for j=3:c
 
    pdf_vw{1,j} = hist3([vSigma(:,j) , wSigma(:,j)],{x_pdf y_pdf}); %standard hist3 (calculated for yours axis)
    pdf_vw_normalize{1,j} = ((pdf_vw{1,j})'./length(vSigma(:,j)));    %normalization means devide it by length of
    
    %This give the joint density plot in 3D. Which can be checked by calculating the integral over the surface with.
    %the following integral should be one:
    integralOverDensityPlot_vw(1,j) = sum(trapz(pdf_vw_normalize{1,j}));
end

pdf_uT = cell(1,c);
pdf_uT_normalize = cell(1,c);
integralOverDensityPlot_uT = zeros(1,c);
for j=3:c
 
    pdf_uT{1,j} = hist3([uSigma(:,j) , TSigma(:,j)],{x_pdf y_pdf}); %standard hist3 (calculated for yours axis)
    pdf_uT_normalize{1,j} = ((pdf_uT{1,j})'./length(uSigma(:,j)));    %normalization means devide it by length of
    
    %This give the joint density plot in 3D. Which can be checked by calculating the integral over the surface with.
    %the following integral should be one:
    integralOverDensityPlot_uT(1,j) = sum(trapz(pdf_uT_normalize{1,j}));
end

pdf_vT = cell(1,c);
pdf_vT_normalize = cell(1,c);
integralOverDensityPlot_vT = zeros(1,c);
for j=3:c

    pdf_vT{1,j} = hist3([vSigma(:,j) , TSigma(:,j)],{x_pdf y_pdf}); %standard hist3 (calculated for yours axis)
    pdf_vT_normalize{1,j} = ((pdf_vT{1,j})'./length(vSigma(:,j)));    %normalization means devide it by length of
    
    %This give the joint density plot in 3D. Which can be checked by calculating the integral over the surface with.
    %the following integral should be one:
    integralOverDensityPlot_vT(1,j) = sum(trapz(pdf_vT_normalize{1,j}));
end

pdf_wT = cell(1,c);
pdf_wT_normalize = cell(1,c);
integralOverDensityPlot_wT = zeros(1,c);
for j=3:c

    pdf_wT{1,j} = hist3([wSigma(:,j) , TSigma(:,j)],{x_pdf y_pdf}); %standard hist3 (calculated for yours axis)
    pdf_wT_normalize{1,j} = ((pdf_wT{1,j})'./length(wSigma(:,j)));    %normalization means devide it by length of
    
    %This give the joint density plot in 3D. Which can be checked by calculating the integral over the surface with.
    %the following integral should be one:
    integralOverDensityPlot_wT(1,j) = sum(trapz(pdf_wT_normalize{1,j}));
end

%% figures

% historam 3D
for j=3:c
    figure('Name',['--> hist3 for probe: ',num2str(j-2)]);
    set(gcf, 'Position',  [50, 50, 1400, 700])
    
    subplot(2,3,1)
    hist3([uSigma(:,j),vSigma(:,j)],'Nbins',[10 10],'CDataMode','auto','FaceColor','interp')
    xlabel('uSigma')
    ylabel('vSigma')
    title(['Histogram of uv for probe: ',num2str(j-2)],'FontSize',12)
    colorbar;
    box on
    
    subplot(2,3,2)
    hist3([uSigma(:,j),wSigma(:,j)],'Nbins',[10 10],'CDataMode','auto','FaceColor','interp')
    xlabel('uSigma')
    ylabel('wSigma')
    title(['Histogram of uw for probe: ',num2str(j-2)],'FontSize',12)
    colorbar;
    box on
    
    subplot(2,3,3)
    hist3([vSigma(:,j),wSigma(:,j)],'Nbins',[10 10],'CDataMode','auto','FaceColor','interp')
    xlabel('vSigma')
    ylabel('wSigma')
    title(['Histogram of vw for probe: ',num2str(j-2)],'FontSize',12)
    colorbar;
    box on
    
    subplot(2,3,4)
    hist3([uSigma(:,j),TSigma(:,j)],'Nbins',[10 10],'CDataMode','auto','FaceColor','interp')
    xlabel('uSigma')
    ylabel('TSigma')
    title(['Histogram of uT for probe: ',num2str(j-2)],'FontSize',12)
    colorbar;
    box on
    
    subplot(2,3,5)
    hist3([vSigma(:,j),TSigma(:,j)],'Nbins',[10 10],'CDataMode','auto','FaceColor','interp')
    xlabel('vSigma')
    ylabel('TSigma')
    title(['Histogram of vT for probe: ',num2str(j-2)],'FontSize',12)
    colorbar;
    box on
    
    subplot(2,3,6)
    hist3([wSigma(:,j),TSigma(:,j)],'Nbins',[10 10],'CDataMode','auto','FaceColor','interp')
    xlabel('wSigma')
    ylabel('TSigma')
    title(['Histogram of wT for probe: ',num2str(j-2)],'FontSize',12)
    colorbar;
    box on
    
end

% top view of normalized 3Dhistoram
for j=3:c
    figure('Name',['--> surf for probe: ',num2str(j-2)]);
    set(gcf, 'Position',  [50, 50, 1400, 700])
    
    subplot(2,3,1)
    surf(X_pdf,Y_pdf,pdf_uv_normalize{1,j}) % plot distribution
    xlabel('uSigma')
    ylabel('vSigma')
    title(['normalized pdf uv for probe: ',num2str(j-2)],'FontSize',12)
    colorbar;
    box on
    view(2)
    
    subplot(2,3,2)
    surf(X_pdf,Y_pdf,pdf_uw_normalize{1,j}) % plot distribution
    xlabel('uSigma')
    ylabel('wSigma')
    title(['normalized pdf uw for probe: ',num2str(j-2)],'FontSize',12)
    colorbar;
    box on
    view(2)
    
    subplot(2,3,3)
    surf(X_pdf,Y_pdf,pdf_vw_normalize{1,j}) % plot distribution
    xlabel('vSigma')
    ylabel('wSigma')
    title(['normalized pdf vw for probe: ',num2str(j-2)],'FontSize',12)
    colorbar;
    box on
    view(2)
    
    subplot(2,3,4)
    surf(X_pdf,Y_pdf,pdf_uT_normalize{1,j}) % plot distribution
    xlabel('uSigma')
    ylabel('TSigma')
    title(['normalized pdf uT for probe: ',num2str(j-2)],'FontSize',12)
    colorbar;
    box on
    view(2)
    
    subplot(2,3,5)
    surf(X_pdf,Y_pdf,pdf_vT_normalize{1,j}) % plot distribution
    xlabel('vSigma')
    ylabel('TSigma')
    title(['normalized pdf vT for probe: ',num2str(j-2)],'FontSize',12)
    colorbar;
    box on
    view(2)
    
    subplot(2,3,6)
    surf(X_pdf,Y_pdf,pdf_wT_normalize{1,j}) % plot distribution
    xlabel('wSigma')
    ylabel('TSigma')
    title(['normalized pdf wT for probe: ',num2str(j-2)],'FontSize',12)
    colorbar;
    box on
    view(2)
    
end


for j=3:c
    figure('Name',['--> contour for probe: ',num2str(j-2)]);
    set(gcf, 'Position',  [50, 50, 1400, 700])
    
    subplot(2,3,1)
    contour(X_pdf,Y_pdf,pdf_uv_normalize{1,j}) % plot distribution
    xlabel('uSigma')
    ylabel('vSigma')
    title(['normalized pdf uv for probe: ',num2str(j-2)],'FontSize',12)
    colorbar;
    box on
    
    subplot(2,3,2)
    contour(X_pdf,Y_pdf,pdf_uw_normalize{1,j}) % plot distribution
    xlabel('uSigma')
    ylabel('wSigma')
    title(['normalized pdf uw for probe: ',num2str(j-2)],'FontSize',12)
    colorbar;
    box on
    
    subplot(2,3,3)
    contour(X_pdf,Y_pdf,pdf_vw_normalize{1,j}) % plot distribution
    xlabel('vSigma')
    ylabel('wSigma')
    title(['normalized pdf vw for probe: ',num2str(j-2)],'FontSize',12)
    colorbar;
    box on
    
    subplot(2,3,4)
    contour(X_pdf,Y_pdf,pdf_uT_normalize{1,j}) % plot distribution
    xlabel('uSigma')
    ylabel('TSigma')
    title(['normalized pdf uT for probe: ',num2str(j-2)],'FontSize',12)
    colorbar;
    box on
    
    subplot(2,3,5)
    contour(X_pdf,Y_pdf,pdf_vT_normalize{1,j}) % plot distribution
    xlabel('vSigma')
    ylabel('TSigma')
    title(['normalized pdf vT for probe: ',num2str(j-2)],'FontSize',12)
    colorbar;
    box on
    
    subplot(2,3,6)
    contour(X_pdf,Y_pdf,pdf_wT_normalize{1,j}) % plot distribution
    xlabel('wSigma')
    ylabel('TSigma')
    title(['normalized pdf wT for probe: ',num2str(j-2)],'FontSize',12)
    colorbar;
    box on
end