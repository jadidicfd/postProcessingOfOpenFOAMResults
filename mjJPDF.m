%% Header
% @author: M.Jadidi

% @program name: JPDF
% @dependency: OpenFoam data in the format of "probes in OpenFOAM" MUST be availabe in the working Dir. 
% @dependency: mjReadingData.m
% @dependency: mjKDE2d.m ("n" in this fuction can be calibrated for better/smooth graph)

% @task: sort, calculate and dispay JPDF. Compute the magnitude of
% correlation between two signals

% @see: [1]	M. Jadidi, H. K. Param, A. Revell, and Y. Mahmoudi, "Flow leakage and Kelvin–Helmholtz instability of turbulent flow over porous media," Physics of Fluids, vol. 34, no. 10, p. 105114, 2022/10/01 2022, doi: 10.1063/5.0111195.
% @see: [2]	M. Jadidi, H. K. Param, A. Revell, and Y. Mahmoudi, "Large eddy simulations of turbulent heat transfer in packed bed energy storage systems," Journal of Energy Storage, vol. 59, p. 106449, 2023/03/01/ 2023, doi: https://doi.org/10.1016/j.est.2022.106449.
% @see: [3]	M. Jadidi, A. Revell, and Y. Mahmoudi, "Pore-scale large eddy simulation of turbulent flow and heat transfer over porous media," Applied Thermal Engineering, vol. 215, p. 118916, 2022/10/01/ 2022, doi: https://doi.org/10.1016/j.applthermaleng.2022.118916.
% @see: https://uk.mathworks.com/matlabcentral/fileexchange/17204-kernel-density-estimation


% @created: May 2021,
% @version 01 (May 2021)
% @version 02 (June 2021)

close all;

%% Sorting probes data for calculating the JPDF for u'v'

uvProbe = cell(1,c);             % returns an n-by-n cell array of empty matrices for initialization
corruv = zeros(1,c);
for j=4:c
    uvProbe{1,j}(:,1) = uSigma(:,j);
    uvProbe{1,j}(:,2) = vSigma(:,j);
    corruv(1,j)=corr(uSigma(:,j),vSigma(:,j));
end

density_uvProbe = cell(1,c);       % returns an n-by-n cell array of empty matrices for initialization
for j=4:c
   
    figure(j);
   
    set(gcf, 'Position',  [50, 50, 1400, 700])
    subplot(5,3,1)
    hold on;
    
    %call the routine
    [bandwidth,density,X,Y]=mjKDE2d(uvProbe{1,j});
    
    %plot the data and the density estimate
    contour3(X,Y,density,10),
    hold on
    %plot(uvProbe{1,j}(:,1),uvProbe{1,j}(:,2),'r.','MarkerSize',5)
    colorbar;
    xlabel('uSigma'); ylabel('vSigma');
    title(['uv JPDF for probe:',num2str(j-2)])
    
    density_uvProbe{1,j} = density;
end

%% Sorting probes data for calculating the JPDF for u'w'

uwProbe = cell(1,c);             % returns an n-by-n cell array of empty matrices for initialization
corruw = zeros(1,c);
for j=4:c
    uwProbe{1,j}(:,1) = uSigma(:,j);
    uwProbe{1,j}(:,2) = wSigma(:,j);
    corruw(1,j)=corr(uSigma(:,j),wSigma(:,j));
end

density_uwProbe = cell(1,c);       % returns an n-by-n cell array of empty matrices for initialization
for j=4:c
   
    figure(j);
    subplot(5,3,2)
    hold on;
    
    %call the routine
    [bandwidth,density,X,Y]=mjKDE2d(uwProbe{1,j});
    
    %plot the data and the density estimate
    contour(X,Y,density,10), hold on
    %plot(uwProbe{1,j}(:,1),uwProbe{1,j}(:,2),'r.','MarkerSize',5)
    colorbar;
    xlabel('uSigma'); ylabel('wSigma');
    title(['uw JPDF for probe:',num2str(j-2)])
    
    density_uwProbe{1,j} = density;
end

%% Sorting probes data for calculating the JPDF for v'w'

vwProbe = cell(1,c);             % returns an n-by-n cell array of empty matrices for initialization
corrvw = zeros(1,c);
for j=4:c
    vwProbe{1,j}(:,1) = vSigma(:,j);
    vwProbe{1,j}(:,2) = wSigma(:,j);
    corrvw(1,j)=corr(vSigma(:,j),wSigma(:,j));
end

density_vwProbe = cell(1,c);       % returns an n-by-n cell array of empty matrices for initialization
for j=4:c
   
    figure(j);
    subplot(5,3,3)
    hold on;
    
    %call the routine
    [bandwidth,density,X,Y]=mjKDE2d(vwProbe{1,j});
    
    %plot the data and the density estimate
    contour(X,Y,density,10), hold on
    %plot(uwProbe{1,j}(:,1),uwProbe{1,j}(:,2),'r.','MarkerSize',5)
    colorbar;
    xlabel('vSigma'); ylabel('wSigma');
    title(['vw JPDF for probe:',num2str(j-2)])
    
    density_vwProbe{1,j} = density;
end
%% Sorting probes data for calculating the JPDF for u'T'

uTProbe = cell(1,c);             % returns an n-by-n cell array of empty matrices for initialization
corruT = zeros(1,c);
for j=4:c
    uTProbe{1,j}(:,1) = uSigma(:,j);
    uTProbe{1,j}(:,2) = TSigma(:,j);
    corruT(1,j)=corr(uSigma(:,j),TSigma(:,j));
end

density_uTProbe = cell(1,c);       % returns an n-by-n cell array of empty matrices for initialization
for j=4:c
   
    figure(j);
    subplot(5,3,4)
    hold on;
    
    %call the routine
    [bandwidth,density,X,Y]=mjKDE2d(uTProbe{1,j});
    
    %plot the data and the density estimate
    contour(X,Y,density,10), hold on
    %plot(uvProbe{1,j}(:,1),uvProbe{1,j}(:,2),'r.','MarkerSize',5)
    colorbar;
    xlabel('uSigma'); ylabel('TSigma');
    title(['uT JPDF for probe:',num2str(j-2)])
    
    density_uTProbe{1,j} = density;
end

%% Sorting probes data for calculating the JPDF for v'T'
vTProbe = cell(1,c);             % returns an n-by-n cell array of empty matrices for initialization
corrvT = zeros(1,c);
for j=4:c
    vTProbe{1,j}(:,1) = vSigma(:,j);
    vTProbe{1,j}(:,2) = TSigma(:,j);
    corrvT(1,j)=corr(vSigma(:,j),TSigma(:,j));
end

density_vTProbe = cell(1,c);       % returns an n-by-n cell array of empty matrices for initialization
for j=4:c

    figure(j);
    subplot(5,3,5)
    hold on;
    
    %call the routine
    [bandwidth,density,X,Y]=mjKDE2d(vTProbe{1,j});
    colorbar;
    %plot the data and the density estimate
    contour(X,Y,density,10), hold on
    %plot(vTProbe{1,j}(:,1),vTProbe{1,j}(:,2),'b.','MarkerSize',2)
    colorbar;
    xlabel('TSigma'); ylabel('vSigma');
    title(['vT JPDF for probe:',num2str(j-2)])
    hold off;
    density_vTProbe{1,j} = density;
end

%% Sorting probes data for calculating the JPDF for w'T'
wTProbe = cell(1,c);             % returns an n-by-n cell array of empty matrices for initialization
corrwT = zeros(1,c);
for j=4:c
    wTProbe{1,j}(:,1) = wSigma(:,j);
    wTProbe{1,j}(:,2) = TSigma(:,j);
    corrwT(1,j)=corr(wSigma(:,j),TSigma(:,j));
end

density_wTProbe = cell(1,c);       % returns an n-by-n cell array of empty matrices for initialization
for j=4:c

    figure(j);
    subplot(5,3,6)
    hold on;
    
    %call the routine
    [bandwidth,density,X,Y]=mjKDE2d(wTProbe{1,j});
    colorbar;
    %plot the data and the density estimate
    contour(X,Y,density,10), hold on
    %plot(vTProbe{1,j}(:,1),vTProbe{1,j}(:,2),'b.','MarkerSize',2)
    colorbar;
    xlabel('TSigma'); ylabel('wSigma');
    title(['wT JPDF for probe:',num2str(j-2)])
    hold off;
    density_wTProbe{1,j} = density;
end
%% Sorting probes data for calculating the JPDF for u'v', u'T' and v'T' when T'>0

uSigmaPosetive = cell (1,c);
vSigmaPosetive = cell (1,c);
TSigmaPosetive = cell (1,c);

for j=4:c
    ind = find(TSigma(:,j)>0);
    TSigmaPosetive{1,j} = [ind,TSigma(ind,j)];
    uSigmaPosetive{1,j} = [ind,uSigma(ind,j)];
    vSigmaPosetive{1,j} = [ind,vSigma(ind,j)];
end

%JPDF for u'v' when T'>0

probeuSigmaPosetivevSigmaPosetive = cell(1,c);             % returns an n-by-n cell array of empty matrices for initialization
corruSigmaPosetivevSigmaPosetive = zeros(1,c);

for j=4:c
    probeuSigmaPosetivevSigmaPosetive{1,j}(:,1) = uSigmaPosetive{1,j}(:,2);
    probeuSigmaPosetivevSigmaPosetive{1,j}(:,2) = vSigmaPosetive{1,j}(:,2);
    corruSigmaPosetivevSigmaPosetive(1,j)=corr(uSigmaPosetive{1,j}(:,2),vSigmaPosetive{1,j}(:,2));
end

densityProbeuSigmaPosetivevSigmaPosetive = cell(1,c);       % returns an n-by-n cell array of empty matrices for initialization
for j=4:c
    
    figure(j);
    subplot(5,3,7)
    hold on;
    
    %call the routine
    [bandwidth,density,X,Y]=mjKDE2d(probeuSigmaPosetivevSigmaPosetive{1,j});
    
    %plot the data and the density estimate
    contour3(X,Y,density,10), hold on
    %plot(uSigmaPosetive{1,j}(:,2),vSigmaPosetive{1,j}(:,2),'g.','MarkerSize',2)
    colorbar;
    xlabel('uSigma'); ylabel('vSigma');
    title(['uv JPDF when TPrime>0 for probe:',num2str(j-2)])
    hold off;
    densityProbeuSigmaPosetivevSigmaPosetive{1,j} = density;
end



% JPDF for u'T' when T'>0

probeTSigmaPosetiveuSigmaPosetive = cell(1, c);             % returns an n-by-n cell array of empty matrices for initialization
corrTSigmaPosetiveuSigmaPosetive = zeros(1,c);
for j=4:c
    probeTSigmaPosetiveuSigmaPosetive{1,j}(:,1) = TSigmaPosetive{1,j}(:,2);
    probeTSigmaPosetiveuSigmaPosetive{1,j}(:,2) = uSigmaPosetive{1,j}(:,2);
    corrTSigmaPosetiveuSigmaPosetive(1,j)=corr(TSigmaPosetive{1,j}(:,2),uSigmaPosetive{1,j}(:,2));
end

densityProbeTSigmaPosetiveuSigmaPosetive = cell(1,c);       % returns an n-by-n cell array of empty matrices for initialization
for j=4:c

    figure(j);
    subplot(5,3,8)
    hold on; 
    
    %call the routine
    [bandwidth,density,X,Y]=mjKDE2d(probeTSigmaPosetiveuSigmaPosetive{1,j});
    
    %plot the data and the density estimate
    contour3(X,Y,density,10), hold on
    %plot(TSigmaPosetive{1,j}(:,2),uSigmaPosetive{1,j}(:,2),'r.','MarkerSize',2)
    colorbar;
    xlabel('TSigmaPosetive'); ylabel('uSigma');
    title([" JPDF of u'T' when T'>0  for probe:",num2str(j-1)])
    
    densityProbeTSigmaPosetiveuSigmaPosetive{1,j} = density;
end


%JPDF for v'T' when T'>0

probeTSigmaPosetivevSigmaPosetive = cell(1, c);             % returns an n-by-n cell array of empty matrices for initialization
corrTSigmaPosetivevSigmaPosetive = zeros(1,c);
for j=4:c
    probeTSigmaPosetivevSigmaPosetive{1,j}(:,1) = TSigmaPosetive{1,j}(:,2);
    probeTSigmaPosetivevSigmaPosetive{1,j}(:,2) = vSigmaPosetive{1,j}(:,2);
    corrTSigmaPosetivevSigmaPosetive(1,j)=corr(TSigmaPosetive{1,j}(:,2),vSigmaPosetive{1,j}(:,2));

end

densityProbeTSigmaPosetivevSigmaPosetive = cell(1,c);       % returns an n-by-n cell array of empty matrices for initialization
for j=4:c
    figure(j);
    subplot(5,3,9)
    hold on; 
    %call the routine
    [bandwidth,density,X,Y]=mjKDE2d(probeTSigmaPosetivevSigmaPosetive{1,j});
    
    %plot the data and the density estimate
    contour3(X,Y,density,10), hold on
    %plot(TSigmaPosetive{1,j}(:,2),vSigmaPosetive{1,j}(:,2),'r.','MarkerSize',2)
    colorbar;
    xlabel('TSigmaPosetive'); ylabel('vSigma');
    title([" JPDF of u'T' when T'>0  for probe:",num2str(j-1)])
    
    densityProbeTSigmaPosetivevSigmaPosetive{1,j} = density;
end


%% Sorting probes data for calculating the JPDF for u'v', u'T' and v'T' when T'<0

TSigmaNegative = cell (1,c);
uSigmaNegative = cell (1,c);
vSigmaNegative = cell (1,c);

for j=4:c
    ind = find(TSigma(:,j)<0);
    TSigmaNegative{1,j} = [ind,TSigma(ind,j)];
    uSigmaNegative{1,j} = [ind,uSigma(ind,j)];
    vSigmaNegative{1,j} = [ind,vSigma(ind,j)];
end


%JPDF for u'v' when T'<0 
probeuSigmaNegativevSigmaNegative = cell(1,c);             % returns an n-by-n cell array of empty matrices for initialization
corruSigmaNegativevSigmaNegative = zeros(1,c);

for j=4:c
    probeuSigmaNegativevSigmaNegative{1,j}(:,1) = uSigmaNegative{1,j}(:,2);
    probeuSigmaNegativevSigmaNegative{1,j}(:,2) = vSigmaNegative{1,j}(:,2);
    corruSigmaNegativevSigmaNegative(1,j)=corr(uSigmaNegative{1,j}(:,2),vSigmaNegative{1,j}(:,2));
end

densityProbeuSigmaNegativevSigmaNegative = cell(1,c);       % returns an n-by-n cell array of empty matrices for initialization
for j=4:c
    figure(j);
    subplot(5,3,10)
    hold on;
    %call the routine
    [bandwidth,density,X,Y]=mjKDE2d(probeuSigmaNegativevSigmaNegative{1,j});
    
    %plot the data and the density estimate
    contour3(X,Y,density,10), hold on
    %plot(uSigmaNegative{1,j}(:,2),vSigmaNegative{1,j}(:,2),'k.','MarkerSize',2)
    colorbar;
    xlabel('uSigma'); ylabel('vSigma');
    title(['uv JPDF when TPrime<0 for probe:',num2str(j-2)])
    hold off;
    densityProbeuSigmaNegativevSigmaNegative{1,j} = density;
end



%JPDF for u'T' when T'<0 
probeTSigmaNegativeuSigmaNegative = cell(1, c);             % returns an n-by-n cell array of empty matrices for initialization
corrTSigmaNegativeuSigmaNegative = zeros(1,c);
for j=4:c
    probeTSigmaNegativeuSigmaNegative{1,j}(:,1) = TSigmaNegative{1,j}(:,2);
    probeTSigmaNegativeuSigmaNegative{1,j}(:,2) = uSigmaNegative{1,j}(:,2);
    corrTSigmaNegativeuSigmaNegative(1,j)=corr(TSigmaNegative{1,j}(:,2),uSigmaNegative{1,j}(:,2));
end

densityProbeTSigmaNegativeuSigmaNegative = cell(1,c);       % returns an n-by-n cell array of empty matrices for initialization
for j=4:c
    figure(j);
    subplot(5,3,11)
    hold on;  
    %call the routine
    [bandwidth,density,X,Y]=mjKDE2d(probeTSigmaNegativeuSigmaNegative{1,j});
    
    %plot the data and the density estimate
    contour3(X,Y,density,10), hold on
    %plot(TSigmaNegative{1,j}(:,2),uSigmaNegative{1,j}(:,2),'r.','MarkerSize',2)
    colorbar;
    xlabel('TSigmaNegative'); ylabel('uSigma');
    title([" JPDF of u'T' when T'<0  for probe:",num2str(j-1)])
    
    densityProbeTSigmaNegativeuSigmaNegative{1,j} = density;
end
% 
%JPDF for v'T' when T'<0
probeTSigmaNegativevSigmaNegative = cell(1, c);             % returns an n-by-n cell array of empty matrices for initialization
corrTSigmaNegativevSigmaNegative = zeros(1,c);
for j=4:c
    probeTSigmaNegativevSigmaNegative{1,j}(:,1) = TSigmaNegative{1,j}(:,2);
    probeTSigmaNegativevSigmaNegative{1,j}(:,2) = vSigmaNegative{1,j}(:,2);
    corrTSigmaNegativevSigmaNegative(1,j)=corr(TSigmaNegative{1,j}(:,2),vSigmaNegative{1,j}(:,2));
end

densityProbeTSigmaPosetivevSigmaPosetive = cell(1,c);       % returns an n-by-n cell array of empty matrices for initialization
for j=4:c
    figure(j);
    subplot(5,3,12)
    hold on;    
    %call the routine
    [bandwidth,density,X,Y]=mjKDE2d(probeTSigmaNegativevSigmaNegative{1,j});
    
    %plot the data and the density estimate
    contour3(X,Y,density,10), hold on
    %plot(TSigmaNegative{1,j}(:,2),vSigmaNegative{1,j}(:,2),'r.','MarkerSize',2)
    colorbar;
    xlabel('TSigmaNegative'); ylabel('vSigma');
    title([" JPDF of v'T' when T'<0  for probe:",num2str(j-1)])
    
    densityProbeTSigmaPosetivevSigmaPosetive{1,j} = density;
end

%% Sorting probes data for calculating the JPDF for u'p', v'p' and v'p'

%JPDF for u'p'

upProbe = cell(1,c);             % returns an n-by-n cell array of empty matrices for initialization
corrup = zeros(1,c);
for j=4:c
    upProbe{1,j}(:,1) = uSigma(:,j);
    upProbe{1,j}(:,2) = pSigma(:,j);
    corrup(1,j)=corr(uSigma(:,j),pSigma(:,j));
end

density_upProbe = cell(1,c);       % returns an n-by-n cell array of empty matrices for initialization
for j=4:c
   
    figure(j);
    subplot(5,3,13)
    hold on;
    
    %call the routine
    [bandwidth,density,X,Y]=mjKDE2d(upProbe{1,j});
    
    %plot the data and the density estimate
    contour3(X,Y,density,10), hold on
    %plot(upProbe{1,j}(:,1),upProbe{1,j}(:,2),'r.','MarkerSize',5)
    colorbar;
    xlabel('uSigma'); ylabel('pSigma');
    title(['up JPDF for probe:',num2str(j-2)])
    
    density_upProbe{1,j} = density;
end

%JPDF for v'p'

vpProbe = cell(1,c);             % returns an n-by-n cell array of empty matrices for initialization
corrvp = zeros(1,c);
for j=4:c
    vpProbe{1,j}(:,1) = vSigma(:,j);
    vpProbe{1,j}(:,2) = pSigma(:,j);
    corrvp(1,j)=corr(vSigma(:,j),pSigma(:,j));
end

density_vpProbe = cell(1,c);       % returns an n-by-n cell array of empty matrices for initialization
for j=4:c
   
    figure(j);
    subplot(5,3,14)
    hold on;
    
    %call the routine
    [bandwidth,density,X,Y]=mjKDE2d(vpProbe{1,j});
    
    %plot the data and the density estimate
    contour3(X,Y,density,10), hold on
    %plot(vpProbe{1,j}(:,1),vpProbe{1,j}(:,2),'r.','MarkerSize',5)
    colorbar;
    xlabel('vSigma'); ylabel('pSigma');
    title(['vp JPDF for probe:',num2str(j-2)])
    
    density_vpProbe{1,j} = density;
end

%JPDF for T'p'

TpProbe = cell(1,c);             % returns an n-by-n cell array of empty matrices for initialization
corrTp = zeros(1,c);
for j=4:c
    TpProbe{1,j}(:,1) = TSigma(:,j);
    TpProbe{1,j}(:,2) = pSigma(:,j);
    corrTp(1,j)=corr(TSigma(:,j),pSigma(:,j));
end

density_TpProbe = cell(1,c);       % returns an n-by-n cell array of empty matrices for initialization
for j=4:c
   
    figure(j);
    subplot(5,3,15)
    hold on;
    
    %call the routine
    [bandwidth,density,X,Y]=mjKDE2d(TpProbe{1,j});
    
    %plot the data and the density estimate
    contour3(X,Y,density,10), hold on
    %plot(TpProbe{1,j}(:,1),TpProbe{1,j}(:,2),'r.','MarkerSize',5)
    colorbar;
    xlabel('TSigma'); ylabel('pSigma');
    title(['Tp JPDF for probe:',num2str(j-2)])
    
    density_TpProbe{1,j} = density;
end