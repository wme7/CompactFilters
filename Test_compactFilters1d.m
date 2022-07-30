%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Testing compact (Finite-Difference based) filters with holes
%         by Manuel A. Diaz, Institut Pprime | Univ-Poitiers 2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, close all;

% Setup figure
h = figure(1);
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultTextFontName','Times',...
'DefaultTextFontSize',20,...
'DefaultAxesFontName','Times',...
'DefaultAxesFontSize',20,...
'DefaultLineLineWidth',1.5,...
'DefaultAxesBox','on',...
'defaultAxesLineWidth',1.5,...
'DefaultFigureColor','w',...
'DefaultLineMarkerSize',7.75)

%% Controling parameters
order = 4; % filter polynomial order
P=order/2; % stencils points [-p,..,-1,0,1,...p]

% Build a discrete mesh
nx=121; xa=-1; xb=1; x=linspace(xa-0.1,xb+0.1,nx); dx=1.0/(nx-1);

% Build filter's mask
mask = (x>=xa)&(x<=xb);
mask_holes = not(mask);

% Build Matlab's logo
M = pi + cos(3*pi*x/2);
M(mask_holes) = NaN; % Hack! Do not touch this cells!

% Load a high-pass filter
CF = compactFiltersWithHoles('TaylorFilter',(nx),mask_holes,order);
Fx = CF.filter_x; % Low-pass filter
CFmask = not(CF.Mask);

% Add noise to surface
Amp = 1*pi; 
M_noisy = M + Amp*cos(15*pi*x/dx(1)); % grid-to-grid noise

% Apply filter
tic; M_filtered = Fx*M_noisy(:); toc

% Visualize original M-array
f1=figure(1);
plot(x,M,'-k',x,M_noisy,'-c',x,M_filtered,'-.m'); 
xlabel('$x$'); hold off; ylabel('$f(x)$');
legend({'Original $f(x)$','Noised $f(x)$','Filtered $f(x)$'},'location','best','orientation','horizontal');
legend boxoff;
print(['figures/Test_',CF.name,'1d'],'-dpng');

% Compute difference
Err_norm = (abs(M(:)-M_filtered).^2)./(abs(M(:)).^2)*100;

f2=figure(2); plot(x,Err_norm,'-m'); 
xlabel('$x$'); ylabel('Error [\%]');
legend({'Filtered $f(x)$ Error [\%]'},'location','best');
print(f2,['figures/Test_',CF.name,'1d_error'],'-dpng');
