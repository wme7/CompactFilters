%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Testing explicit (Finite-Difference based) filters with holes
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
nx=31; xa=0; xb=1; x=linspace(xa,xb,nx); dx=1.0/(nx-1);
ny=31; ya=0; yb=1; y=linspace(ya,yb,ny); dy=1.0/(ny-1);
[X,Y] = meshgrid(x,y);

% Build filter's mask
mask = not((X<=0.5)&(Y>=0.5));
mask_holes = not(mask);

% Build Matlab's logo
M = 100 + 160*membrane(1,(nx-1)/2);
M(mask_holes) = NaN; % Hack! Do not touch this cells!

% Load CF to build mask and coefs to implement our filter explicity:
CF = compactFiltersWithHoles('TaylorFilter',[nx,ny],mask_holes,order);
C = compactFiltersWithHoles.TaylorFilterCoefs(order);
C1 =-C(1); C2 =-C(2); C3 =1.0-C(3); % only key coeficients
CFmask = not(CF.Mask);

% Visualize original M-array
subplot(221); plot_me(X,Y,M); view(3); title('Original $f(x,y)$');
xlabel('x'); ylabel('y');

% Add noise to surface
Amp = 20*pi; 
%M_noisy = M + Amp*(cos(pi*X/dx(1)) ); % x-direction
%M_noisy = M + Amp*(cos(pi*Y/dy(1)) ); % y-direction
M_noisy = M + Amp*(cos(pi*X/dx(1)) + cos(pi*Y/dy(1))); % Both directions

% Visualize original M-array
subplot(222); plot_me(X,Y,M_noisy); view(3); title('Noised $f(x,y)$');

% Apply filter explicitly
M_temporal = zeros(size(M_noisy));
M_filtered = zeros(size(M_noisy));
tic;
% filter in x-direction
for j=1:ny
    for i=1:nx
        if CFmask(j,i)
            M_temporal(j,i) = ... 
                C1*M_noisy( j ,i-2) + ...
                C2*M_noisy( j ,i-1) + ...
                C3*M_noisy( j , i ) + ...
                C2*M_noisy( j ,i+1) + ...
                C1*M_noisy( j ,i+2) ;
        else
            M_temporal(j,i) = M_noisy(j,i);
        end
    end
end
% filter in y-direction
for i=1:nx
    for j=1:ny
        if CFmask(j,i)
            M_filtered(j,i) = ...
                C1*M_temporal(j-2, i ) + ...
                C2*M_temporal(j-1, i ) + ...
                C3*M_temporal( j , i ) + ...
                C2*M_temporal(j+1, i ) + ...
                C1*M_temporal(j+2, i ) ;
        else
            M_filtered(j,i) = M_temporal(j,i);
        end
    end
end
toc;

% Visualize original M-array
subplot(223); plot_me(X,Y,M_filtered); view(3); title('Filtered $f(x,y)$');

% Compute difference
Err_norm = (abs(M-M_filtered).^2)./(abs(M).^2)*100;
Err_norm(not(CFmask)) = NaN; % Hack! Do not count this cells!

subplot(224); plot_me(X,Y,Err_norm); view(3); title('Error [\%]');
print('figures/Test_explicitFilters2d','-dpng');

%% Visualization tool
function plot_me(x,y,data)
    surface(x,y,data,'EdgeColor','none');
    colormap cool

    l1 = light;
    l1.Position = [160 400 80];
    l1.Style = 'local';
    l1.Color = [0 0.8 0.8];

    l2 = light;
    l2.Position = [.5 -1 .4];
    l2.Color = [0.8 0.8 0];
end