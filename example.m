
%
% This script attemps to reproduce the figures in the paper
% Sugihara et al. Science 2012
% 20120408 - Leonardo Barbosa
%

%% define dynamical system and generate samples

update = @(v, u, r, beta)(v * (r - r*v - beta.*u));

r_x = 3.8; beta_x_y = 0.02;
r_y = 3.5; beta_y_x = 0.1;

D = 0;
T = 1000;

% D = 5000;
% T = 8500;

% D = 0;
% T = 3500;

% D = 10;
% T = 3510;

x = zeros(1,T);
y = zeros(1,T);
x(1) = 0.4;
y(1) = 0.2;

for t = 1:T-1
    x(t+1) = update(x(t), y(t), r_x, beta_x_y);
    y(t+1) = update(y(t), x(t), r_y, beta_y_x);
end

if D
    x = x(D+1:end);
    y = y(D+1:end);
end

[rho, p] = corr(x',y');

%% quickly inspect the dynamics of the system

figure
plot(y, 'r', 'linewidth', 2);
hold on
plot(x, 'b', 'linewidth', 2);
ylim([.2 1])
xlim([100 200])
title(sprintf('Correlation: %1.4f (%1.3f)', rho, p));
% hold on
% pause(.1)


% 
% figure
% subplot(1,3,1)
% plot(x, y, 'k.');
% 
% subplot(1,3,2)
% plot(x, 'r.');
% 
% subplot(1,3,3)
% plot(y, 'b.');
% 
% % hold on
% % pause(.1)
% 
% 

%% caculate CCM

fprintf('\nCalculating CCM...')

tau = 1;
E = 2;
t0 = 1+(E-1)*tau;

% step = 1;
% Ls = 1:step:T;

step = 10;
Ls = 30:step:T-D;

rhos_x = zeros(1, length(Ls));
rhos_y = zeros(1, length(Ls));

for iL = 1:length(Ls)
% parfor iL = 1:length(Ls)

    L = Ls(iL);
    fprintf('finding ccm for L = %d (%d of %d)\n', L, iL, length(Ls))

    % This did not work, otherwise it is too noisy (very few samples for small L)
%     points = 1:L;
%     [x_R, y_R] = ccm(x(points), y(points), tau, E);
%     rhos_x(iL) = corr(x(points(t0:end))',x_R');
%     rhos_y(iL) = corr(y(points(t0:end))',y_R');

    % This reproduces Figure 1A, so I'm guessing that is what they did
    intervals = L:step:Ls(end);
    n_corr = length(intervals);
    x_cors = zeros(1, n_corr);
    y_cors = zeros(1, n_corr);
    for iiL = 1:n_corr
        l = intervals(iiL);
        points = 1+l-L:l;
        [x_R, y_R] = ccm(x(points), y(points), tau, E);
        x_cors(iiL) = corr(x(points(t0:end))',x_R');
        y_cors(iiL) = corr(y(points(t0:end))',y_R');
    end
    rhos_x(iL) = mean(x_cors);
    rhos_y(iL) = mean(y_cors);
end

fprintf('done.\n\n')

%% Reproduce Figure 1 A

figure;
plot(Ls, rhos_x, 'b', 'linewidth', 2);
hold on
plot(Ls, rhos_y, 'r', 'linewidth', 2);
ylim([0 1])
title('Figure 1A');


%% Reproduce Figure 1 C and D

fprintf('\n Reproducing figures 1 C and D...')

T = 3500;

% for some reason the figure 1C looks very different, with an empty region
beta_x_y = 0; 

x1cd = zeros(1,T); y1cd = zeros(1,T);
x1cd(1) = 0.4; y1cd(1) = 0.2;
for t = 1:T-1
    x1cd(t+1) = update(x1cd(t), y1cd(t), r_x, beta_x_y);
    y1cd(t+1) = update(y1cd(t), x1cd(t), r_y, beta_y_x);
end

tau = 1; E = 2;
t0 = 1+(E-1)*tau;
[x1cd_R, y1cd_R, x1cd_o, y1cd_o] = ccm(x1cd, y1cd, tau, E);

% , 'markersize', 10
[r, p] = corr(y1cd(t0:end)', y1cd_R');
figure; plot(y1cd(t0:end), y1cd_R, 'r.'); ylim([0 1]); xlim([ 0 1]); title(sprintf('Figure 1C (%1.4f [%1.4f])', r, p));
[r, p] = corr(x1cd(t0:end)', x1cd_R');
figure; plot(x1cd(t0:end), x1cd_R, 'b.'); ylim([0 1]); xlim([ 0 1]); title(sprintf('Figure 1D (%1.4f [%1.4f])', r, p));

fprintf('done.\n\n')


