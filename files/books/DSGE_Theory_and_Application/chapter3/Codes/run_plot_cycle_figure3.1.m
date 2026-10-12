%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 根据解析解模拟经济变量走势            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 作者: Deng, Guichuan                  %
% 邮箱: denggch@mail.sysu.edu.cn        %
% 最新更新时间: 2026/10/11              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;

% 参数值
beta   =  0.99;
sigma  =  2.10;
varphi =  3.00;
epsilon=  10.0;
rho    =  0.88;
sigmaa =  0.01;

% 模拟
T = 21;
etaa = normrnd(0, sigmaa, [T-1, 1]);
a = zeros(T, 1); w = zeros(T, 1); y = zeros(T, 1); n  = zeros(T, 1); r  = zeros(T, 1);
for tt = 2:T;
    a(tt) = rho*a(tt-1) + etaa(tt-1);
    w(tt) = a(tt);
    y(tt) = (1+varphi)/(sigma+varphi)*a(tt);
    n(tt) = (1-sigma)/(sigma+varphi)*a(tt);
    r(tt) = sigma*(1+varphi)*(rho-1)/(sigma+varphi)*a(tt);
end

% 时间序列图
t = 1:T-1; fsize = 12
h0 = figure 
set(h0, 'Units','centimeters','Position', [2 2 14 8]);
subplot(2, 1, 1)
plot(t, etaa*100, 'Color', [0,0,0], 'LineStyle', '-', 'LineWidth', 0.5, 'marker', '^')
xlim([1, T-1]);
ylabel('$\eta_t^a\times100$', 'interpreter', 'latex', 'Fontsize', fsize);

subplot(2, 1, 2)
plot(t, y(t+1)*100, 'Color', [0,0,0], 'LineStyle', '-', 'LineWidth', 1.5)
hold on;
plot(t, w(t+1)*100, 'Color', [0.2,0.2,0.2], 'LineStyle', '--', 'LineWidth', 1.5);
plot(t, n(t+1)*100, 'Color', [0.5,0.5,0.5], 'LineStyle', ':', 'LineWidth', 1.5);
plot(t, r(t+1)*100, 'Color', [0.7,0.7,0.7], 'LineStyle', '-.', 'LineWidth', 1.5)
hold off;
xlim([1, T-1]);
ylabel('%');
legend({'$y_t$', '$\hat{w}_t$', '$n_t$', '$\hat{r}_t$'}, 'interpreter', 'latex', 'Fontsize', fsize, 'position', [0.45, 0.03, 0.10, 0.25], 'Orientation','horizontal');
legend('boxoff')




