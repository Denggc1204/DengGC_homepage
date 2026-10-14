%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 第6章：引入工资黏性的DSGE模型求解        
% 两种情形：                       
% （1）无指数化调整；              
% （2）有指数化调整；             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 作者: Deng, Guichuan             
% 邮箱: denggch@mail.sysu.edu.cn   
% 最新更新时间: 2026/10/15         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 解析解
clear; 
close all; 
clc;

% 脉冲反应期数
T = 30;  

% 初始化
IRFpiw_flex = zeros([1, T]); IRFw_flex = zeros([1, T]); 
IRFpiw_stic = zeros([1, T]); IRFw_stic = zeros([1, T]); 
IRFpi_flex = zeros([1, T]);  IRFp_flex = zeros([1, T]); 
IRFpi_stic = zeros([1, T]);  IRFp_stic = zeros([1, T]); 

% 参数值
bet  = 0.99;  sig  = 2.10;  vph  = 3.00;
epw  = 10.0;  thtw = 0.75;  phipi= 1.50;  
rho  = 0.88;  siga = 0.01;  
kapw = (1-thtw)*(1-bet*thtw)/thtw;
kap1 = kapw*(1+vph)/(1+epw*vph) + 1 + bet*(1-rho);

% 通胀政策函数系数
tht2 = sig*(1+epw*vph)/(kapw*(sig+vph)*phipi + sig*(1+epw*vph));
tht1 = -(kapw*(sig + vph)/(sig*(1+epw*vph))*(phipi - rho) + (1 - bet*rho)*(1-rho))^(-1) *...
    (kapw*(sig + vph)/(sig*(1+epw*vph))*(phipi-1) + (kap1/tht2 - bet)*(1-rho))*tht2;

% 脉冲反应函数
IRFpi_flex(1) = -siga; IRFpi_flex(2) = -siga*rho; 
IRFp_flex(1) = IRFpi_flex(1); IRFp_flex(2) = IRFp_flex(1) + IRFpi_flex(2); 
IRFpi_stic(1) = tht1*siga; 
IRFpi_stic(2) = tht1*rho*siga+tht2*siga;
IRFp_stic(1) = IRFpi_stic(1); 
IRFp_stic(2) = IRFp_stic(1)+IRFpi_stic(2);

IRFpiw_flex(2) = -siga; 
IRFw_flex(2) = IRFw_flex(1) + IRFpiw_flex(2);
IRFpiw_stic(1) = siga+IRFpi_stic(1); 
IRFpiw_stic(2) = (rho-1)*siga + IRFpi_stic(2);
IRFw_stic(1) = IRFpiw_stic(1); 
IRFw_stic(2) = IRFw_stic(1)+IRFpiw_stic(2);

for ii = 3:T
   IRFpi_flex(ii) = IRFpi_flex(ii-1)*rho;
   IRFp_flex(ii) = IRFp_flex(ii-1)+IRFpi_flex(ii);
   IRFpi_stic(ii) = IRFpi_stic(ii-1)*rho; 
   IRFp_stic(ii) = IRFp_stic(ii-1)+IRFpi_stic(ii); 
   IRFpiw_flex(ii) = IRFpiw_flex(ii-1)*rho; 
   IRFw_flex(ii) = IRFw_flex(ii-1)+IRFpiw_flex(ii); 
   IRFpiw_stic(ii) = (rho-1)*rho^(ii-1)*siga + IRFpi_stic(ii);
   IRFw_stic(ii) = IRFw_stic(ii-1)+IRFpiw_stic(ii);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% 图6.2：脉冲反应函数图
%%%%%%%%%%%%%%%%%%%%%%%%%
X = 10; x = 1:X;
h0 = figure;
set(h0, 'unit', 'centimeters', 'position', [2 2 14 10]);
subplot(2, 2, 1)
plot(x, IRFw_stic(x)*100, 'Color', [0.3, 0.3, 0.3], 'LineStyle', '-', 'LineWidth', 1.5);
xlim([1, X]);
hold on;
plot(x, IRFw_flex(x)*100, 'Color', [0, 0, 0], 'LineStyle', '--', 'LineWidth', 1.5);
hold off;
ylabel('(%)', 'interpreter', 'latex', 'Fontsize', 12);
title('名义工资', 'Fontname', 'kaishu', 'fontsize', 12)
xlabel('(a)', 'fontsize', 12)
legend('$w_t$', '$\bar{w}_t$', 'interpreter', 'latex', ...
    'position', [0.30, 0.75, 0.10, 0.15], 'Fontsize', 12);
legend('boxoff')

subplot(2,2,2)
plot(x, IRFp_stic(x)*100, 'Color', [0.3, 0.3, 0.3], 'LineStyle', '-', 'LineWidth', 1.5);
xlim([1, X]);
hold on;
plot(x, IRFp_flex(x)*100, 'Color', [0, 0, 0], 'LineStyle', '--', 'LineWidth', 1.5);
hold off;
ylabel('(%)', 'interpreter', 'latex', 'Fontsize', 12);
legend('$p_t$', '$\bar{p}_t$', 'interpreter', 'latex', ...
    'position', [0.75, 0.72, 0.10, 0.15], 'Fontsize', 12);
legend('boxoff')
title('名义价格', 'Fontname', 'kaishu', 'fontsize', 12)
xlabel('(b)', 'fontsize', 12)

subplot(2, 2, 3)
plot(x, IRFpiw_stic(x)*100, 'Color', [0.3, 0.3, 0.3], 'LineStyle', '-', 'LineWidth', 1.5);
xlim([1, X]);
hold on;
plot(x, IRFpiw_flex(x)*100, 'Color', [0, 0, 0], 'LineStyle', '--', 'LineWidth', 1.5);
hold off;
ylabel('(%)', 'interpreter', 'latex', 'Fontsize', 12);
title('工资通胀', 'Fontname', 'kaishu', 'fontsize', 12)
xlabel('(c)', 'fontsize', 12)
legend('$\pi_t^w$', '$\bar{\pi}_t^w$', 'interpreter', 'latex', ...
    'position', [0.22, 0.27, 0.10, 0.15], 'Fontsize', 12);
legend('boxoff')

subplot(2, 2, 4)
plot(x, IRFpi_stic(x)*100, 'Color', [0.3, 0.3, 0.3], 'LineStyle', '-', 'LineWidth', 1.5);
xlim([1, X]);
hold on;
plot(x, IRFpi_flex(x)*100, 'Color', [0, 0, 0], 'LineStyle', '--', 'LineWidth', 1.5);
plot([0, X], [0, 0], 'Color', [1, 0, 0], 'LineStyle', '--', 'LineWidth', 0.75);
hold off;
ylabel('(%)', 'interpreter', 'latex', 'Fontsize', 12);
legend('$\pi_t$', '$\bar{\pi}_t$', 'interpreter', 'latex', ...
    'position', [0.66, 0.27, 0.10, 0.15], 'Fontsize', 12);
legend('boxoff')
title('价格通胀', 'Fontname', 'kaishu', 'fontsize', 12)
xlabel('(d)', 'fontsize', 12)


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 计算不同黏性程度下模型解和脉冲反应
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; 
close all;
clc; 

% 参数值
bet  = 0.99;  sig  = 2.10;  vph  = 3.00;
epw  = 10.0;  phipi= 1.50;  rho  = 0.88;  
thtw = 0.01:0.01:0.99;

% 通胀和产出缺口政策函数值
for ii = 1:length(thtw)
    kapw(ii) = (1-thtw(ii))*(1-bet*thtw(ii))/thtw(ii);
    kap1(ii) = kapw(ii)*(1+vph)/(1+epw*vph)+1+bet*(1-rho);
    tht2(ii) = sig*(1+epw*vph)/(kapw(ii)*(sig+vph)*phipi + sig*(1+epw*vph));
    tht1(ii) = - (kapw(ii)*(sig + vph)/(sig*(1+epw*vph))*(phipi - rho) + (1-bet*rho)*(1-rho))^(-1)*...
        tht2(ii)*((kap1(ii)/tht2(ii)-bet)*(1-rho) + kapw(ii)*(sig + vph)/(sig*(1+epw*vph))*(phipi-1));
    tht1y(ii) = ((1 - bet*rho)*(1+tht1(ii)) + bet*(1-tht2(ii)))*(1+epw*vph)/(kapw(ii)*(sig+vph));
    tht2y(ii) = (1+epw*vph)*(1-tht2(ii))/(kapw(ii)*(sig+vph));
end

% 脉冲反应函数
T = 30; IRFpi = zeros([T, 2]); IRFpiw = zeros([T, 2]);  
IRFygap = zeros([T, 2]); siga = 0.01; 
clear thtw kapw kap1; 
thtw = [0.25, 0.75];  
kapw = (1-thtw).*(1-bet*thtw)./thtw;
kap1 = kapw*(1+vph)/(1+epw*vph)+1+bet*(1-rho);
tht2s = sig*(1+epw*vph)./(kapw*(sig+vph)*phipi + sig*(1+epw*vph));
tht1s = -(kapw*(sig+vph)/(sig*(1+epw*vph))*(phipi-rho) + (1-bet*rho)*(1-rho)).^(-1).*...
        tht2s.*((kap1./tht2s-bet)*(1-rho) + kapw*(sig+vph)/(sig*(1+epw*vph))*(phipi-1));
IRFpi(1, :) = tht1s*siga; IRFpi(2, :) = IRFpi(1, :)*rho+tht2s*siga;
IRFpiw(1, :) = IRFpi(1, :)+siga; IRFpiw(2, :) = IRFpi(2, :)+(rho-1)*siga;
IRFygap(1,:) = ((1 - bet*rho).*(1+tht1s) + bet*(1-tht2s))*(1+epw*vph)./(kapw*(sig+vph))*siga;
IRFygap(2,:) = IRFygap(1,:)*rho - (1+epw*vph)*(1-tht2s)./(kapw*(sig+vph))*siga;
for ii = 3:T
    IRFpi(ii, :) = IRFpi(ii-1,:)*rho; 
    IRFpiw(ii, :) = IRFpi(ii, :) + (rho-1)*rho^(ii-2)*siga; 
    IRFygap(ii, :) = IRFygap(ii-1, :)*rho;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 图6.3： 工资黏性程度与脉冲反应
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = 15; x = 1:X;
h0 = figure;
set(h0, 'unit', 'centimeters', 'position', [2 2 14 10]);
subplot(2, 2, 1)
plot(x, IRFpi(x, 1)*100, 'Color', [0, 0, 0], 'LineStyle', '--', 'LineWidth', 1.5);
hold on;
plot(x, IRFpi(x, 2)*100, 'Color', [0, 0, 0], 'LineStyle', '-', 'LineWidth', 1.5);
plot([0, X], [0, 0], 'Color', [1, 0, 0], 'LineStyle', '--', 'LineWidth', 0.75);
hold off;
xlim([1, X]);
ylabel('(%)');
xlabel('(a)', 'fontsize', 12)
title('$\pi_t$', 'interpreter', 'latex', 'Fontsize', 12);
legend('$\theta_w=0.25$', '$\theta_w=0.75$', 'interpreter', 'latex', 'orientation', ...
    'vertical', 'position', [0.25, 0.65, 0.10, 0.15], 'Fontsize', 12);
legend('boxoff')

subplot(2, 2, 2)
plot(x, IRFygap(x, 1)*100, 'Color', [0, 0, 0], 'LineStyle', '--', 'LineWidth', 1.5);
hold on;
plot(x, IRFygap(x, 2)*100, 'Color', [0, 0, 0], 'LineStyle', '-', 'LineWidth', 1.5);
plot([0, X], [0, 0], 'Color', [1, 0, 0], 'LineStyle', '--', 'LineWidth', 0.75);
hold off;
xlim([1, X]);
ylabel('(%)');
xlabel('(b)', 'fontsize', 12)
title('$\tilde{y}_t$', 'interpreter', 'latex', 'Fontsize', 12);

subplot(2, 2, 3)
plot(0.01:0.01:0.99, tht1, 'Color', [0, 0, 0], 'LineStyle', '-', 'LineWidth', 1.5);
xlim([0, 1]);
ylabel('$\theta_1^\pi$', 'interpreter', 'latex', 'Fontsize', 12);
xlabel('$(c)\theta_w$', 'interpreter', 'latex', 'Fontsize', 12)
title('$\theta_1^\pi$', 'interpreter', 'latex', 'Fontsize', 12);

subplot(2, 2, 4)
plot(0.01:0.01:0.99, tht1y, 'Color', [0, 0, 0], 'LineStyle', '-', 'LineWidth', 1.5);
hold on;
plot([0, 1], [0, 0], 'Color', [1, 0, 0], 'LineStyle', '--', 'LineWidth', 0.75);
hold off;
xlim([0, 1]);
ylabel('$\theta_1^{\tilde{y}}$', 'interpreter', 'latex', 'Fontsize', 12);
xlabel('$(d)\theta_w$', 'interpreter', 'latex', 'Fontsize', 12)
title('$\theta_1^{\tilde{y}}$', 'interpreter', 'latex', 'Fontsize', 12);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 货币政策强度与脉冲反应
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; 
close all;
clc; 

% 参数值
bet  = 0.99;  sig  = 2.10;  vph  = 3.00;
epw  = 10.0;  thtw = 0.75;  rho  = 0.88; 
kapw = (1-thtw)*(1-bet*thtw)/thtw;
kap1 = kapw*(1+vph)/(1+epw*vph)+1+bet*(1-rho);

% 不同货币政策强度下政策函数系数
phipi= 1.01:0.01:1000;
for ii = 1:length(phipi)
    tht2(ii) = sig*(1+epw*vph)/(kapw*(sig+vph)*phipi(ii) + sig*(1+epw*vph));
    tht1(ii) = - (kapw*(sig + vph)/(sig*(1+epw*vph))*(phipi(ii) - rho) + (1-bet*rho)*(1-rho))^(-1)*...
        tht2(ii)*((kap1/tht2(ii)-bet)*(1-rho) + kapw*(sig + vph)/(sig*(1+epw*vph))*(phipi(ii)-1));
    tht1y(ii) = ((1 - bet*rho)*(1+tht1(ii)) + bet*(1-tht2(ii)))*(1+epw*vph)/(kapw*(sig+vph));
    tht2y(ii) = (1+epw*vph)*(1-tht2(ii))/(kapw*(sig+vph));
end

% 三种货币政策强度与脉冲反应
clear; 
close all; 
clc;

% 脉冲反应初始化
T = 30; IRFpi = zeros([T, 3]); IRFpiw = zeros([T, 3]);  
IRFygap = zeros([T, 3]);  IRFy = zeros([T, 3]); 

% 参数值
bet  = 0.99;  sig  = 2.10;  vph  = 3.00;
epw  = 10.0;  thtw = 0.75;  rho  = 0.88; siga = 0.01; 
phipi = [1.5, 5.0, 8.5];  
kapw = (1-thtw)*(1-bet*thtw)/thtw;
kap1 = kapw*(1+vph)/(1+epw*vph)+1+bet*(1-rho);

% 通胀政策函数值
tht2 = sig*(1+epw*vph)./(kapw*(sig+vph)*phipi + sig*(1+epw*vph));
tht1 = -(kapw*(sig+vph)/(sig*(1+epw*vph))*(phipi-rho) + (1-bet*rho)*(1-rho)).^(-1).*...
        tht2.*((kap1./tht2-bet)*(1-rho) + kapw*(sig+vph)/(sig*(1+epw*vph))*(phipi-1));

% 脉冲反应函数值
IRFpi(1, :) = tht1*siga; IRFpi(2, :) = IRFpi(1, :)*rho+tht2*siga;
IRFpiw(1, :) = IRFpi(1, :)+siga; IRFpiw(2, :) = IRFpi(2, :)+(rho-1)*siga;
IRFygap(1,:) = ((1 - bet*rho).*(1+tht1) + bet*(1-tht2))*(1+epw*vph)./(kapw*(sig+vph))*siga;
IRFygap(2,:) = IRFygap(1,:)*rho - (1+epw*vph)*(1-tht2)./(kapw*(sig+vph))*siga;
IRFy(1,:) = IRFygap(1,:) + (1+vph)/(sig+vph)*siga;
IRFy(2,:) = IRFygap(2,:) + (1+vph)/(sig+vph)*siga*rho;
for ii = 3:T
    IRFpi(ii, :) = IRFpi(ii-1,:)*rho; 
    IRFpiw(ii, :) = IRFpi(ii, :) + (rho-1)*rho^(ii-2)*siga; 
    IRFygap(ii, :) = IRFygap(ii-1, :)*rho;
    IRFy(ii, :) = IRFygap(ii, :) + (1+vph)/(sig+vph)*rho^(ii-1)*siga;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 图6.4：不同政策函数强度与脉冲反应
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = 15; x = 1:X;
h0 = figure;
set(h0, 'unit', 'centimeters', 'position', [2 2 14 10]);
subplot(2, 2, 1)
plot(x, IRFpiw(x, 1)*100, 'Color', [0, 0, 0], 'LineStyle', '--', 'LineWidth', 1.5);
hold on;
plot(x, IRFpiw(x, 2)*100, 'Color', [0, 0, 0], 'LineStyle', ':', 'LineWidth', 1.5);
plot(x, IRFpiw(x, 3)*100, 'Color', [0, 0, 0], 'LineStyle', '-', 'LineWidth', 1.5);
plot([0, X], [0, 0], 'Color', [1, 0, 0], 'LineStyle', '--', 'LineWidth', 0.75);
hold off;
xlim([1, X]);
ylabel('(%)');
title('$\pi_t^w$', 'interpreter', 'latex', 'Fontsize', 12);

subplot(2, 2, 2)
plot(x, IRFygap(x, 1)*100, 'Color', [0, 0, 0], 'LineStyle', '--', 'LineWidth', 1.5);
hold on;
plot(x, IRFygap(x, 2)*100, 'Color', [0, 0, 0], 'LineStyle', ':', 'LineWidth', 1.5);
plot(x, IRFygap(x, 3)*100, 'Color', [0, 0, 0], 'LineStyle', '-', 'LineWidth', 1.5);
plot([0, X], [0, 0], 'Color', [1, 0, 0], 'LineStyle', '--', 'LineWidth', 0.75);
hold off;
xlim([1, X]);
ylabel('(%)');
title('$\tilde{y}_t$', 'interpreter', 'latex', 'Fontsize', 12);

subplot(2, 2, 3)
plot(x, IRFpi(x, 1)*100, 'Color', [0, 0, 0], 'LineStyle', '--', 'LineWidth', 1.5);
hold on;
plot(x, IRFpi(x, 2)*100, 'Color', [0, 0, 0], 'LineStyle', ':', 'LineWidth', 1.5);
plot(x, IRFpi(x, 3)*100, 'Color', [0, 0, 0], 'LineStyle', '-', 'LineWidth', 1.5);
plot([0, X], [0, 0], 'Color', [1, 0, 0], 'LineStyle', '--', 'LineWidth', 0.75);
hold off;
xlim([1, X]);
ylabel('(%)');
title('$\pi_t$', 'interpreter', 'latex', 'Fontsize', 12);
legend('$\phi_\pi=1.5$', '$\phi_\pi=5.0$', '$\phi_\pi=8.5$', 'interpreter', 'latex', 'orientation', ...
    'vertical', 'position', [0.25, 0.15, 0.10, 0.15], 'Fontsize', 12);
legend('boxoff')

subplot(2, 2, 4)
plot(x, IRFy(x, 1)*100, 'Color', [0, 0, 0], 'LineStyle', '--', 'LineWidth', 1.5);
hold on;
plot(x, IRFy(x, 2)*100, 'Color', [0, 0, 0], 'LineStyle', ':', 'LineWidth', 1.5);
plot(x, IRFy(x, 3)*100, 'Color', [0, 0, 0], 'LineStyle', '-', 'LineWidth', 1.5);
plot([0, X], [0, 0], 'Color', [1, 0, 0], 'LineStyle', '--', 'LineWidth', 0.75);
hold off;
xlim([1, X]);
ylabel('(%)');
title('$y_t$', 'interpreter', 'latex', 'Fontsize', 12);


%% 使用数值方法求解工资黏性模型

% 新建子文件夹
mkdir('./result')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 求解模型：
% （1）无指数化调整：高黏性、中黏性、低黏性
% （2）有指数化调整：高黏性、中黏性、低黏性
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
% （1）无指数化调整
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
% （1.1） 高黏性
%%%%%%%%%%%%%%%%%
clear; 
close all; 
clc;

% 调整参数值
thtw = 0.25;

% 求解模型和模拟
dynare wagesticky.mod;

% 保存结果
wo025 = oo_;
save('./result./wo025.mat', 'wo025');

%%%%%%%%%%%%%%%%%
% （1.2） 中黏性
%%%%%%%%%%%%%%%%%
clear; 
close all; 
clc;

% 调整参数值
thtw = 0.50;

% 求解模型和模拟
dynare wagesticky.mod;

% 保存结果
wo050 = oo_;
save('./result./wo050.mat', 'wo050');

%%%%%%%%%%%%%%%%%
% （1.3） 低黏性
%%%%%%%%%%%%%%%%%
clear; 
close all; 
clc;

% 调整参数值
thtw = 0.75;

% 求解模型和模拟
dynare wagesticky.mod;

% 保存结果
wo075 = oo_; wm = M_;
save('./result./wo075.mat', 'wo075');
save('./result./wm.mat', 'wm');


%%%%%%%%%%%%%%%%%%%%%%%%%
% （2）有指数化调整
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
% （2.1） 高黏性
%%%%%%%%%%%%%%%%%
clear; 
close all; 
clc;

% 调整参数值
thtw = 0.25;

% 求解模型和模拟
dynare wagesticky.mod -Dindex=1;

% 保存结果
wio025 = oo_;
save('./result./wio025.mat', 'wio025');

%%%%%%%%%%%%%%%%%
% （2.2） 中黏性
%%%%%%%%%%%%%%%%%
clear; 
close all; 
clc;

% 调整参数值
thtw = 0.50;

% 求解模型和模拟
dynare wagesticky -Dindex=1;

% 保存结果
wio050 = oo_;
save('./result./wio050.mat', 'wio050');

%%%%%%%%%%%%%%%%%
% （2.3） 低黏性
%%%%%%%%%%%%%%%%%
clear; 
close all; 
clc;

% 调整参数值
thtw = 0.75;

% 求解模型和模拟
dynare wagesticky -Dindex=1;

% 保存结果
wio075 = oo_;
save('./result./wio075.mat', 'wio075');

%% 画脉冲反应函数
%%%%%%%%%%%%%%%%%
% Plot the IRF  %
%%%%%%%%%%%%%%%%%
clear; 
close all; 
clc;

% 加载数据
load('./result./wo075.mat');
load('./result./wo025.mat');
load('./result./wio075.mat');
load('./result./wio050.mat');
load('./result./wio025.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 图6.5：工资黏性下无指数化调整 VS 有指数化调整
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = 20; x = 1:T;
% indexation
h0 = figure;
set(h0, 'unit', 'centimeters', 'position', [1 1 14 12]);
subplot(3, 3, 1)
plot(x, wo075.irfs.y_etaa(x)*100, 'k-', x, wio075.irfs.y_etaa(x)*100, 'k--', 'LineWidth', 1.25);
ylabel('(%)'); 
title('$y_t$', 'interpreter', 'latex', 'FontSize', 12);

subplot(3, 3, 2)
plot(x, wo075.irfs.c_etaa(x)*100, 'k-', x, wio075.irfs.c_etaa(x)*100, 'k--', 'LineWidth', 1.25);
ylabel('(%)'); 
title('$c_t$', 'interpreter', 'latex', 'FontSize', 12);

subplot(3, 3, 3)
plot(x, wo075.irfs.n_etaa(x)*100, 'k-', x, wio075.irfs.n_etaa(x)*100, 'k--', 'LineWidth', 1.25);
ylabel('(%)'); 
title('$n_t$', 'interpreter', 'latex', 'FontSize', 12);

subplot(3, 3, 4)
plot(x, wo075.irfs.pi_etaa(x)*100, 'k-', x, wio075.irfs.pi_etaa(x)*100, 'k--', 'LineWidth', 1.25);
ylabel('(%)'); 
title('$\pi_t$', 'interpreter', 'latex', 'FontSize', 12);

subplot(3, 3, 5)
plot(x, wo075.irfs.piw_etaa(x)*100, 'k-', x, wio075.irfs.piw_etaa(x)*100, 'k--', 'LineWidth', 1.25);
ylabel('(%)'); 
title('$\pi_t^w$', 'interpreter', 'latex', 'FontSize', 12);

subplot(3, 3, 6)
plot(x, wo075.irfs.w_etaa(x)*100, 'k-', x, wio075.irfs.w_etaa(x)*100, 'k--', 'LineWidth', 1.25);
ylabel('(%)'); 
title('$w_t$', 'interpreter', 'latex', 'FontSize', 12);

subplot(3, 3, 7)
plot(x, wo075.irfs.r_etaa(x)*100, 'k-', x, wio075.irfs.r_etaa(x)*100, 'k--', 'LineWidth', 1.25);
ylabel('(%)'); 
title('$r_t$', 'interpreter', 'latex', 'FontSize', 12);

subplot(3, 3, 8)
plot(x, wo075.irfs.ygap_etaa(x)*100, 'k-', x, wio075.irfs.ygap_etaa(x)*100, 'k--', 'LineWidth', 1.25);
ylabel('(%)'); title('$\tilde{y}_t$', 'interpreter', 'latex', 'FontSize', 12);
legend({'无指数化调整', '有指数化调整'}, ...
    'position', [0.72 0.10 0.15 0.20], 'FontSize', 10, 'FontName', 'kaiti')
legend('boxoff')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 图6.6：指数化调整下不同工资黏性强度的脉冲反应
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h0 = figure;
set(h0, 'unit', 'centimeters', 'position', [1 1 14 12]);
subplot(3, 3, 1)
plot(x, wio025.irfs.y_etaa(x)*100, 'k:', x, wio050.irfs.y_etaa(x)*100, 'k-', x, wio075.irfs.y_etaa(x)*100, 'k--', 'LineWidth', 1.0);
ylabel('(%)'); 
title('$y_t$', 'interpreter', 'latex', 'FontSize', 12);

subplot(3, 3, 2)
plot(x, wio025.irfs.c_etaa(x)*100, 'k:', x, wio050.irfs.c_etaa(x)*100, 'k-', x, wio075.irfs.c_etaa(x)*100, 'k--', 'LineWidth', 1.0);
ylabel('(%)'); 
title('$c_t$', 'interpreter', 'latex', 'FontSize', 12);

subplot(3, 3, 3)
plot(x, wio025.irfs.n_etaa(x)*100, 'k:', x, wio050.irfs.n_etaa(x)*100, 'k-', x, wio075.irfs.n_etaa(x)*100, 'k--', 'LineWidth', 1.0);
ylabel('(%)'); 
title('$n_t$', 'interpreter', 'latex', 'FontSize', 12);

subplot(3, 3, 4)
plot(x, wio025.irfs.pi_etaa(x)*100, 'k:', x, wio050.irfs.pi_etaa(x)*100, 'k-', x, wio075.irfs.pi_etaa(x)*100, 'k--', 'LineWidth', 1.0);
ylabel('(%)'); 
title('$\pi_t$', 'interpreter', 'latex', 'FontSize', 12);

subplot(3, 3, 5)
plot(x, wio025.irfs.piw_etaa(x)*100, 'k:', x, wio050.irfs.piw_etaa(x)*100, 'k-', x, wio075.irfs.piw_etaa(x)*100, 'k--', 'LineWidth', 1.0);
ylabel('(%)'); 
title('$\pi_t^w$', 'interpreter', 'latex', 'FontSize', 12);

subplot(3, 3, 6)
plot(x, wio025.irfs.w_etaa(x)*100, 'k:', x, wio050.irfs.w_etaa(x)*100, 'k-', x, wio075.irfs.w_etaa(x)*100, 'k--', 'LineWidth', 1.0);
ylabel('(%)'); 
title('$w_t$', 'interpreter', 'latex', 'FontSize', 12);

subplot(3, 3, 7)
plot(x, wio025.irfs.r_etaa(x)*100, 'k:', x, wio050.irfs.r_etaa(x)*100, 'k-', x, wio075.irfs.r_etaa(x)*100, 'k--', 'LineWidth', 1.0);
ylabel('(%)'); 
title('$r_t$', 'interpreter', 'latex', 'FontSize', 12);

subplot(3, 3, 8)
plot(x, wio025.irfs.ygap_etaa(x)*100, 'k:', x, wio050.irfs.ygap_etaa(x)*100, 'k-', x, wio075.irfs.ygap_etaa(x)*100, 'k--', 'LineWidth', 1.0);
ylabel('(%)'); title('$\tilde{y}_t$', 'interpreter', 'latex', 'FontSize', 12);
legend({'$\theta_w = 0.25$', '$\theta_w = 0.50$', '$\theta_w = 0.75$'}, 'interpreter', 'latex', ...
    'position', [0.72 0.10 0.15 0.20], 'FontSize', 10, 'FontName', 'kaiti')
legend('boxoff')



%% 计算矩    
clear; 
close all; 
clc;

% 加载数据
load('./result./wo075.mat');
load('./result./wo050.mat');
load('./result./wo025.mat');
load('./result./wio075.mat');
load('./result./wio050.mat');
load('./result./wio025.mat');
load('./result./wm.mat');

% 提取变量位置
yp     = strmatch('y', wm.endo_names ,'exact');
cp     = strmatch('c', wm.endo_names ,'exact');
pip    = strmatch('pi', wm.endo_names ,'exact');
piwp   = strmatch('piw', wm.endo_names ,'exact');
np     = strmatch('n', wm.endo_names ,'exact');
wp     = strmatch('w', wm.endo_names ,'exact');
rp     = strmatch('r', wm.endo_names ,'exact');
ygp    = strmatch('ygap', wm.endo_names ,'exact');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 表6.5：波动率：无条件标准差
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vol = round(...
      sqrt(...
       [wo025.var(yp, yp) wo025.var(cp, cp) wo025.var(pip, pip) wo025.var(piwp, piwp) ... 
        wo025.var(np, np) wo025.var(wp, wp) wo025.var(rp, rp) wo025.var(ygp, ygp);
        wio025.var(yp, yp) wio025.var(cp, cp) wio025.var(pip, pip) wio025.var(piwp, piwp) ... 
        wio025.var(np, np) wio025.var(wp, wp) wio025.var(rp, rp) wio025.var(ygp, ygp);
        wo050.var(yp, yp) wo050.var(cp, cp) wo050.var(pip, pip) wo050.var(piwp, piwp) ... 
        wo050.var(np, np) wo050.var(wp, wp) wo050.var(rp, rp) wo050.var(ygp, ygp);
        wio050.var(yp, yp) wio050.var(cp, cp) wio050.var(pip, pip) wio050.var(piwp, piwp) ... 
        wio050.var(np, np) wio050.var(wp, wp) wio050.var(rp, rp) wio050.var(ygp, ygp);
        wo075.var(yp, yp) wo075.var(cp, cp) wo075.var(pip, pip) wo075.var(piwp, piwp) ...
        wo075.var(np, np) wo075.var(wp, wp) wo075.var(rp, rp) wo075.var(ygp, ygp);
        wio075.var(yp, yp) wio075.var(cp, cp) wio075.var(pip, pip) wio075.var(piwp, piwp) ...
        wio075.var(np, np) wio075.var(wp, wp) wio075.var(rp, rp) wio075.var(ygp, ygp)])*100, 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 表6.6：与产出的相关系数
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cor = round(...
    [1 wo025.var(yp, cp)/sqrt(wo025.var(yp, yp)*wo025.var(cp, cp)) ...
       wo025.var(yp, pip)/sqrt(wo025.var(yp, yp)*wo025.var(pip, pip)) ...
       wo025.var(yp, piwp)/sqrt(wo025.var(yp, yp)*wo025.var(piwp, piwp)) ... 
       wo025.var(yp, np)/sqrt(wo025.var(yp, yp)*wo025.var(np, np)) ...
       wo025.var(yp, wp)/sqrt(wo025.var(yp, yp)*wo025.var(wp, wp)) ...
       wo025.var(yp, rp)/sqrt(wo025.var(yp, yp)*wo025.var(rp, rp)) ...
       wo025.var(yp, ygp)/sqrt(wo025.var(yp, yp)*wo025.var(ygp, ygp));
     1 wio025.var(yp, cp)/sqrt(wio025.var(yp, yp)*wio025.var(cp, cp)) ...
       wio025.var(yp, pip)/sqrt(wio025.var(yp, yp)*wio025.var(pip, pip)) ...
       wio025.var(yp, piwp)/sqrt(wio025.var(yp, yp)*wio025.var(piwp, piwp)) ... 
       wio025.var(yp, np)/sqrt(wio025.var(yp, yp)*wio025.var(np, np)) ...
       wio025.var(yp, wp)/sqrt(wio025.var(yp, yp)*wio025.var(wp, wp)) ...
       wio025.var(yp, rp)/sqrt(wio025.var(yp, yp)*wio025.var(rp, rp)) ...
       wio025.var(yp, ygp)/sqrt(wio025.var(yp, yp)*wio025.var(ygp, ygp));
     1 wo050.var(yp, cp)/sqrt(wo050.var(yp, yp)*wo050.var(cp, cp)) ...
       wo050.var(yp, pip)/sqrt(wo050.var(yp, yp)*wo050.var(pip, pip)) ...
       wo050.var(yp, piwp)/sqrt(wo050.var(yp, yp)*wo050.var(piwp, piwp)) ... 
       wo050.var(yp, np)/sqrt(wo050.var(yp, yp)*wo050.var(np, np)) ...
       wo050.var(yp, wp)/sqrt(wo050.var(yp, yp)*wo050.var(wp, wp)) ...
       wo050.var(yp, rp)/sqrt(wo050.var(yp, yp)*wo050.var(rp, rp)) ...
       wo050.var(yp, ygp)/sqrt(wo050.var(yp, yp)*wo050.var(ygp, ygp));
     1 wio050.var(yp, cp)/sqrt(wio050.var(yp, yp)*wio050.var(cp, cp)) ...
       wio050.var(yp, pip)/sqrt(wio050.var(yp, yp)*wio050.var(pip, pip)) ...
       wio050.var(yp, piwp)/sqrt(wio050.var(yp, yp)*wio050.var(piwp, piwp)) ... 
       wio050.var(yp, np)/sqrt(wio050.var(yp, yp)*wio050.var(np, np)) ...
       wio050.var(yp, wp)/sqrt(wio050.var(yp, yp)*wio050.var(wp, wp)) ...
       wio050.var(yp, rp)/sqrt(wio050.var(yp, yp)*wio050.var(rp, rp)) ...
       wio050.var(yp, ygp)/sqrt(wio050.var(yp, yp)*wio050.var(ygp, ygp));
     1 wo075.var(yp, cp)/sqrt(wo075.var(yp, yp)*wo075.var(cp, cp)) ...
       wo075.var(yp, pip)/sqrt(wo075.var(yp, yp)*wo075.var(pip, pip)) ...
       wo075.var(yp, piwp)/sqrt(wo075.var(yp, yp)*wo075.var(piwp, piwp)) ... 
       wo075.var(yp, np)/sqrt(wo075.var(yp, yp)*wo075.var(np, np)) ...
       wo075.var(yp, wp)/sqrt(wo075.var(yp, yp)*wo075.var(wp, wp)) ...
       wo075.var(yp, rp)/sqrt(wo075.var(yp, yp)*wo075.var(rp, rp)) ...
       wo075.var(yp, ygp)/sqrt(wo075.var(yp, yp)*wo075.var(ygp, ygp));
     1 wio075.var(yp, cp)/sqrt(wio075.var(yp, yp)*wio075.var(cp, cp)) ...
       wio075.var(yp, pip)/sqrt(wio075.var(yp, yp)*wio075.var(pip, pip)) ...
       wio075.var(yp, piwp)/sqrt(wio075.var(yp, yp)*wio075.var(piwp, piwp)) ... 
       wio075.var(yp, np)/sqrt(wio075.var(yp, yp)*wio075.var(np, np)) ...
       wio075.var(yp, wp)/sqrt(wio075.var(yp, yp)*wio075.var(wp, wp)) ...
       wio075.var(yp, rp)/sqrt(wio075.var(yp, yp)*wio075.var(rp, rp)) ...
       wio075.var(yp, ygp)/sqrt(wio075.var(yp, yp)*wio075.var(ygp, ygp))], 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 表6.7：一阶自相关系数
%%%%%%%%%%%%%%%%%%%%%%%%%%
Autocor = round(...
    [wo025.autocorr{1,1}(yp,yp) wo025.autocorr{1,1}(cp,cp) wo025.autocorr{1,1}(pip,pip) ...
     wo025.autocorr{1,1}(piwp,piwp) wo025.autocorr{1,1}(np,np) wo025.autocorr{1,1}(wp,wp) ...
     wo025.autocorr{1,1}(rp,rp) wo025.autocorr{1,1}(ygp,ygp);
     wio025.autocorr{1,1}(yp,yp) wio025.autocorr{1,1}(cp,cp) wio025.autocorr{1,1}(pip,pip) ...
     wio025.autocorr{1,1}(piwp,piwp) wio025.autocorr{1,1}(np,np) wio025.autocorr{1,1}(wp,wp) ...
     wio025.autocorr{1,1}(rp,rp) wio025.autocorr{1,1}(ygp,ygp);
     wo050.autocorr{1,1}(yp,yp) wo050.autocorr{1,1}(cp,cp) wo050.autocorr{1,1}(pip,pip) ...
     wo050.autocorr{1,1}(piwp,piwp) wo050.autocorr{1,1}(np,np) wo050.autocorr{1,1}(wp,wp) ...
     wo050.autocorr{1,1}(rp,rp) wo050.autocorr{1,1}(ygp,ygp);
     wio050.autocorr{1,1}(yp,yp) wio050.autocorr{1,1}(cp,cp) wio050.autocorr{1,1}(pip,pip) ...
     wio050.autocorr{1,1}(piwp,piwp) wio050.autocorr{1,1}(np,np) wio050.autocorr{1,1}(wp,wp) ...
     wio050.autocorr{1,1}(rp,rp) wio050.autocorr{1,1}(ygp,ygp);
     wo075.autocorr{1,1}(yp,yp) wo075.autocorr{1,1}(cp,cp) wo075.autocorr{1,1}(pip,pip) ...
     wo075.autocorr{1,1}(piwp,piwp) wo075.autocorr{1,1}(np,np) wo075.autocorr{1,1}(wp,wp) ...
     wo075.autocorr{1,1}(rp,rp) wo075.autocorr{1,1}(ygp,ygp)
     wio075.autocorr{1,1}(yp,yp) wio075.autocorr{1,1}(cp,cp) wio075.autocorr{1,1}(pip,pip) ...
     wio075.autocorr{1,1}(piwp,piwp) wio075.autocorr{1,1}(np,np) wio075.autocorr{1,1}(wp,wp) ...
     wio075.autocorr{1,1}(rp,rp) wio075.autocorr{1,1}(ygp,ygp)], 3);

save('./result./Vol.mat', 'Vol');
save('./result./Cor.mat', 'Cor');
save('./result./Autocor.mat', 'Autocor');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 工资黏性对矩条件的影响 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; 
close all; 
clc; 
N = 49; Volw = zeros(N, 16); Corw = zeros(N, 16); Autocorw = zeros(N, 16);
for ii = 1:N
    thtw = ii/50
    dynare wagesticky;            % 工资指数化调整
    wo = oo_;
    clear oo_;
    dynare wagesticky -Dindex=1;  % 无指数调整
    wio = oo_;
    wm = M_;
    yp = strmatch('y', wm.endo_names ,'exact');
    cp = strmatch('c', wm.endo_names ,'exact');
    pip = strmatch('pi', wm.endo_names ,'exact');
    piwp = strmatch('piw', wm.endo_names ,'exact');
    np = strmatch('n', wm.endo_names ,'exact');
    wp = strmatch('w', wm.endo_names ,'exact');
    rp = strmatch('r', wm.endo_names ,'exact');
    ygp = strmatch('ygap', wm.endo_names ,'exact');
    % 标准差
    Volw(ii, :) = round(sqrt(...
       [wo.var(yp, yp) wo.var(cp, cp) wo.var(pip, pip) wo.var(piwp, piwp) ... 
        wo.var(np, np) wo.var(wp, wp) wo.var(rp, rp) wo.var(ygp, ygp) ...
        wio.var(yp, yp) wio.var(cp, cp) wio.var(pip, pip) wio.var(piwp, piwp) ... 
        wio.var(np, np) wio.var(wp, wp) wio.var(rp, rp) wio.var(ygp, ygp)])*100, 3);
    % 与产出相关系数
    Corw(ii, :) = round(...
    [1 wo.var(yp, cp)/sqrt(wo.var(yp, yp)*wo.var(cp, cp)) ...
       wo.var(yp, pip)/sqrt(wo.var(yp, yp)*wo.var(pip, pip)) ...
       wo.var(yp, piwp)/sqrt(wo.var(yp, yp)*wo.var(piwp, piwp)) ... 
       wo.var(yp, np)/sqrt(wo.var(yp, yp)*wo.var(np, np)) ...
       wo.var(yp, wp)/sqrt(wo.var(yp, yp)*wo.var(wp, wp)) ...
       wo.var(yp, rp)/sqrt(wo.var(yp, yp)*wo.var(rp, rp)) ...
       wo.var(yp, ygp)/sqrt(wo.var(yp, yp)*wo.var(ygp, ygp)) ...
     1 wio.var(yp, cp)/sqrt(wio.var(yp, yp)*wio.var(cp, cp)) ...
       wio.var(yp, pip)/sqrt(wio.var(yp, yp)*wio.var(pip, pip)) ...
       wio.var(yp, piwp)/sqrt(wio.var(yp, yp)*wio.var(piwp, piwp)) ... 
       wio.var(yp, np)/sqrt(wio.var(yp, yp)*wio.var(np, np)) ...
       wio.var(yp, wp)/sqrt(wio.var(yp, yp)*wio.var(wp, wp)) ...
       wio.var(yp, rp)/sqrt(wio.var(yp, yp)*wio.var(rp, rp)) ...
       wio.var(yp, ygp)/sqrt(wio.var(yp, yp)*wio.var(ygp, ygp))], 3);
    % 一阶自相关系数
    Autocorw(ii, :) = round(...
    [wo.autocorr{1,1}(yp,yp) wo.autocorr{1,1}(cp,cp) wo.autocorr{1,1}(pip,pip) ...
     wo.autocorr{1,1}(piwp,piwp) wo.autocorr{1,1}(np,np) wo.autocorr{1,1}(wp,wp) ...
     wo.autocorr{1,1}(rp,rp) wo.autocorr{1,1}(ygp,ygp) ...
     wio.autocorr{1,1}(yp,yp) wio.autocorr{1,1}(cp,cp) wio.autocorr{1,1}(pip,pip) ...
     wio.autocorr{1,1}(piwp,piwp) wio.autocorr{1,1}(np,np) wio.autocorr{1,1}(wp,wp) ...
     wio.autocorr{1,1}(rp,rp) wio.autocorr{1,1}(ygp,ygp)], 3);
end

save('./result./Volw.mat', 'Volw');
save('./result./Corw.mat', 'Corw');
save('./result./Autocorw.mat', 'Autocorw');


% 画矩条件随工资黏性的变动轨迹
clear; 
close all; 
clc;

% 加载数据
load('./result./Volw.mat');
load('./result./Corw.mat');
load('./result./Autocorw.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 图6.7：工资黏性程度对矩条件的影响
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h0 = figure;
set(h0, 'unit', 'centimeters', 'position', [1 1 14 8]);

% 标准差
subplot(2, 3, 1)
plot(0.02:0.02:0.98, Volw(:, 3), 'k-', 0.02:0.02:0.98, Volw(:, 11), 'k--', 'Linewidth', 1.25)
title('$std(\pi_t)$', 'interpreter', 'latex', 'Fontsize', 12)
xlim([0, 1]); % xlabel('$\theta_w$', 'interpreter', 'latex', 'Fontsize', 12);

subplot(2, 3, 4)
plot(0.02:0.02:0.98, Volw(:, 8), 'k-', 0.02:0.02:0.98, Volw(:, 16), 'k--', 'Linewidth', 1.25);
title('$std(\tilde{y}_t)$', 'interpreter', 'latex', 'Fontsize', 12)
xlim([0, 1]); % xlabel('$\theta_w$', 'interpreter', 'latex', 'Fontsize', 12);

% 与产出相关系数
subplot(2, 3, 2)
plot(0.02:0.02:0.98, Corw(:, 3), 'k-', 0.02:0.02:0.98, Corw(:, 11), 'k--', 'Linewidth', 1.25)
title('$corr(y_t,\pi_t)$', 'interpreter', 'latex', 'Fontsize', 12)
xlim([0, 1]); % xlabel('$\theta_w$', 'interpreter', 'latex', 'Fontsize', 12);

subplot(2, 3, 5)
plot(0.02:0.02:0.98, Corw(:, 8), 'k-', 0.02:0.02:0.98, Corw(:, 16), 'k--', 'Linewidth', 1.25);
title('$corr(y_t, \tilde{y}_t)$', 'interpreter', 'latex', 'Fontsize', 12)
xlim([0, 1]); % xlabel('$\theta_w$', 'interpreter', 'latex', 'Fontsize', 12);

% 自相关系数
subplot(2, 3, 3)
plot(0.02:0.02:0.98, Autocorw(:, 3), 'k-', 0.02:0.02:0.98, Autocorw(:, 11), 'k--', 'Linewidth', 1.25)
title('$corr(\pi_t,\pi_{t-1})$', 'interpreter', 'latex', 'Fontsize', 12)
xlim([0, 1]); xlabel('$\theta_w$', 'interpreter', 'latex', 'Fontsize', 12);
subplot(2, 3, 6)

plot(0.02:0.02:0.98, Autocorw(:, 8), 'k-', 0.02:0.02:0.98, Autocorw(:, 16), 'k--', 'Linewidth', 1.25);
title('$corr(\tilde{y}_t, \tilde{y}_{t-1})$', 'interpreter', 'latex', 'Fontsize', 12)
xlim([0, 1]); xlabel('$\theta_w$', 'interpreter', 'latex', 'Fontsize', 12);
legend({'无指数化', '有指数化'}, 'Orientation', 'horizontal', ...
    'position', [0.43 0.78 0.15 0.10], 'FontSize', 10, 'FontName', 'kaiti')
legend('boxoff')




