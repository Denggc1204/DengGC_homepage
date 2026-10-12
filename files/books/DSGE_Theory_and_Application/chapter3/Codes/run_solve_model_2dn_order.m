%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 二阶近似求解模型，四种情形：       %
% （1）原始方程形式                  %
% （2）非零稳态指数形式              %
% （3）零稳态指数形式                %
% （4）对数线性化形式                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 作者: Deng, Guichuan               %
% 邮箱: denggch@mail.sysu.edu.cn     %
% 最新更新时间: 2026/10/11           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; 
close all; 
clc;
mkdir('./results_sec');

%%%%%%%%%%%%%%%%%%%%
% 原始方程形式
%%%%%%%%%%%%%%%%%%%%
clear; 
close all; 
clc;

% 求解和模拟
dynare primitive_sec;

% 保存结果
po = oo_; pm = M_;
save('./results_sec./po.mat', 'po');
save('./results_sec./pm.mat', 'pm');

%%%%%%%%%%%%%%%%%%%%
% 非零稳态指数形式
%%%%%%%%%%%%%%%%%%%%
clear; 
close all; 
clc;

% 求解和模拟
dynare exp_non_zero_ss_sec;

% 保存结果
eno = oo_; enm = M_;
save('./results_sec./eno.mat', 'eno');
save('./results_sec./enm.mat', 'enm');

%%%%%%%%%%%%%%%%%%%%
% 零稳态指数形式
%%%%%%%%%%%%%%%%%%%%
clear; 
close all;
clc; 

% 求解和模拟
dynare exp_zero_ss_sec;

% 保存结果
eo = oo_; em = M_;
save('./results_sec./eo.mat', 'eo');
save('./results_sec./em.mat', 'em');

%%%%%%%%%%%%%%%%%%%%
% 对数线性化形式
%%%%%%%%%%%%%%%%%%%%
clear; 
close all;
clc; 

% 求解和模拟
dynare loglinear_sec;

% 保存结果
lo = oo_; lm = M_;
save('./results_sec./lo.mat', 'lo');
save('./results_sec./lm.mat', 'lm');

%%
%%%%%%%%%%%%%%%%
% 模型分析     %
%%%%%%%%%%%%%%%%
clear; 
close all; 
clc;

% 加载结果
load('./results_sec./po.mat');
load('./results_sec./pm.mat');
load('./results_sec./eno.mat');
load('./results_sec./enm.mat');
load('./results_sec./eo.mat');
load('./results_sec./em.mat');
load('./results_sec./lo.mat');
load('./results_sec./lm.mat');

% 提取变量位置信息
uposp = strmatch('U',    pm.endo_names,  'exact');
yposp = strmatch('Y',    pm.endo_names,  'exact');
cposp = strmatch('C',    pm.endo_names,  'exact');
nposp = strmatch('N',    pm.endo_names,  'exact');
rposp = strmatch('R',    pm.endo_names,  'exact');
wposp = strmatch('W',    pm.endo_names,  'exact');
uposen = strmatch('U',   enm.endo_names, 'exact');
yposen = strmatch('lnY', enm.endo_names, 'exact');
cposen = strmatch('lnC', enm.endo_names, 'exact');
nposen = strmatch('lnN', enm.endo_names, 'exact');
rposen = strmatch('lnR', enm.endo_names, 'exact');
wposen = strmatch('lnW', enm.endo_names, 'exact');
upose = strmatch('U',    em.endo_names,  'exact');
ypose = strmatch('y',    em.endo_names,  'exact');
cpose = strmatch('c',    em.endo_names,  'exact');
npose = strmatch('n',    em.endo_names,  'exact');
rpose = strmatch('r',    em.endo_names,  'exact');
wpose = strmatch('w',    em.endo_names,  'exact');
uposl = strmatch('U',    lm.endo_names,  'exact');
yposl = strmatch('y',    lm.endo_names,  'exact');
cposl = strmatch('c',    lm.endo_names,  'exact');
nposl = strmatch('n',    lm.endo_names,  'exact');
rposl = strmatch('r',    lm.endo_names,  'exact');
wposl = strmatch('w',    lm.endo_names,  'exact');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 图3.6: 脉冲反应函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = 30; t = 1:T; fsize = 12
h0 = figure 
set(h0, 'Units','centimeters','Position', [2 2 14 8]);
subplot(2, 3, 1)
plot(t, po.irfs.Y_etaa(t)/po.dr.ys(yposp)*100, 'Color', [0,0,0], 'LineStyle', ':', 'LineWidth', 1.0);
hold on;
plot(t, eno.irfs.lnY_etaa(t)*100, 'Color', [0.2,0.2,0.2], 'LineStyle', '-.', 'LineWidth', 1.0);
plot(t, eo.irfs.y_etaa(t)*100, 'Color', [0.4,0.4,0.4], 'LineStyle', '--', 'LineWidth', 1.0);
plot(t, lo.irfs.y_etaa(t)*100, 'Color', [0.6,0.6,0.6], 'LineStyle', '-', 'LineWidth', 1.0);
hold off;
xlim([1, T]);
ylabel('$10^{-2}$', 'interpreter', 'latex', 'Fontsize', 10);
title('$y_t$', 'interpreter', 'latex', 'Fontsize', fsize);

subplot(2, 3, 2)
plot(t, po.irfs.C_etaa(t)/po.dr.ys(cposp)*100, 'Color', [0,0,0], 'LineStyle', ':', 'LineWidth', 1.0);
hold on;
plot(t, eno.irfs.lnC_etaa(t)*100, 'Color', [0.2,0.2,0.2], 'LineStyle', '-.', 'LineWidth', 1.0);
plot(t, eo.irfs.c_etaa(t)*100, 'Color', [0.4,0.4,0.4], 'LineStyle', '--', 'LineWidth', 1.0);
plot(t, lo.irfs.c_etaa(t)*100, 'Color', [0.6,0.6,0.6], 'LineStyle', '-', 'LineWidth', 1.0);
hold off;
xlim([1, T]);
ylabel('$10^{-2}$', 'interpreter', 'latex', 'Fontsize', 10);
title('$c_t$', 'interpreter', 'latex', 'Fontsize', fsize);

subplot(2, 3, 3)
plot(t, po.irfs.N_etaa(t)/po.dr.ys(nposp)*100, 'Color', [0,0,0], 'LineStyle', ':', 'LineWidth', 1.0);
hold on;
plot(t, eno.irfs.lnN_etaa(t)*100, 'Color', [0.2,0.2,0.2], 'LineStyle', '-.', 'LineWidth', 1.0);
plot(t, eo.irfs.n_etaa(t)*100, 'Color', [0.4,0.4,0.4], 'LineStyle', '--', 'LineWidth', 1.0);
plot(t, lo.irfs.n_etaa(t)*100, 'Color', [0.6,0.6,0.6], 'LineStyle', '-', 'LineWidth', 1.0);
hold off;
xlim([1, T]);
ylabel('$10^{-2}$', 'interpreter', 'latex', 'Fontsize', 10);
title('$n_t$', 'interpreter', 'latex', 'Fontsize', fsize);

subplot(2, 3, 4)
plot(t, po.irfs.R_etaa(t)/po.dr.ys(rposp)*100, 'Color', [0,0,0], 'LineStyle', ':', 'LineWidth', 1.0);
hold on;
plot(t, eno.irfs.lnR_etaa(t)*100, 'Color', [0.2,0.2,0.2], 'LineStyle', '-.', 'LineWidth', 1.0);
plot(t, eo.irfs.r_etaa(t)*100, 'Color', [0.4,0.4,0.4], 'LineStyle', '--', 'LineWidth', 1.0);
plot(t, lo.irfs.r_etaa(t)*100, 'Color', [0.6,0.6,0.6], 'LineStyle', '-', 'LineWidth', 1.0);
hold off;
xlim([1, T]);
ylabel('$10^{-2}$', 'interpreter', 'latex', 'Fontsize', 10);
title('$\hat{r}_t$', 'interpreter', 'latex', 'Fontsize', fsize);

subplot(2, 3, 5)
plot(t, po.irfs.W_etaa(t)/po.dr.ys(wposp)*100, 'Color', [0,0,0], 'LineStyle', ':', 'LineWidth', 1.0);
hold on;
plot(t, eno.irfs.lnW_etaa(t)*100, 'Color', [0.2,0.2,0.2], 'LineStyle', '-.', 'LineWidth', 1.0);
plot(t, eo.irfs.w_etaa(t)*100, 'Color', [0.4,0.4,0.4], 'LineStyle', '--', 'LineWidth', 1.0);
plot(t, lo.irfs.w_etaa(t)*100, 'Color', [0.6,0.6,0.6], 'LineStyle', '-', 'LineWidth', 1.0);
hold off;
xlim([1, T]);
ylabel('$10^{-2}$', 'interpreter', 'latex', 'Fontsize', 10);
title('$\hat{w}_t$', 'interpreter', 'latex', 'Fontsize', fsize);

subplot(2, 3, 6)
plot(t, po.irfs.U_etaa(t)*100, 'Color', [0,0,0], 'LineStyle', ':', 'LineWidth', 1.0);
hold on;
plot(t, eno.irfs.U_etaa(t)*100, 'Color', [0.2,0.2,0.2], 'LineStyle', '-.', 'LineWidth', 1.0);
plot(t, eo.irfs.U_etaa(t)*100, 'Color', [0.4,0.4,0.4], 'LineStyle', '--', 'LineWidth', 1.0);
plot(t, lo.irfs.U_etaa(t)*100, 'Color', [0.6,0.6,0.6], 'LineStyle', '-', 'LineWidth', 1.0);
hold off;
xlim([1, T]);
ylabel('$10^{-2}$', 'interpreter', 'latex', 'Fontsize', 10);
title('$U_t$', 'interpreter', 'latex', 'Fontsize', fsize);
legend({'原始方程', '非零稳态', '零稳态', '线性方程'}, 'FontSize', 10, ...
    'Position', [0.45 -0.10 0.15 0.25], 'Orientation', 'horizontal');
legend('boxoff')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 表3.15: 无条件均值和条件均值
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mean = round([
    po.mean(yposp), po.mean(cposp), po.mean(nposp), po.mean(rposp), po.mean(wposp), po.mean(uposp);
    po.dr.ys(yposp)+0.5*po.dr.ghs2(po.dr.inv_order_var(yposp)), ...
    po.dr.ys(cposp)+0.5*po.dr.ghs2(po.dr.inv_order_var(cposp)), ...
    po.dr.ys(nposp)+0.5*po.dr.ghs2(po.dr.inv_order_var(nposp)), ...
    po.dr.ys(rposp)+0.5*po.dr.ghs2(po.dr.inv_order_var(rposp)), ...
    po.dr.ys(wposp)+0.5*po.dr.ghs2(po.dr.inv_order_var(wposp)), ...
    po.dr.ys(uposp)+0.5*po.dr.ghs2(po.dr.inv_order_var(uposp));
    po.dr.ys(yposp), po.dr.ys(cposp), po.dr.ys(nposp), po.dr.ys(rposp), po.dr.ys(wposp), po.dr.ys(uposp);
    %
    eno.mean(yposen), eno.mean(cposen), eno.mean(nposen),... 
    eno.mean(rposen), eno.mean(wposen), eno.mean(uposen);
    eno.dr.ys(yposen)+0.5*eno.dr.ghs2(eno.dr.inv_order_var(yposen)), ...
    eno.dr.ys(cposen)+0.5*eno.dr.ghs2(eno.dr.inv_order_var(cposen)), ...
    eno.dr.ys(nposen)+0.5*eno.dr.ghs2(eno.dr.inv_order_var(nposen)), ...
    eno.dr.ys(rposen)+0.5*eno.dr.ghs2(eno.dr.inv_order_var(rposen)), ...
    eno.dr.ys(wposen)+0.5*eno.dr.ghs2(eno.dr.inv_order_var(wposen)), ...
    eno.dr.ys(uposen)+0.5*eno.dr.ghs2(eno.dr.inv_order_var(uposen));
    eno.dr.ys(yposen), eno.dr.ys(cposen), eno.dr.ys(nposen),... 
    eno.dr.ys(rposen), eno.dr.ys(wposen), eno.dr.ys(uposen);
    %
    eo.mean(ypose), eo.mean(cpose), eo.mean(npose), eo.mean(rpose), eo.mean(wpose), eo.mean(upose);
    eo.dr.ys(ypose)+0.5*eo.dr.ghs2(eo.dr.inv_order_var(ypose)), ...
    eo.dr.ys(cpose)+0.5*eo.dr.ghs2(eo.dr.inv_order_var(cpose)), ...
    eo.dr.ys(npose)+0.5*eo.dr.ghs2(eo.dr.inv_order_var(npose)), ...
    eo.dr.ys(rpose)+0.5*eo.dr.ghs2(eo.dr.inv_order_var(rpose)), ...
    eo.dr.ys(wpose)+0.5*eo.dr.ghs2(eo.dr.inv_order_var(wpose)), ...
    eo.dr.ys(upose)+0.5*eo.dr.ghs2(eo.dr.inv_order_var(upose));
    eo.dr.ys(ypose), eo.dr.ys(cpose), eo.dr.ys(npose), eo.dr.ys(rpose), eo.dr.ys(wpose), eo.dr.ys(upose);
    %
    lo.mean(yposl), lo.mean(cposl), lo.mean(nposl), lo.mean(rposl), lo.mean(wposl), lo.mean(uposl);
    lo.dr.ys(yposl), lo.dr.ys(cposl), lo.dr.ys(nposl), lo.dr.ys(rposl), lo.dr.ys(wposl), lo.dr.ys(uposl);
    lo.dr.ys(yposl), lo.dr.ys(cposl), lo.dr.ys(nposl), lo.dr.ys(rposl), lo.dr.ys(wposl), lo.dr.ys(uposl)], 4);

%%%%%%%%%%%%%
% 解析解
%%%%%%%%%%%%%
bet = 0.99; sig = 2.10; phi = 3.00; eps = 10; mu = (eps-1)/eps; rho = 0.88; siga = 0.01;
Uss = (mu^((1-sig)/(sig+phi))/(1-sig) - mu^((1+phi)/(sig+phi))/(1+phi))/(1-bet);
U_un = Uss*(1 + (1-sig)^2*(1+phi)^2/(2*(1-rho^2)*(sig+phi)^2)*siga^2);
U_con = Uss*(1 + (1-bet)*(1-sig)^2*(1+phi)^2/(2*(1-rho^2)*(sig+phi)^2)*(1/(1-bet)-rho^2/(1-bet*rho^2))*siga^2);
sdA = siga/sqrt(1-rho^2);
sdY = (1+phi)/(sig+phi)*sdA;
sdN = (sig-1)/(sig+phi)*sdA;
sdR = sig*(1+phi)*(1-rho)/(sig+phi)*sdA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 表3.16：标准差
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vol = round(sqrt([
    po.var(yposp, yposp)/po.dr.ys(yposp)^2,    po.var(cposp, cposp)/po.dr.ys(cposp)^2,...    
    po.var(nposp, nposp)/po.dr.ys(nposp)^2,    po.var(rposp, rposp)/po.dr.ys(rposp)^2,...    
    po.var(wposp, wposp)/po.dr.ys(wposp)^2,    po.var(uposp, uposp);
    eno.var(yposen, yposen), eno.var(cposen, cposen), eno.var(nposen, nposen), ...
    eno.var(rposen, rposen), eno.var(wposen, wposen), eno.var(uposen, uposen);
    eo.var(ypose, ypose),    eo.var(cpose, cpose),    eo.var(npose, npose),...    
    eo.var(rpose, rpose),    eo.var(wpose, wpose),    eo.var(upose, upose);
    lo.var(yposl, yposl),    lo.var(cposl, cposl),    lo.var(nposl, nposl),...    
    lo.var(rposl, rposl),    lo.var(wposl, wposl),    lo.var(uposl, uposl)])*100, 4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 表3.17：与产出相关系数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cor = round([
    1, po.var(yposp, cposp)/sqrt(po.var(yposp, yposp)*po.var(cposp, cposp)),... 
       po.var(yposp, nposp)/sqrt(po.var(yposp, yposp)*po.var(nposp, nposp)),... 
       po.var(yposp, rposp)/sqrt(po.var(yposp, yposp)*po.var(rposp, rposp)),...
       po.var(yposp, wposp)/sqrt(po.var(yposp, yposp)*po.var(wposp, wposp));
    1, eno.var(yposen, cposen)/sqrt(eno.var(yposen, yposen)*eno.var(cposen, cposen)),... 
       eno.var(yposen, nposen)/sqrt(eno.var(yposen, yposen)*eno.var(nposen, nposen)),... 
       eno.var(yposen, rposen)/sqrt(eno.var(yposen, yposen)*eno.var(rposen, rposen)),...
       eno.var(yposen, wposen)/sqrt(eno.var(yposen, yposen)*eno.var(wposen, wposen));
    1, eo.var(ypose, cpose)/sqrt(eo.var(ypose, ypose)*eo.var(cpose, cpose)),... 
       eo.var(ypose, npose)/sqrt(eo.var(ypose, ypose)*eo.var(npose, npose)),... 
       eo.var(ypose, rpose)/sqrt(eo.var(ypose, ypose)*eo.var(rpose, rpose)),...
       eo.var(ypose, wpose)/sqrt(eo.var(ypose, ypose)*eo.var(wpose, wpose));
    1, lo.var(yposl, cposl)/sqrt(lo.var(yposl, yposl)*lo.var(cposl, cposl)),... 
       lo.var(yposl, nposl)/sqrt(lo.var(yposl, yposl)*lo.var(nposl, nposl)),... 
       lo.var(yposl, rposl)/sqrt(lo.var(yposl, yposl)*lo.var(rposl, rposl)),...
       lo.var(yposl, wposl)/sqrt(lo.var(yposl, yposl)*lo.var(wposl, wposl))], 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 表3.18： 自相关系数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Autocor = round([
    po.autocorr{1, 1}(yposp, yposp), po.autocorr{1, 1}(cposp, cposp),... 
    po.autocorr{1, 1}(nposp, nposp), po.autocorr{1, 1}(rposp, rposp), po.autocorr{1, 1}(wposp, wposp);
    eno.autocorr{1, 1}(yposen, yposen), eno.autocorr{1, 1}(cposen, cposen),... 
    eno.autocorr{1, 1}(nposen, nposen), eno.autocorr{1, 1}(rposen, rposen), eno.autocorr{1, 1}(wposen, wposen);
    eo.autocorr{1, 1}(ypose, ypose), eo.autocorr{1, 1}(cpose, cpose),... 
    eo.autocorr{1, 1}(npose, npose), eo.autocorr{1, 1}(rpose, rpose), eo.autocorr{1, 1}(wpose, wpose);
    lo.autocorr{1, 1}(yposl, yposl), lo.autocorr{1, 1}(cposl, cposl),... 
    lo.autocorr{1, 1}(nposl, nposl), lo.autocorr{1, 1}(rposl, rposl), lo.autocorr{1, 1}(wposl, wposl)], 3);


