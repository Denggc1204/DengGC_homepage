%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 一阶近似求解模型，四种情形：       %
% （1）原始方程形式                  %
% （2）非零稳态指数形式              %
% （3）零稳态指数形式                %
% （4）对数线性化形式                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 作者: Deng, Guichuan               %
% 邮箱: denggch@mail.sysu.edu.cn     %
% 最新更新时间: 2026/10/11           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
% 原始方程形式 %
%%%%%%%%%%%%%%%%
clear; 
close all;
clc; 

% 使用 dynare 函数进行求解和模拟
dynare primitive;

% 提取政策函数
ypos = strmatch('Y', M_.endo_names, 'exact');
cpos = strmatch('C', M_.endo_names, 'exact');
npos = strmatch('N', M_.endo_names, 'exact');
rpos = strmatch('R', M_.endo_names, 'exact');
wpos = strmatch('W', M_.endo_names, 'exact');
apos = strmatch('A', M_.endo_names, 'exact');

% 提取政策函数和转移函数系数
pfun = [oo_.dr.ys(ypos), oo_.dr.ghx(oo_.dr.inv_order_var(ypos)), oo_.dr.ghu(oo_.dr.inv_order_var(ypos));
        oo_.dr.ys(cpos), oo_.dr.ghx(oo_.dr.inv_order_var(cpos)), oo_.dr.ghu(oo_.dr.inv_order_var(cpos));
        oo_.dr.ys(npos), oo_.dr.ghx(oo_.dr.inv_order_var(npos)), oo_.dr.ghu(oo_.dr.inv_order_var(npos));
        oo_.dr.ys(rpos), oo_.dr.ghx(oo_.dr.inv_order_var(rpos)), oo_.dr.ghu(oo_.dr.inv_order_var(rpos));
        oo_.dr.ys(wpos), oo_.dr.ghx(oo_.dr.inv_order_var(wpos)), oo_.dr.ghu(oo_.dr.inv_order_var(wpos));
        oo_.dr.ys(apos), oo_.dr.ghx(oo_.dr.inv_order_var(apos)), oo_.dr.ghu(oo_.dr.inv_order_var(apos))];

% 系数保留4位有效数字
pfun = round(pfun, 4);

% 图3.2：脉冲反应函数
T = 30; t = 1:T; fsize = 12
h0 = figure 
set(h0, 'Units','centimeters','Position', [2 2 14 6]);
subplot(2, 3, 1)
plot(t, oo_.irfs.Y_etaa(t)*100, 'Color', [0,0,0], 'LineStyle', '-', 'LineWidth', 1.0)
xlim([0, T]);
ylabel('$10^{-2}$', 'interpreter', 'latex', 'Fontsize', 10);
title('$Y_t$', 'interpreter', 'latex', 'Fontsize', fsize);

subplot(2, 3, 2)
plot(t, oo_.irfs.C_etaa(t)*100, 'Color', [0,0,0], 'LineStyle', '-', 'LineWidth', 1.0)
xlim([0, T]);
ylabel('$10^{-2}$', 'interpreter', 'latex', 'Fontsize', 10);
title('$C_t$', 'interpreter', 'latex', 'Fontsize', fsize);

subplot(2, 3, 3)
plot(t, oo_.irfs.N_etaa(t)*100, 'Color', [0,0,0], 'LineStyle', '-', 'LineWidth', 1.0)
xlim([0, T]);
ylabel('$10^{-2}$', 'interpreter', 'latex', 'Fontsize', 10);
title('$N_t$', 'interpreter', 'latex', 'Fontsize', fsize);

subplot(2, 3, 4)
plot(t, oo_.irfs.W_etaa(t)*100, 'Color', [0,0,0], 'LineStyle', '-', 'LineWidth', 1.0)
xlim([0, T]);
ylabel('$10^{-2}$', 'interpreter', 'latex', 'Fontsize', 10);
title('$\hat{W}_t$', 'interpreter', 'latex', 'Fontsize', fsize);

subplot(2, 3, 5)
plot(t, oo_.irfs.R_etaa(t)*100, 'Color', [0,0,0], 'LineStyle', '-', 'LineWidth', 1.0)
xlim([0, T]);
ylabel('$10^{-2}$', 'interpreter', 'latex', 'Fontsize', 10);
title('$\hat{R}_t$', 'interpreter', 'latex', 'Fontsize', fsize);

subplot(2, 3, 6)
plot(t, oo_.irfs.A_etaa(t)*100, 'Color', [0,0,0], 'LineStyle', '-', 'LineWidth', 1.0)
xlim([0, T]);
ylabel('$10^{-2}$', 'interpreter', 'latex', 'Fontsize', 10);
title('$A_t$', 'interpreter', 'latex', 'Fontsize', fsize);


% 方差-协方差
co_var = oo_.var*10000;

% 表3.7：标准差
sd = round(sqrt(diag(co_var)),4);

% 表 3.8 第一行:与产出相关系数
coef = [1, co_var(ypos, cpos)/sqrt(co_var(ypos, ypos)*co_var(cpos, cpos)), ...
        co_var(ypos, npos)/sqrt(co_var(ypos, ypos)*co_var(npos, npos)), ...
        co_var(ypos, rpos)/sqrt(co_var(ypos, ypos)*co_var(rpos, rpos)), ...
        co_var(ypos, wpos)/sqrt(co_var(ypos, ypos)*co_var(wpos, wpos)), ...
        co_var(ypos, apos)/sqrt(co_var(ypos, ypos)*co_var(apos, apos))];

% 表3.8 第二至五行：自相关系数
auto_coef = [oo_.autocorr{1, 1}(ypos, ypos), oo_.autocorr{1, 1}(cpos, cpos), ...
             oo_.autocorr{1, 1}(npos, npos), oo_.autocorr{1, 1}(rpos, rpos), ...
             oo_.autocorr{1, 1}(wpos, wpos), oo_.autocorr{1, 1}(apos, apos);
             oo_.autocorr{1, 2}(ypos, ypos), oo_.autocorr{1, 1}(cpos, cpos), ...
             oo_.autocorr{1, 2}(npos, npos), oo_.autocorr{1, 1}(rpos, rpos), ...
             oo_.autocorr{1, 2}(wpos, wpos), oo_.autocorr{1, 1}(apos, apos);
             oo_.autocorr{1, 3}(ypos, ypos), oo_.autocorr{1, 1}(cpos, cpos), ...
             oo_.autocorr{1, 3}(npos, npos), oo_.autocorr{1, 1}(rpos, rpos), ...
             oo_.autocorr{1, 3}(wpos, wpos), oo_.autocorr{1, 1}(apos, apos);
             oo_.autocorr{1, 4}(ypos, ypos), oo_.autocorr{1, 1}(cpos, cpos), ...
             oo_.autocorr{1, 4}(npos, npos), oo_.autocorr{1, 1}(rpos, rpos), ...
             oo_.autocorr{1, 4}(wpos, wpos), oo_.autocorr{1, 1}(apos, apos);
             oo_.autocorr{1, 5}(ypos, ypos), oo_.autocorr{1, 1}(cpos, cpos), ...
             oo_.autocorr{1, 5}(npos, npos), oo_.autocorr{1, 1}(rpos, rpos), ...
             oo_.autocorr{1, 5}(wpos, wpos), oo_.autocorr{1, 1}(apos, apos)];

%%
%%%%%%%%%%%%%%%%%%%%
% 非零稳态指数形式 %
%%%%%%%%%%%%%%%%%%%%
clear; 
close all; 
clc;

% 求解和模拟
dynare exp_non_zero_ss;

% 提取核心变量位置
ypos = strmatch('lnY', M_.endo_names, 'exact');
cpos = strmatch('lnC', M_.endo_names, 'exact');
npos = strmatch('lnN', M_.endo_names, 'exact');
rpos = strmatch('lnR', M_.endo_names, 'exact');
wpos = strmatch('lnW', M_.endo_names, 'exact');
apos = strmatch('lnA', M_.endo_names, 'exact');

% 提取政策函数和转移函数系数
pfun = [oo_.dr.ys(ypos), oo_.dr.ghx(oo_.dr.inv_order_var(ypos)), oo_.dr.ghu(oo_.dr.inv_order_var(ypos));
        oo_.dr.ys(cpos), oo_.dr.ghx(oo_.dr.inv_order_var(cpos)), oo_.dr.ghu(oo_.dr.inv_order_var(cpos));
        oo_.dr.ys(npos), oo_.dr.ghx(oo_.dr.inv_order_var(npos)), oo_.dr.ghu(oo_.dr.inv_order_var(npos));
        oo_.dr.ys(rpos), oo_.dr.ghx(oo_.dr.inv_order_var(rpos)), oo_.dr.ghu(oo_.dr.inv_order_var(rpos));
        oo_.dr.ys(wpos), oo_.dr.ghx(oo_.dr.inv_order_var(wpos)), oo_.dr.ghu(oo_.dr.inv_order_var(wpos));
        oo_.dr.ys(apos), oo_.dr.ghx(oo_.dr.inv_order_var(apos)), oo_.dr.ghu(oo_.dr.inv_order_var(apos))];

% 系数保留4位有效数字
pfun = round(pfun, 4);

% 图3.3: 脉冲反应函数
T = 30; 
t = 1:T; 
fsize = 12
h0 = figure 
set(h0, 'Units','centimeters','Position', [2 2 14 6]);
subplot(2, 3, 1)
plot(t, oo_.irfs.lnY_etaa(t)*100, 'Color', [0,0,0], 'LineStyle', '-', 'LineWidth', 1.0)
xlim([0, T]);
ylabel('$10^{-2}$', 'interpreter', 'latex', 'Fontsize', 10);
title('$\ln Y_t$', 'interpreter', 'latex', 'Fontsize', fsize);

subplot(2, 3, 2)
plot(t, oo_.irfs.lnC_etaa(t)*100, 'Color', [0,0,0], 'LineStyle', '-', 'LineWidth', 1.0)
xlim([0, T]);
ylabel('$10^{-2}$', 'interpreter', 'latex', 'Fontsize', 10);
title('$\ln C_t$', 'interpreter', 'latex', 'Fontsize', fsize);

subplot(2, 3, 3)
plot(t, oo_.irfs.lnN_etaa(t)*100, 'Color', [0,0,0], 'LineStyle', '-', 'LineWidth', 1.0)
xlim([0, T]);
ylabel('$10^{-2}$', 'interpreter', 'latex', 'Fontsize', 10);
title('$\ln N_t$', 'interpreter', 'latex', 'Fontsize', fsize);

subplot(2, 3, 4)
plot(t, oo_.irfs.lnW_etaa(t)*100, 'Color', [0,0,0], 'LineStyle', '-', 'LineWidth', 1.0)
xlim([0, T]);
ylabel('$10^{-2}$', 'interpreter', 'latex', 'Fontsize', 10);
title('$\ln \hat{W}_t$', 'interpreter', 'latex', 'Fontsize', fsize);

subplot(2, 3, 5)
plot(t, oo_.irfs.lnR_etaa(t)*100, 'Color', [0,0,0], 'LineStyle', '-', 'LineWidth', 1.0)
xlim([0, T]);
ylabel('$10^{-2}$', 'interpreter', 'latex', 'Fontsize', 10);
title('$\ln \hat{R}_t$', 'interpreter', 'latex', 'Fontsize', fsize);

subplot(2, 3, 6)
plot(t, oo_.irfs.lnA_etaa(t)*100, 'Color', [0,0,0], 'LineStyle', '-', 'LineWidth', 1.0)
xlim([0, T]);
ylabel('$10^{-2}$', 'interpreter', 'latex', 'Fontsize', 10);
title('$\ln A_t$', 'interpreter', 'latex', 'Fontsize', fsize);


% 方差-协方差
co_var = oo_.var*10000;

% 表3.9: 标准差
sd = round(sqrt(diag(co_var)),4);

% 表3.10 第一行：与产出相关系数
coef = [1, co_var(ypos, cpos)/sqrt(co_var(ypos, ypos)*co_var(cpos, cpos)), ...
        co_var(ypos, npos)/sqrt(co_var(ypos, ypos)*co_var(npos, npos)), ...
        co_var(ypos, rpos)/sqrt(co_var(ypos, ypos)*co_var(rpos, rpos)), ...
        co_var(ypos, wpos)/sqrt(co_var(ypos, ypos)*co_var(wpos, wpos)), ...
        co_var(ypos, apos)/sqrt(co_var(ypos, ypos)*co_var(apos, apos))];

% 表3.10第二至五行：自相关系数
auto_coef = [oo_.autocorr{1, 1}(ypos, ypos), oo_.autocorr{1, 1}(cpos, cpos), ...
             oo_.autocorr{1, 1}(npos, npos), oo_.autocorr{1, 1}(rpos, rpos), ...
             oo_.autocorr{1, 1}(wpos, wpos), oo_.autocorr{1, 1}(apos, apos);
             oo_.autocorr{1, 2}(ypos, ypos), oo_.autocorr{1, 1}(cpos, cpos), ...
             oo_.autocorr{1, 2}(npos, npos), oo_.autocorr{1, 1}(rpos, rpos), ...
             oo_.autocorr{1, 2}(wpos, wpos), oo_.autocorr{1, 1}(apos, apos);
             oo_.autocorr{1, 3}(ypos, ypos), oo_.autocorr{1, 1}(cpos, cpos), ...
             oo_.autocorr{1, 3}(npos, npos), oo_.autocorr{1, 1}(rpos, rpos), ...
             oo_.autocorr{1, 3}(wpos, wpos), oo_.autocorr{1, 1}(apos, apos);
             oo_.autocorr{1, 4}(ypos, ypos), oo_.autocorr{1, 1}(cpos, cpos), ...
             oo_.autocorr{1, 4}(npos, npos), oo_.autocorr{1, 1}(rpos, rpos), ...
             oo_.autocorr{1, 4}(wpos, wpos), oo_.autocorr{1, 1}(apos, apos);
             oo_.autocorr{1, 5}(ypos, ypos), oo_.autocorr{1, 1}(cpos, cpos), ...
             oo_.autocorr{1, 5}(npos, npos), oo_.autocorr{1, 1}(rpos, rpos), ...
             oo_.autocorr{1, 5}(wpos, wpos), oo_.autocorr{1, 1}(apos, apos)];
 
%%
%%%%%%%%%%%%%%%%%%
% 零稳态指数形式 %
%%%%%%%%%%%%%%%%%%
clear; 
close all; 
clc;

% 求解和模拟
dynare exp_zero_ss;

% 提取核心变量位置
ypos = strmatch('y', M_.endo_names, 'exact');
cpos = strmatch('c', M_.endo_names, 'exact');
npos = strmatch('n', M_.endo_names, 'exact');
rpos = strmatch('r', M_.endo_names, 'exact');
wpos = strmatch('w', M_.endo_names, 'exact');
apos = strmatch('a', M_.endo_names, 'exact');

% 提取政策函数和转移函数系数
pfun = [oo_.dr.ys(ypos), oo_.dr.ghx(oo_.dr.inv_order_var(ypos)), oo_.dr.ghu(oo_.dr.inv_order_var(ypos));
        oo_.dr.ys(cpos), oo_.dr.ghx(oo_.dr.inv_order_var(cpos)), oo_.dr.ghu(oo_.dr.inv_order_var(cpos));
        oo_.dr.ys(npos), oo_.dr.ghx(oo_.dr.inv_order_var(npos)), oo_.dr.ghu(oo_.dr.inv_order_var(npos));
        oo_.dr.ys(rpos), oo_.dr.ghx(oo_.dr.inv_order_var(rpos)), oo_.dr.ghu(oo_.dr.inv_order_var(rpos));
        oo_.dr.ys(wpos), oo_.dr.ghx(oo_.dr.inv_order_var(wpos)), oo_.dr.ghu(oo_.dr.inv_order_var(wpos));
        oo_.dr.ys(apos), oo_.dr.ghx(oo_.dr.inv_order_var(apos)), oo_.dr.ghu(oo_.dr.inv_order_var(apos))];

% 系数保留4位有效数字
pfun = round(pfun, 4);

% 图3.4: 脉冲反应函数
T = 30; t = 1:T; fsize = 12
h0 = figure 
set(h0, 'Units','centimeters','Position', [2 2 14 6]);
subplot(2, 3, 1)
plot(t, oo_.irfs.y_etaa(t)*100, 'Color', [0,0,0], 'LineStyle', '-', 'LineWidth', 1.0)
xlim([0, T]);
ylabel('$10^{-2}$', 'interpreter', 'latex', 'Fontsize', 10);
title('$y_t$', 'interpreter', 'latex', 'Fontsize', fsize);

subplot(2, 3, 2)
plot(t, oo_.irfs.c_etaa(t)*100, 'Color', [0,0,0], 'LineStyle', '-', 'LineWidth', 1.0)
xlim([0, T]);
ylabel('$10^{-2}$', 'interpreter', 'latex', 'Fontsize', 10);
title('$c_t$', 'interpreter', 'latex', 'Fontsize', fsize);

subplot(2, 3, 3)
plot(t, oo_.irfs.n_etaa(t)*100, 'Color', [0,0,0], 'LineStyle', '-', 'LineWidth', 1.0)
xlim([0, T]);
ylabel('$10^{-2}$', 'interpreter', 'latex', 'Fontsize', 10);
title('$n_t$', 'interpreter', 'latex', 'Fontsize', fsize);

subplot(2, 3, 4)
plot(t, oo_.irfs.w_etaa(t)*100, 'Color', [0,0,0], 'LineStyle', '-', 'LineWidth', 1.0)
xlim([0, T]);
ylabel('$10^{-2}$', 'interpreter', 'latex', 'Fontsize', 10);
title('$\hat{w}_t$', 'interpreter', 'latex', 'Fontsize', fsize);

subplot(2, 3, 5)
plot(t, oo_.irfs.r_etaa(t)*100, 'Color', [0,0,0], 'LineStyle', '-', 'LineWidth', 1.0)
xlim([0, T]);
ylabel('$10^{-2}$', 'interpreter', 'latex', 'Fontsize', 10);
title('$\hat{r}_t$', 'interpreter', 'latex', 'Fontsize', fsize);

subplot(2, 3, 6)
plot(t, oo_.irfs.a_etaa(t)*100, 'Color', [0,0,0], 'LineStyle', '-', 'LineWidth', 1.0)
xlim([0, T]);
ylabel('$10^{-2}$', 'interpreter', 'latex', 'Fontsize', 10);
title('$a_t$', 'interpreter', 'latex', 'Fontsize', fsize);


% 方差-协方差
co_var = oo_.var*10000;

% 表3.11: 标准差
sd = round(sqrt(diag(co_var)),4);

% 表3.12 第一行：与产出相关系数
coef = [1, co_var(ypos, cpos)/sqrt(co_var(ypos, ypos)*co_var(cpos, cpos)), ...
        co_var(ypos, npos)/sqrt(co_var(ypos, ypos)*co_var(npos, npos)), ...
        co_var(ypos, rpos)/sqrt(co_var(ypos, ypos)*co_var(rpos, rpos)), ...
        co_var(ypos, wpos)/sqrt(co_var(ypos, ypos)*co_var(wpos, wpos)), ...
        co_var(ypos, apos)/sqrt(co_var(ypos, ypos)*co_var(apos, apos))];
    
% 表3.12 第二至五行：自相关系数
auto_coef = [oo_.autocorr{1, 1}(ypos, ypos), oo_.autocorr{1, 1}(cpos, cpos), ...
             oo_.autocorr{1, 1}(npos, npos), oo_.autocorr{1, 1}(rpos, rpos), ...
             oo_.autocorr{1, 1}(wpos, wpos), oo_.autocorr{1, 1}(apos, apos);
             oo_.autocorr{1, 2}(ypos, ypos), oo_.autocorr{1, 1}(cpos, cpos), ...
             oo_.autocorr{1, 2}(npos, npos), oo_.autocorr{1, 1}(rpos, rpos), ...
             oo_.autocorr{1, 2}(wpos, wpos), oo_.autocorr{1, 1}(apos, apos);
             oo_.autocorr{1, 3}(ypos, ypos), oo_.autocorr{1, 1}(cpos, cpos), ...
             oo_.autocorr{1, 3}(npos, npos), oo_.autocorr{1, 1}(rpos, rpos), ...
             oo_.autocorr{1, 3}(wpos, wpos), oo_.autocorr{1, 1}(apos, apos);
             oo_.autocorr{1, 4}(ypos, ypos), oo_.autocorr{1, 1}(cpos, cpos), ...
             oo_.autocorr{1, 4}(npos, npos), oo_.autocorr{1, 1}(rpos, rpos), ...
             oo_.autocorr{1, 4}(wpos, wpos), oo_.autocorr{1, 1}(apos, apos);
             oo_.autocorr{1, 5}(ypos, ypos), oo_.autocorr{1, 1}(cpos, cpos), ...
             oo_.autocorr{1, 5}(npos, npos), oo_.autocorr{1, 1}(rpos, rpos), ...
             oo_.autocorr{1, 5}(wpos, wpos), oo_.autocorr{1, 1}(apos, apos)];
  

%%         
%%%%%%%%%%%%%%%%%%
% 对数线性化形式 %
%%%%%%%%%%%%%%%%%%
clear; 
close all; 
clc;

% 求解和模拟
dynare loglinear;

% 提取政策方程
ypos = strmatch('y', M_.endo_names, 'exact');
cpos = strmatch('c', M_.endo_names, 'exact');
npos = strmatch('n', M_.endo_names, 'exact');
rpos = strmatch('r', M_.endo_names, 'exact');
wpos = strmatch('w', M_.endo_names, 'exact');
apos = strmatch('a', M_.endo_names, 'exact');

% 提取政策函数和转移函数系数
pfun = [oo_.dr.ys(ypos), oo_.dr.ghx(oo_.dr.inv_order_var(ypos)), oo_.dr.ghu(oo_.dr.inv_order_var(ypos));
        oo_.dr.ys(cpos), oo_.dr.ghx(oo_.dr.inv_order_var(cpos)), oo_.dr.ghu(oo_.dr.inv_order_var(cpos));
        oo_.dr.ys(npos), oo_.dr.ghx(oo_.dr.inv_order_var(npos)), oo_.dr.ghu(oo_.dr.inv_order_var(npos));
        oo_.dr.ys(rpos), oo_.dr.ghx(oo_.dr.inv_order_var(rpos)), oo_.dr.ghu(oo_.dr.inv_order_var(rpos));
        oo_.dr.ys(wpos), oo_.dr.ghx(oo_.dr.inv_order_var(wpos)), oo_.dr.ghu(oo_.dr.inv_order_var(wpos));
        oo_.dr.ys(apos), oo_.dr.ghx(oo_.dr.inv_order_var(apos)), oo_.dr.ghu(oo_.dr.inv_order_var(apos))];

% 系数保留4位有效数字
pfun = round(pfun, 4);

% 图3.5: 脉冲反应函数
T = 30; t = 1:T; fsize = 12
h0 = figure 
set(h0, 'Units','centimeters','Position', [2 2 14 6]);
subplot(2, 3, 1)
plot(t, oo_.irfs.y_etaa(t)*100, 'Color', [0,0,0], 'LineStyle', '-', 'LineWidth', 1.0)
xlim([0, T]);
ylabel('$10^{-2}$', 'interpreter', 'latex', 'Fontsize', 10);
title('$y_t$', 'interpreter', 'latex', 'Fontsize', fsize);

subplot(2, 3, 2)
plot(t, oo_.irfs.c_etaa(t)*100, 'Color', [0,0,0], 'LineStyle', '-', 'LineWidth', 1.0)
xlim([0, T]);
ylabel('$10^{-2}$', 'interpreter', 'latex', 'Fontsize', 10);
title('$c_t$', 'interpreter', 'latex', 'Fontsize', fsize);

subplot(2, 3, 3)
plot(t, oo_.irfs.n_etaa(t)*100, 'Color', [0,0,0], 'LineStyle', '-', 'LineWidth', 1.0)
xlim([0, T]);
ylabel('$10^{-2}$', 'interpreter', 'latex', 'Fontsize', 10);
title('$n_t$', 'interpreter', 'latex', 'Fontsize', fsize);

subplot(2, 3, 4)
plot(t, oo_.irfs.w_etaa(t)*100, 'Color', [0,0,0], 'LineStyle', '-', 'LineWidth', 1.0)
xlim([0, T]);
ylabel('$10^{-2}$', 'interpreter', 'latex', 'Fontsize', 10);
title('$\hat{w}_t$', 'interpreter', 'latex', 'Fontsize', fsize);

subplot(2, 3, 5)
plot(t, oo_.irfs.r_etaa(t)*100, 'Color', [0,0,0], 'LineStyle', '-', 'LineWidth', 1.0)
xlim([0, T]);
ylabel('$10^{-2}$', 'interpreter', 'latex', 'Fontsize', 10);
title('$\hat{r}_t$', 'interpreter', 'latex', 'Fontsize', fsize);

subplot(2, 3, 6)
plot(t, oo_.irfs.a_etaa(t)*100, 'Color', [0,0,0], 'LineStyle', '-', 'LineWidth', 1.0)
xlim([0, T]);
ylabel('$10^{-2}$', 'interpreter', 'latex', 'Fontsize', 10);
title('$a_t$', 'interpreter', 'latex', 'Fontsize', fsize);


% 方差-协方差
co_var = oo_.var*10000;

% 表3.13: 标准差
sd = round(sqrt(diag(co_var)),4);

% 表3.14第一行：与产出相关系数
coef = [1, co_var(ypos, cpos)/sqrt(co_var(ypos, ypos)*co_var(cpos, cpos)), ...
        co_var(ypos, npos)/sqrt(co_var(ypos, ypos)*co_var(npos, npos)), ...
        co_var(ypos, rpos)/sqrt(co_var(ypos, ypos)*co_var(rpos, rpos)), ...
        co_var(ypos, wpos)/sqrt(co_var(ypos, ypos)*co_var(wpos, wpos)), ...
        co_var(ypos, apos)/sqrt(co_var(ypos, ypos)*co_var(apos, apos))];

% 表3.14第二至五行：自相关系数
auto_coef = [oo_.autocorr{1, 1}(ypos, ypos), oo_.autocorr{1, 1}(cpos, cpos), ...
             oo_.autocorr{1, 1}(npos, npos), oo_.autocorr{1, 1}(rpos, rpos), ...
             oo_.autocorr{1, 1}(wpos, wpos), oo_.autocorr{1, 1}(apos, apos);
             oo_.autocorr{1, 2}(ypos, ypos), oo_.autocorr{1, 1}(cpos, cpos), ...
             oo_.autocorr{1, 2}(npos, npos), oo_.autocorr{1, 1}(rpos, rpos), ...
             oo_.autocorr{1, 2}(wpos, wpos), oo_.autocorr{1, 1}(apos, apos);
             oo_.autocorr{1, 3}(ypos, ypos), oo_.autocorr{1, 1}(cpos, cpos), ...
             oo_.autocorr{1, 3}(npos, npos), oo_.autocorr{1, 1}(rpos, rpos), ...
             oo_.autocorr{1, 3}(wpos, wpos), oo_.autocorr{1, 1}(apos, apos);
             oo_.autocorr{1, 4}(ypos, ypos), oo_.autocorr{1, 1}(cpos, cpos), ...
             oo_.autocorr{1, 4}(npos, npos), oo_.autocorr{1, 1}(rpos, rpos), ...
             oo_.autocorr{1, 4}(wpos, wpos), oo_.autocorr{1, 1}(apos, apos);
             oo_.autocorr{1, 5}(ypos, ypos), oo_.autocorr{1, 1}(cpos, cpos), ...
             oo_.autocorr{1, 5}(npos, npos), oo_.autocorr{1, 1}(rpos, rpos), ...
             oo_.autocorr{1, 5}(wpos, wpos), oo_.autocorr{1, 1}(apos, apos)];

