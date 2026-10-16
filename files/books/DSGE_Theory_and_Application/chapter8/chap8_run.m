%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 第8章：求解非对称两国模型和小国模型                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 作者: Deng, Guichuan             
% 邮箱: denggch@mail.sysu.edu.cn   
% 最新更新时间: 2026/10/16         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 新建子文件夹
mkdir('./result')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 求解非对称两国模型：
% 国际金融市场完备程度（omega）：从0 至 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; 
close all; 
clc;

% 循环求解
for ii = 1:21
    % 设定参数值
    ome = (ii-1)/20; 
    % 求解和模拟
    dynare asy_two_country.mod
    % 保存求解结果
    asyo = ['asy', num2str(ii)];
    assignin('base', asyo, oo_);
    filepath = ['./result./asyo', num2str(ii), '.mat']
    save(filepath, asyo);
end

% 保存 M_ 文件
asym = M_;
save('./result./asym.mat', 'asym');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 非对称两国模型动态分析
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; 
close all; 
clc;

% 加载数据
load('./result./asyo1')  % omega = 0.00;
load('./result./asyo21') % omega = 1.00;

% 设置脉冲反应期数
T = 10; t = 1:T;

% 图8.1：本国生产率冲击
h0 = figure
set(h0,'Units','centimeters','Position',[2 2 14 12]);
subplot(3, 3, 1);
plot(t, asy1.irfs.y_etaa(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, asy21.irfs.y_etaa(t)*100, 'Color', [0.4,0.4,0.4], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$y_t$', 'interpreter', 'latex', 'fontsize', 12);

subplot(3, 3, 2);
plot(t, asy1.irfs.yst_etaa(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, asy21.irfs.yst_etaa(t)*100, 'Color', [0.4,0.4,0.4], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$y_t^\ast$', 'interpreter', 'latex', 'fontsize', 12);

subplot(3, 3, 3)
plot(t, asy1.irfs.c_etaa(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, asy21.irfs.c_etaa(t)*100, 'Color', [0.4,0.4,0.4], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$c_t$', 'interpreter', 'latex', 'fontsize', 12);

subplot(3, 3, 4)
plot(t, asy1.irfs.cst_etaa(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, asy21.irfs.cst_etaa(t)*100, 'Color', [0.4,0.4,0.4], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$c_t^\ast$', 'interpreter', 'latex', 'fontsize', 12);

subplot(3, 3, 5)
plot(t, asy1.irfs.piH_etaa(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, asy21.irfs.piH_etaa(t)*100, 'Color', [0.4,0.4,0.4], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$\pi_{Ht}$', 'interpreter', 'latex', 'fontsize', 12);

subplot(3, 3, 6)
plot(t, asy1.irfs.piFst_etaa(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, asy21.irfs.piFst_etaa(t)*100, 'Color', [0.4,0.4,0.4], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$\pi_{Ft}^\ast$', 'interpreter', 'latex', 'fontsize', 12);

subplot(3, 3, 7)
plot(t, asy1.irfs.pi_etaa(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, asy21.irfs.pi_etaa(t)*100, 'Color', [0.4,0.4,0.4], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$\pi_t$', 'interpreter', 'latex', 'fontsize', 12);

subplot(3, 3, 8)
plot(t, asy1.irfs.pist_etaa(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, asy21.irfs.pist_etaa(t)*100, 'Color', [0.4,0.4,0.4], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$\pi_t^\ast$', 'interpreter', 'latex', 'fontsize', 12);

subplot(3, 3, 9)
plot(t, asy1.irfs.q_etaa(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, asy21.irfs.q_etaa(t)*100, 'Color', [0.4,0.4,0.4], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$q_t$', 'interpreter', 'latex', 'fontsize', 12);

legend({'$\omega=0.0$', '$\omega=1.0$'}, 'interpreter', 'latex', 'position', [0.47 0.10 0.20 0.15], ...
    'fontname', 'kaiti', 'fontsize', 10.5);
legend('boxoff');

% 图8.2：外国生产率冲击
h0 = figure
set(h0,'Units','centimeters','Position',[2 2 14 12]);
subplot(3, 3, 1);
plot(t, asy1.irfs.y_etaast(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, asy21.irfs.y_etaast(t)*100, 'Color', [0.4,0.4,0.4], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$y_t$', 'interpreter', 'latex', 'fontsize', 12);

subplot(3, 3, 2);
plot(t, asy1.irfs.yst_etaast(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, asy21.irfs.yst_etaast(t)*100, 'Color', [0.4,0.4,0.4], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$y_t^\ast$', 'interpreter', 'latex', 'fontsize', 12);

subplot(3, 3, 3)
plot(t, asy1.irfs.c_etaast(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, asy21.irfs.c_etaast(t)*100, 'Color', [0.4,0.4,0.4], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$c_t$', 'interpreter', 'latex', 'fontsize', 12);

subplot(3, 3, 4)
plot(t, asy1.irfs.cst_etaast(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, asy21.irfs.cst_etaast(t)*100, 'Color', [0.4,0.4,0.4], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$c_t^\ast$', 'interpreter', 'latex', 'fontsize', 12);

subplot(3, 3, 5)
plot(t, asy1.irfs.piH_etaast(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, asy21.irfs.piH_etaast(t)*100, 'Color', [0.4,0.4,0.4], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$\pi_{Ht}$', 'interpreter', 'latex', 'fontsize', 12);

subplot(3, 3, 6)
plot(t, asy1.irfs.piFst_etaast(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, asy21.irfs.piFst_etaast(t)*100, 'Color', [0.4,0.4,0.4], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$\pi_{Ft}^\ast$', 'interpreter', 'latex', 'fontsize', 12);

subplot(3, 3, 7)
plot(t, asy1.irfs.pi_etaast(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, asy21.irfs.pi_etaast(t)*100, 'Color', [0.4,0.4,0.4], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$\pi_t$', 'interpreter', 'latex', 'fontsize', 12);

subplot(3, 3, 8)
plot(t, asy1.irfs.pist_etaast(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, asy21.irfs.pist_etaast(t)*100, 'Color', [0.4,0.4,0.4], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$\pi_t^\ast$', 'interpreter', 'latex', 'fontsize', 12);

subplot(3, 3, 9)
plot(t, asy1.irfs.q_etaast(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, asy21.irfs.q_etaast(t)*100, 'Color', [0.4,0.4,0.4], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$q_t$', 'interpreter', 'latex', 'fontsize', 12);

legend({'$\omega=0.0$', '$\omega=1.0$'}, 'interpreter', 'latex', 'position', [0.47 0.10 0.20 0.15], ...
    'fontname', 'kaiti', 'fontsize', 10.5);
legend('boxoff');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 讨论参数 {sigma eta varphi theta phipi o}
% 的变化对小国模型唯一均衡解的影响  
% 换句话讲，这些参数会影响 BK 条件
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; 
close all; 
clc;

% 参数值初始化
sig = 2.11; eta = 2.50; o = 0.37; tht = 0.75; 
bet = 0.99; vphi=3;     phipi=1.5;

% sigma 变动对唯一均衡解的影响
clear sig;
for ii = 1:51
    for jj = 1:51
        ome = (ii-1)/50;
        sig = 0.5+(jj-1)/10;
        Psi = eta*o*(2-o)*((sig+o)*ome - o)/((1-2*eta+eta*o)*o+(1-2*o+2*eta*o-eta*o^2)*ome);
        Phi = eta*(2-o)*(1-o+Psi)/(sig*eta*(2-o)+Psi);
        Gam = (1-tht)*(1-bet*tht) * (1 + vphi*Phi);
        if Gam < -2*tht*(1+bet)/(1+phipi) || Gam > 0
            sigcon(ii, jj) = 1;
        else
            sigcon(ii, jj) = 0;
        end
    end
end

% eta 变动对唯一均衡解的影响
clear ome eta Psi Phi Gam;
sig = 2.11;
for ii = 1:51
    for jj = 1:51
        ome = (ii-1)/50;
        eta = 1+(jj-1)/10;
        Psi = eta*o*(2-o)*((sig+o)*ome - o)/((1-2*eta+eta*o)*o+(1-2*o+2*eta*o-eta*o^2)*ome);
        Phi = eta*(2-o)*(1-o+Psi)/(sig*eta*(2-o)+Psi);
        Gam = (1-tht)*(1-bet*tht) * (1 + vphi*Phi);
        if Gam < -2*tht*(1+bet)/(1+phipi) || Gam > 0
            etacon(ii, jj) = 1;
        else
            etacon(ii, jj) = 0;
        end
    end
end

% theta 变动对唯一均衡解的影响
clear ome tht Psi Phi Gam;
eta = 2.50;
for ii = 1:51
    for jj = 1:51
        ome = (ii-1)/50;
        tht = (jj-1)/50;
        Psi = eta*o*(2-o)*((sig+o)*ome - o)/((1-2*eta+eta*o)*o+(1-2*o+2*eta*o-eta*o^2)*ome);
        Phi = eta*(2-o)*(1-o+Psi)/(sig*eta*(2-o)+Psi);
        Gam = (1-tht)*(1-bet*tht) * (1 + vphi*Phi);
        if Gam < -2*tht*(1+bet)/(1+phipi) || Gam > 0
            thtcon(ii, jj) = 1;
        else
            thtcon(ii, jj) = 0;
        end
    end
end

% varphi 变动对唯一均衡解的影响
clear ome vphi Psi Phi Gam;
tht = 0.75;
for ii = 1:51
    for jj = 1:51
        ome = (ii-1)/50;
        vphi = (jj-1)/5;
        Psi = eta*o*(2-o)*((sig+o)*ome - o)/((1-2*eta+eta*o)*o+(1-2*o+2*eta*o-eta*o^2)*ome);
        Phi = eta*(2-o)*(1-o+Psi)/(sig*eta*(2-o)+Psi);
        Gam = (1-tht)*(1-bet*tht) * (1 + vphi*Phi);
        if Gam < -2*tht*(1+bet)/(1+phipi) || Gam > 0
            vphcon(ii, jj) = 1;
        else
            vphcon(ii, jj) = 0;
        end
    end
end

% phipi 变动对唯一均衡解的影响
clear ome phipi Psi Phi Gam;
vphi = 3.00;
for ii = 1:50
    for jj = 1:51
        ome = (ii-1)/50;
        phipi = 1+jj/10;
        Psi = eta*o*(2-o)*((sig+o)*ome - o)/((1-2*eta+eta*o)*o+(1-2*o+2*eta*o-eta*o^2)*ome);
        Phi = eta*(2-o)*(1-o+Psi)/(sig*eta*(2-o)+Psi);
        Gam = (1-tht)*(1-bet*tht) * (1 + vphi*Phi);
        if Gam < -2*tht*(1+bet)/(1+phipi) || Gam > 0
            phipicon(ii, jj) = 1;
        else
            phipicon(ii, jj) = 0;
        end
    end
end

% o 变动对唯一均衡解的影响
clear ome o Psi Phi Gam;
phipi = 1.50;
for ii = 1:51
    for jj = 1:51
        ome = (ii-1)/50;
        o = (ii-1)/50;
        Psi = eta*o*(2-o)*((sig+o)*ome - o)/((1-2*eta+eta*o)*o+(1-2*o+2*eta*o-eta*o^2)*ome);
        Phi = eta*(2-o)*(1-o+Psi)/(sig*eta*(2-o)+Psi);
        Gam = (1-tht)*(1-bet*tht) * (1 + vphi*Phi);
        if Gam < -2*tht*(1+bet)/(1+phipi) || Gam > 0
            ocon(ii, jj) = 1;
        else
            ocon(ii, jj) = 0;
        end
    end
end

% 图8.3：各参数变动对唯一均衡解的影响
h0 = figure
set(h0,'Units','centimeters','Position',[2 2 14 18]);
subplot(3,2,1)
spy(sigcon, 'ko', 3);
ylabel('$\omega$', 'interpreter', 'latex', 'fontsize',12);
xlabel('$\sigma$', 'interpreter', 'latex', 'fontsize',12);
set(gca, 'xlim', [1, 51]);set(gca, 'xtick', 1:10:51);set(gca, 'xticklabel', 0.5:1:5.5);
set(gca, 'ylim', [1, 51]);set(gca, 'ytick', 1:10:51);set(gca, 'yticklabel', 0:0.2:1);

subplot(3,2,2)
spy(etacon, 'ko', 3);
ylabel('$\omega$', 'interpreter', 'latex', 'fontsize',12);
xlabel('$\eta$', 'interpreter', 'latex', 'fontsize',12);
set(gca, 'xlim', [1, 51]);set(gca, 'xtick', 1:10:51);set(gca, 'xticklabel', 1:1:6);
set(gca, 'ylim', [1, 51]);set(gca, 'ytick', 1:10:51);set(gca, 'yticklabel', 0:0.2:1);

subplot(3,2,3)
spy(thtcon, 'ko', 3);
ylabel('$\omega$', 'interpreter', 'latex', 'fontsize',12);
xlabel('$\theta$', 'interpreter', 'latex', 'fontsize',12);
set(gca, 'xlim', [1, 51]);set(gca, 'xtick', 1:10:51);set(gca, 'xticklabel', 0:0.2:1);
set(gca, 'ylim', [1, 51]);set(gca, 'ytick', 1:10:51);set(gca, 'yticklabel', 0:0.2:1);

subplot(3,2,4)
spy(vphcon, 'ko', 3);
xlabel('$\varphi$', 'interpreter', 'latex', 'fontsize',12);
ylabel('$\omega$', 'interpreter', 'latex', 'fontsize',12);
set(gca, 'xlim', [1, 51]);set(gca, 'xtick', 1:10:51);set(gca, 'xticklabel', 0:2:10);
set(gca, 'ylim', [1, 51]);set(gca, 'ytick', 1:10:51);set(gca, 'yticklabel', 0:0.2:1);

subplot(3,2,5)
spy(phipicon, 'ko', 3);
xlabel('$\phi_\pi$', 'interpreter', 'latex', 'fontsize',12);
ylabel('$\omega$', 'interpreter', 'latex', 'fontsize',12);
set(gca, 'xlim', [1, 50]);set(gca, 'xtick', 1:10:50);set(gca, 'xticklabel', 1.1:1:6);
set(gca, 'ylim', [1, 51]);set(gca, 'ytick', 1:10:51);set(gca, 'yticklabel', 0:0.2:1);

subplot(3,2,6)
spy(vphcon, 'ko', 3);
xlabel('$o$', 'interpreter', 'latex', 'fontsize',12);
ylabel('$\omega$', 'interpreter', 'latex', 'fontsize',12);
set(gca, 'xlim', [1, 51]);set(gca, 'xtick', 1:10:51);set(gca, 'xticklabel', 0:0.2:1);
set(gca, 'ylim', [1, 51]);set(gca, 'ytick', 1:10:51);set(gca, 'yticklabel', 0:0.2:1);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 不同参数（omega phiq）组合下求解小国开放经济模型
% （1）omega: 0:0.2:1;   
% （2）phiq: {-1, 0, 1};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; 
close all; 
clc;

% 循环求解
for ii = 1:6
    for jj = 1:3
        % 参数组合
        ome = (ii-1)/5; 
        phiq = jj-2;
        % 求解和模拟
        dynare small_country_nonlinear.mod;
        % 保存结果
        smo = ['sm', num2str(ii), num2str(jj)];
        assignin('base', smo, oo_);
        filepath = ['./result./smo', num2str(ii), num2str(jj), '.mat']
        save(filepath, smo);
    end
end

% 保存 M_ 文件
smm = M_;
save('./result./smm.mat', 'smm');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 脉冲反应图：
%（1）资本管制
%（2）有管理的浮动汇率制
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; 
close all; 
clc;

% 加载数据
load('./result./smo12')    % omega=0.0, phiq=0;
load('./result./smo22')    % omega=0.2, phiq=0;
load('./result./smo42')    % omega=0.6, phiq=0;
load('./result./smo52')    % omega=0.8, phiq=0;
load('./result./smo62')    % omega=1.0, phiq=0;
load('./result./smo61')    % omega=1.0, phiq=-1;
load('./result./smo63')    % omega=1.0, phiq=1;

% 脉冲反应期数
T = 10; 
t = 1:T;

%%%%%%%%%%%%%%%%%%%%
% （1）资本管制
%%%%%%%%%%%%%%%%%%%%
% （1.1）图8.4：本国生产率冲击
h0 = figure
set(h0,'Units','centimeters','Position',[2 2 14 8]);

subplot(2, 3, 1);
plot(t, sm12.irfs.y_etaa(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, sm42.irfs.y_etaa(t)*100, 'Color', [0.3,0.3,0.3], 'Linestyle', ':', 'LineWidth', 1);
plot(t, sm62.irfs.y_etaa(t)*100, 'Color', [0.5,0.5,0.5], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$y_t$', 'interpreter', 'latex', 'fontsize', 12);

subplot(2, 3, 2);
plot(t, sm12.irfs.ybar_etaa(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, sm42.irfs.ybar_etaa(t)*100, 'Color', [0.3,0.3,0.3], 'Linestyle', ':', 'LineWidth', 1);
plot(t, sm62.irfs.ybar_etaa(t)*100, 'Color', [0.5,0.5,0.5], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$\bar{y}_t$', 'interpreter', 'latex', 'fontsize', 12);

subplot(2, 3, 3);
plot(t, sm12.irfs.ygap_etaa(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, sm42.irfs.ygap_etaa(t)*100, 'Color', [0.3,0.3,0.3], 'Linestyle', ':', 'LineWidth', 1);
plot(t, sm62.irfs.ygap_etaa(t)*100, 'Color', [0.5,0.5,0.5], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$\tilde{y}_t$', 'interpreter', 'latex', 'fontsize', 12);

subplot(2, 3, 4)
plot(t, sm12.irfs.c_etaa(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, sm42.irfs.c_etaa(t)*100, 'Color', [0.3,0.3,0.3], 'Linestyle', ':', 'LineWidth', 1);
plot(t, sm62.irfs.c_etaa(t)*100, 'Color', [0.5,0.5,0.5], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$c_t$', 'interpreter', 'latex', 'fontsize', 12);

subplot(2, 3, 5)
plot(t, sm12.irfs.piH_etaa(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, sm42.irfs.piH_etaa(t)*100, 'Color', [0.3,0.3,0.3], 'Linestyle', ':', 'LineWidth', 1);
plot(t, sm62.irfs.piH_etaa(t)*100, 'Color', [0.5,0.5,0.5], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$\pi_{Ht}$', 'interpreter', 'latex', 'fontsize', 12);

subplot(2, 3, 6)
plot(t, sm12.irfs.q_etaa(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, sm42.irfs.q_etaa(t)*100, 'Color', [0.3,0.3,0.3], 'Linestyle', ':', 'LineWidth', 1);
plot(t, sm62.irfs.q_etaa(t)*100, 'Color', [0.5,0.5,0.5], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$q_t$', 'interpreter', 'latex', 'fontsize', 12);

legend({'$\omega=0.0$', '$\omega=0.6$', '$\omega=1.0$'}, 'Orientation', 'horizontal', ...
    'interpreter', 'latex', 'position', [0.45 0.02 0.20 0.04], ...
    'fontname', 'kaiti', 'fontsize', 10.5);
legend('boxoff');


% （1.2）图8.5：外国生产率冲击
h0 = figure
set(h0,'Units','centimeters','Position',[2 2 14 8]);
subplot(2, 3, 1);
plot(t, sm12.irfs.y_etaast(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, sm42.irfs.y_etaast(t)*100, 'Color', [0.3,0.3,0.3], 'Linestyle', ':', 'LineWidth', 1);
plot(t, sm62.irfs.y_etaast(t)*100, 'Color', [0.5,0.5,0.5], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$y_t$', 'interpreter', 'latex', 'fontsize', 12);

subplot(2, 3, 2);
plot(t, sm12.irfs.ybar_etaast(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, sm42.irfs.ybar_etaast(t)*100, 'Color', [0.3,0.3,0.3], 'Linestyle', ':', 'LineWidth', 1);
plot(t, sm62.irfs.ybar_etaast(t)*100, 'Color', [0.5,0.5,0.5], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$\bar{y}_t$', 'interpreter', 'latex', 'fontsize', 12);

subplot(2, 3, 3);
plot(t, sm12.irfs.ygap_etaast(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, sm42.irfs.ygap_etaast(t)*100, 'Color', [0.3,0.3,0.3], 'Linestyle', ':', 'LineWidth', 1);
plot(t, sm62.irfs.ygap_etaast(t)*100, 'Color', [0.5,0.5,0.5], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$\tilde{y}_t$', 'interpreter', 'latex', 'fontsize', 12);

subplot(2, 3, 4)
plot(t, sm12.irfs.c_etaast(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, sm42.irfs.c_etaast(t)*100, 'Color', [0.3,0.3,0.3], 'Linestyle', ':', 'LineWidth', 1);
plot(t, sm62.irfs.c_etaast(t)*100, 'Color', [0.5,0.5,0.5], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$c_t$', 'interpreter', 'latex', 'fontsize', 12);

subplot(2, 3, 5)
plot(t, sm12.irfs.piH_etaast(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, sm42.irfs.piH_etaast(t)*100, 'Color', [0.3,0.3,0.3], 'Linestyle', ':', 'LineWidth', 1);
plot(t, sm62.irfs.piH_etaast(t)*100, 'Color', [0.5,0.5,0.5], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$\pi_{Ht}$', 'interpreter', 'latex', 'fontsize', 12);

subplot(2, 3, 6)
plot(t, sm12.irfs.q_etaast(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, sm42.irfs.q_etaast(t)*100, 'Color', [0.3,0.3,0.3], 'Linestyle', ':', 'LineWidth', 1);
plot(t, sm62.irfs.q_etaast(t)*100, 'Color', [0.5,0.5,0.5], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$q_t$', 'interpreter', 'latex', 'fontsize', 12);

legend({'$\omega=0.0$', '$\omega=0.6$', '$\omega=1.0$'}, 'Orientation', 'horizontal', ...
    'interpreter', 'latex', 'position', [0.45 0.02 0.20 0.04], ...
    'fontname', 'kaiti', 'fontsize', 10.5);
legend('boxoff');


% 图8.6：本国产出的脉冲反应随资本管制程度的变化
h0 = figure
set(h0,'Units','centimeters','Position',[2 2 14 6]);

subplot(1, 2, 1)
plot(t, sm12.irfs.y_etaa(t)*100, 'Color', [0 0 0], 'Linestyle', ':', 'Linewidth', 1);
hold on
plot(t, sm22.irfs.y_etaa(t)*100, 'Color', [0.2 0.2 0.2], 'Linestyle', '-.', 'Linewidth', 1);
plot(t, sm42.irfs.y_etaa(t)*100, 'Color', [0.4 0.4 0.4], 'Linestyle', '--', 'Linewidth', 1);
plot(t, sm52.irfs.y_etaa(t)*100, 'Color', [0.6 0.6 0.6], 'Linestyle', '-', 'Linewidth', 1);
plot(t, sm62.irfs.y_etaa(t)*100, 'Color', [0.8 0.8 0.8], 'Linestyle', '-', 'Linewidth', 1);
hold off
xlim([0 T]); ylabel('%'); title('本国生产率冲击', 'fontsize', 12);

subplot(1, 2, 2)
plot(t, sm12.irfs.y_etaast(t)*100, 'Color', [0 0 0], 'Linestyle', ':', 'Linewidth', 1);
hold on
plot(t, sm22.irfs.y_etaast(t)*100, 'Color', [0.2 0.2 0.2], 'Linestyle', '-.', 'Linewidth', 1);
plot(t, sm42.irfs.y_etaast(t)*100, 'Color', [0.4 0.4 0.4], 'Linestyle', '--', 'Linewidth', 1);
plot(t, sm52.irfs.y_etaast(t)*100, 'Color', [0.6 0.6 0.6], 'Linestyle', '-', 'Linewidth', 1);
plot(t, sm62.irfs.y_etaast(t)*100, 'Color', [0.8 0.8 0.8], 'Linestyle', '-', 'Linewidth', 1);
hold off
xlim([0 T]); ylabel('%'); title('外国生产率冲击', 'fontsize', 12);

legend({'$\omega=0.0$', '$\omega=0.2$', '$\omega=0.6$', '$\omega=0.8$', '$\omega=1.0$'}, ...
    'Orientation', 'horizontal', 'interpreter', 'latex', ...
    'position', [0.35 0.01 0.20 0.04], 'fontsize', 12);
legend('boxoff');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% （2）货币政策对汇率作出反应
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% （2.1）图8.8：本国生产率冲击
h0 = figure
set(h0,'Units','centimeters','Position',[2 2 14 8]);

subplot(2, 3, 1);
plot(t, sm61.irfs.y_etaa(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, sm62.irfs.y_etaa(t)*100, 'Color', [0.3,0.3,0.3], 'Linestyle', ':', 'LineWidth', 1);
plot(t, sm63.irfs.y_etaa(t)*100, 'Color', [0.5,0.5,0.5], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$y_t$', 'interpreter', 'latex', 'fontsize', 12);

subplot(2, 3, 2);
plot(t, sm61.irfs.ygap_etaa(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, sm62.irfs.ygap_etaa(t)*100, 'Color', [0.3,0.3,0.3], 'Linestyle', ':', 'LineWidth', 1);
plot(t, sm63.irfs.ygap_etaa(t)*100, 'Color', [0.5,0.5,0.5], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$\tilde{y}_t$', 'interpreter', 'latex', 'fontsize', 12);

subplot(2, 3, 3)
plot(t, sm61.irfs.c_etaa(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, sm62.irfs.c_etaa(t)*100, 'Color', [0.3,0.3,0.3], 'Linestyle', ':', 'LineWidth', 1);
plot(t, sm63.irfs.c_etaa(t)*100, 'Color', [0.5,0.5,0.5], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$c_t$', 'interpreter', 'latex', 'fontsize', 12);

subplot(2, 3, 4)
plot(t, sm61.irfs.piH_etaa(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, sm62.irfs.piH_etaa(t)*100, 'Color', [0.3,0.3,0.3], 'Linestyle', ':', 'LineWidth', 1);
plot(t, sm63.irfs.piH_etaa(t)*100, 'Color', [0.5,0.5,0.5], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$\pi_{Ht}$', 'interpreter', 'latex', 'fontsize', 12);

subplot(2, 3, 5)
plot(t, sm61.irfs.pi_etaa(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, sm62.irfs.pi_etaa(t)*100, 'Color', [0.3,0.3,0.3], 'Linestyle', ':', 'LineWidth', 1);
plot(t, sm63.irfs.pi_etaa(t)*100, 'Color', [0.5,0.5,0.5], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$\pi_t$', 'interpreter', 'latex', 'fontsize', 12);

subplot(2, 3, 6)
plot(t, sm61.irfs.q_etaa(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, sm62.irfs.q_etaa(t)*100, 'Color', [0.3,0.3,0.3], 'Linestyle', ':', 'LineWidth', 1);
plot(t, sm63.irfs.q_etaa(t)*100, 'Color', [0.5,0.5,0.5], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$q_t$', 'interpreter', 'latex', 'fontsize', 12);

legend({'$\phi_q=-1$', '$\phi_q=0.0$', '$\phi_q=1.0$'}, 'interpreter', 'latex', 'position', [0.75 0.25 0.20 0.15], ...
    'fontname', 'kaiti', 'fontsize', 10.5);
legend('boxoff');

% （2.2）图8.9：外国生产率冲击
h0 = figure
set(h0,'Units','centimeters','Position',[2 2 14 8]);

subplot(2, 3, 1);
plot(t, sm61.irfs.y_etaast(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, sm62.irfs.y_etaast(t)*100, 'Color', [0.3,0.3,0.3], 'Linestyle', ':', 'LineWidth', 1);
plot(t, sm63.irfs.y_etaast(t)*100, 'Color', [0.5,0.5,0.5], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$y_t$', 'interpreter', 'latex', 'fontsize', 12);

subplot(2, 3, 2);
plot(t, sm61.irfs.ygap_etaast(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, sm62.irfs.ygap_etaast(t)*100, 'Color', [0.3,0.3,0.3], 'Linestyle', ':', 'LineWidth', 1);
plot(t, sm63.irfs.ygap_etaast(t)*100, 'Color', [0.5,0.5,0.5], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$\tilde{y}_t$', 'interpreter', 'latex', 'fontsize', 12);

subplot(2, 3, 3)
plot(t, sm61.irfs.c_etaast(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, sm62.irfs.c_etaast(t)*100, 'Color', [0.3,0.3,0.3], 'Linestyle', ':', 'LineWidth', 1);
plot(t, sm63.irfs.c_etaast(t)*100, 'Color', [0.5,0.5,0.5], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$c_t$', 'interpreter', 'latex', 'fontsize', 12);

subplot(2, 3, 4)
plot(t, sm61.irfs.piH_etaast(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, sm62.irfs.piH_etaast(t)*100, 'Color', [0.3,0.3,0.3], 'Linestyle', ':', 'LineWidth', 1);
plot(t, sm63.irfs.piH_etaast(t)*100, 'Color', [0.5,0.5,0.5], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$\pi_{Ht}$', 'interpreter', 'latex', 'fontsize', 12);

subplot(2, 3, 5)
plot(t, sm61.irfs.pi_etaast(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, sm62.irfs.pi_etaast(t)*100, 'Color', [0.3,0.3,0.3], 'Linestyle', ':', 'LineWidth', 1);
plot(t, sm63.irfs.pi_etaast(t)*100, 'Color', [0.5,0.5,0.5], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$\pi_t$', 'interpreter', 'latex', 'fontsize', 12);

subplot(2, 3, 6)
plot(t, sm61.irfs.q_etaast(t)*100, 'k-', 'LineWidth', 1);
hold on;
plot(t, sm62.irfs.q_etaast(t)*100, 'Color', [0.3,0.3,0.3], 'Linestyle', ':', 'LineWidth', 1);
plot(t, sm63.irfs.q_etaast(t)*100, 'Color', [0.5,0.5,0.5], 'Linestyle', '--', 'LineWidth', 1);
hold off;
xlim([0 T]); ylabel('%'); title('$q_t$', 'interpreter', 'latex', 'fontsize', 12);

legend({'$\phi_q=-1$', '$\phi_q=0.0$', '$\phi_q=1.0$'}, 'interpreter', 'latex', 'position', [0.75 0.12 0.20 0.15], ...
    'fontname', 'kaiti', 'fontsize', 10.5);
legend('boxoff');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 图8.7：资本管制对产出的影响
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; 
close all; 
clc;

% 参数值
sig = 2.11; eta = 2.50; o = 0.37; bet = 0.99; vphi=3; phipi=1.5;
ome = 0:0.05:1;
Psi = eta*o*(2-o)*((sig+o)*ome - o)./((1-2*eta+eta*o)*o+(1-2*o+2*eta*o-eta*o^2)*ome);
Xi = (1+vphi)*eta*(2-o)*(1-o+Psi)./(sig*eta*(2-o) + Psi + vphi*eta*(2-o)*(1 - o + Psi));
h0 = figure
set(h0,'Units','centimeters','Position',[2 2 14 6]);
subplot(1,2,1)
plot(ome(1:7), Xi(1:7), 'k-', 'LineWidth', 1);
xlim([0 1]); ylim([-1.5 5]); ylabel('%'); 
xlabel('$\omega$', 'interpreter', 'latex', 'fontsize', 12);
hold on
plot(ome(8:end), Xi(8:end), 'k-', 'LineWidth', 1);
plot([0.325, 0.325], [-1.5, 5], 'r--', 'LineWidth', 0.5);
text(0.3, -2.02, '$\omega^\ast$', 'interpreter', 'latex');
hold off
title('$\partial y_t/\partial a_t$', 'interpreter', 'latex', 'fontsize', 12);
subplot(1,2,2)
plot(ome(1:13), Psi(1:13), 'k-', 'LineWidth', 1);
xlim([0 1]); ylim([-24 22]);
xlabel('$\omega$', 'interpreter', 'latex', 'fontsize', 12);
hold on
plot(ome(15:end), Psi(15:end), 'k-', 'LineWidth', 1);
plot([0.65, 0.65], [-24, 22], 'r--', 'LineWidth', 0.5);
hold off
title('$\Psi$', 'interpreter', 'latex', 'fontsize', 12);
