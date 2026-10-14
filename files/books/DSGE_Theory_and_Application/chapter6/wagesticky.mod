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

// 定义宏变量 index
// 如果 index = 1，表示有指数化调整
// 如果 index = 0，表示无指数化调整，默认值
@#ifndef index
    @#define index = 0
@#endif

// ------------
// 声明内生变量
// ------------
var y c n w pi piw r a ybar ygap;

// -------------
//  声明外生变量
// -------------
varexo etaa;

// --------
// 声明参数
// --------
parameters bet sig vph epw thtw gaw phipi kapw rho siga;

bet          =      0.99;
sig          =      2.10;
vph          =      3.00;
epw          =      10.0;
set_param_value('thtw', thtw);
gaw          =      0.58;
phipi        =      1.50;
kapw         =      (1-thtw)*(1-bet*thtw)/thtw;
rho          =      0.88;
siga         =      0.01;

// --------
// 声明模型
// --------
model(linear);
r = sig * (c(+1) - c) + pi(+1);
@#if index == 0
   piw = kapw/(1+vph*epw) * (sig * c + vph * n - w) + bet*piw(+1);
@#else
   piw = kapw/(1+vph*epw) * (sig * c + vph * n - w) + bet*piw(+1) - gaw*(bet*pi-pi(-1));
@#endif
piw = w - w(-1) + pi;
y = a + n;
w = a;
y = c;
r = phipi * pi;
a = rho * a(-1) + etaa; 
ybar = (1+vph)/(sig+vph)*a;
ygap = y - ybar;
end;

// ----------
// 声明初始值
// ----------
initval;
y     =   0; 
c     =   0;  
n     =   0;
w     =   0;
pi    =   0;
piw   =   0;
r     =   0; 
w     =   0; 
a     =   0; 
end;

steady;
check;

shocks;
var etaa = siga^2;
end;

stoch_simul(irf = 40, order = 1, nograph, noprint, periods = 0);
