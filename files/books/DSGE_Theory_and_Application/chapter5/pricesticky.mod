%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 引入价格黏性的DSGE模型求解         %
% 两种情形：                         %
% （1）无指数化调整；                %
% （2）有指数化调整；                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 作者: Deng, Guichuan               %
% 邮箱: denggch@mail.sysu.edu.cn     %
% 最新更新时间: 2026/10/12           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// 定义宏变量 index，区分有指数化调整和没有指数分调整两种情况
// 当 index == 0, 表示模型无指数化调整，此为默认值
// 当 index == 1, 表示模型有指数化调整
@#ifndef index
    @#define index = 0
@#endif

// ------------
// 声明内生变量
// ------------
var y c n mc pi r w a ybar ytilde;

// -------------
//  声明外生变量
// -------------
varexo etaa;

// --------
// 声明参数
// --------
parameters bet sig phi eps tht gap phipi mu rho siga;
bet        =      0.99;
sig        =      2.10;
phi        =      3.00;
eps        =      10.0;
tht        =      0.75;
set_param_value('gap', gap);
set_param_value('phipi', phipi);
mu         =      (eps-1)/eps;
rho        =      0.88;
siga       =      0.01;

// --------
// 声明模型
// --------
model(linear);
sig*(c(+1) - c) + pi(+1) = r;
phi*n + sig*c = w;
mc = w - a;
y = a + n;
y = c;
@#if index == 0
   // 无指数化调整的情形
   pi = bet*pi(+1) + (1-tht)*(1-bet*tht)/tht*mc;
@#else
   // 存在指数化调整的情形
   (1+bet*gap)*pi = bet*pi(+1) + (1-tht)*(1-bet*tht)/tht*mc + gap*pi(-1);
@#endif
r = phipi * pi;
a = rho*a(-1) + etaa; 
ybar = (1+phi)/(sig+phi) * a; 
ytilde = y - ybar; 
end;

// ----------
// 声明初始值
// ----------
initval;
y     =   0; 
c     =   0;  
n     =   0;
pi    =   0;
mc    =   0;
r     =   0; 
w     =   0; 
a     =   0; 
end;

steady;
check;

shocks;
var etaa = siga^2;
end;

stoch_simul(irf = 100, order = 1, noprint, nograph, periods = 20000, drop = 400);