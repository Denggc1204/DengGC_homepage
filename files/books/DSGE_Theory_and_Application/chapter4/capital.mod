%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 引入生产性资本后的RBC模型求解      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 作者: Deng, Guichuan               %
% 邮箱: denggch@mail.sysu.edu.cn     %
% 最新更新时间: 2026/10/12           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// ------------
// 声明内生变量
// ------------
var y c n a k i w r rk;

// -------------
//  声明外生变量
// -------------
varexo etaa;

// --------
// 声明参数
// --------
parameters bet sig phi eps alp del mu rho siga rss rkss Thet;
bet   =      0.99;
sig   =      2.10;
phi   =      3.00;
eps   =      10.0;
alp   =      0.40;
del   =      0.025;
mu    =      (eps-1)/eps;
// set_param_value 函数可以实现在 mod 文件外给参数赋值
set_param_value('rho', rho); 
siga  =      0.01;
rss   =      1/bet; 
rkss  =      rss+del; 
Thet  =      alp*mu/(rkss-1);

// --------
// 声明模型
// --------
model(linear);
sig * (c(+1) - c) = r;
phi * n + sig * c = w;
rkss * rk(+1) = rss*r;
k = (1 - del) * k(-1) + del * i; 
y = a + alp * k(-1) + (1-alp) * n;
w + n = y;
k(-1) + rkss/(rkss - 1)*rk= y;
y = (1-del*Thet) * c + del*Thet * i;
a = rho * a(-1) + etaa; 
end;

// ----------
// 声明初始值
// ----------
initval;
r     =   0; 
w     =   0; 
a     =   0; 
y     =   0; 
c     =   0;  
n     =   0;
k     =   0; 
i     =   0;  
rk    =   0;
end;

steady;
check;

shocks;
var etaa = siga^2;
end;

stoch_simul(irf = 60, order = 1, periods = 200000, drop = 400);
