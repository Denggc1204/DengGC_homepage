/////////////////////////////////////////////////
// 第3章: 该Mod文件用于求解基准RBC模型         //
// 原始方程形式且二阶近似                      //
// Mod文件的模型与动态系统中的变量一致         //
/////////////////////////////////////////////////
// 作者: Deng, Guichuan                        //
// 邮箱: denggch@mail.sysu.edu.cn              //
// 最新更新时间: 2026/10/11                    //
/////////////////////////////////////////////////

// -------------
// 声明内生变量
// -------------
var Y C N R W A U;

// ------------
// 声明外生变量
// ------------
varexo etaa;

// --------------
// 声明参数并赋值
// --------------
parameters bet sig phi eps mu rho siga Uss;
bet    =      0.99;
sig    =      2.10;
phi    =      3.00;
eps    =      10;
mu     =      (eps-1)/eps;
rho    =      0.88;
siga   =      0.01;
Uss = (mu^((1-sig)/(sig+phi))/(1-sig) - mu^((1+phi)/(sig+phi))/(1+phi))/(1-bet);

// --------
// 声明模型
// --------
model;
bet*C(+1)^(-sig)/(C^(-sig)) = 1/R;
N^phi * C^sig = W; 
W = mu*A; 
Y = A * N;
Y = C;
log(A)= rho * log(A(-1)) + etaa;
U = C^(1-sig)/(1-sig) - N^(1+phi)/(1+phi) + bet*U(+1);
end;

// ----------
// 声明初始值
// ----------
initval;
A     =   1; 
R     =   1/bet; 
W     =   mu*A; 
Y     =   mu^(1/(sig+phi)); 
C     =   mu^(1/(sig+phi));  
N     =   mu^(1/(sig+phi));
U     =   Uss;
end;

steady;
check;

shocks;
var etaa = siga^2;
end;

stoch_simul(irf = 100, order = 2, noprint, nograph, periods = 200000);