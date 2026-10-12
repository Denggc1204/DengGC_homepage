/////////////////////////////////////////////////
// 第3章: 该Mod文件用于求解基准RBC模型         //
// 非零稳态指数形式且二阶近似                  //
// Mod文件的模型与动态系统中的变量一致         //
/////////////////////////////////////////////////
// 作者: Deng, Guichuan                        //
// 邮箱: denggch@mail.sysu.edu.cn              //
// 最新更新时间: 2026/10/11                    //
/////////////////////////////////////////////////

// -------------
// 声明内生变量
// -------------
var lnY lnC lnN lnR lnW lnA U; 

// ------------
// 声明外生变量
// ------------
varexo etaa;

// --------------
// 声明参数并赋值
// --------------
parameters bet sig phi eps mu rho siga Uss;
bet = 0.99; 
sig = 2.10;
phi = 3.00; 
eps = 10.0;    
mu  = (eps-1)/eps; 
rho = 0.88; 
siga= 0.01; 
Uss = (mu^((1-sig)/(sig+phi))/(1-sig) - mu^((1+phi)/(sig+phi))/(1+phi))/(1-bet);

// --------
// 声明模型
// --------
model;
bet*exp(-sig*(lnC(+1) - lnC)) = exp(-lnR);
exp(phi*lnN + sig*lnC)= exp(lnW);
exp(lnW) = mu*exp(lnA);
exp(lnY) = exp(lnA + lnN);
exp(lnY) = exp(lnC);
lnA = rho*lnA(-1) + etaa;
U = exp((1-sig)*lnC)/(1-sig) - exp((1+phi)*lnN)/(1+phi) + bet*U(+1);
end;

// ----------
// 声明初始值
// ----------
initval;
lnA = 0;
lnR = -log(bet);
lnW = log(mu);
lnY = 1/(sig+phi)*log(mu);
lnC = 1/(sig+phi)*log(mu);
lnN = 1/(sig+phi)*log(mu);
U   = Uss;
end;

steady;
check;

shocks;
var etaa = siga^2;
end;

stoch_simul(order=2, irf=100, noprint, nograph, periods = 200000);
