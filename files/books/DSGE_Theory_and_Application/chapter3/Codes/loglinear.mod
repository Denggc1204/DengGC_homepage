/////////////////////////////////////////////////
// 第3章: 该Mod文件用于求解基准RBC模型         //
// 对数线性化形式且一阶近似                    //
// Mod文件的模型与动态系统中的变量一致         //
/////////////////////////////////////////////////
// 作者: Deng, Guichuan                        //
// 邮箱: denggch@mail.sysu.edu.cn              //
// 最新更新时间: 2026/10/11                    //
/////////////////////////////////////////////////

// -------------
// 声明内生变量
// -------------
var y c n r w a;

// ------------
// 声明外生变量
// ------------
varexo etaa;

// --------------
// 声明参数并赋值
// --------------
parameters bet sig phi eps mu rho siga;
bet = 0.99; 
sig = 2.10;
phi = 3.00; 
eps = 10.0;    
mu  = (eps-1)/eps; 
rho = 0.88; 
siga= 0.01; 

// --------
// 声明模型
// --------
model(linear);
-sig * c(+1) + sig * c = - r;
phi * n + sig * c = w; 
w = a; 
y = a + n;
y = c;
a = rho*a(-1) + etaa;
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
end;

steady;
check;

shocks;
var etaa = siga^2;
end;

stoch_simul(order=1, irf=40, periods = 200000);

