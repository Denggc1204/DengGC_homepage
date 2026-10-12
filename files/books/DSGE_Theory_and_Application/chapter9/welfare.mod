%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 一阶近似求解标准NK模型，并计算福利水平
% 比较四种货币政策规则福利效果：
% （1）通胀目标制
% （2）标准泰勒规则
% （3）带利率平滑的泰勒规则
% （4）对产出作出反应的泰勒规则
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 作者: Deng, Guichuan               
% 邮箱: denggch@mail.sysu.edu.cn    
% 最新更新时间: 2026/10/13           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// 定义宏变量向量 policy
@#define policy = ["r = phipi * pi;", \\
                   "r = phiy * ygap + phipi * pi;", \\
                   "r = rhor * r(-1) + (1-rhor) * (phiy * ygap + phipi * pi);",\\
                   "r = phiy * y + phipi * pi;"]

// 定义宏变量 SEC, 如果 SEC == 0 则一阶近似，默认值
@#ifndef SEC
    @#define SEC = 0
@#endif

// ------------------
// 声明内生变量
// ------------------
var c n y yi w r pi p s1 s2 v a ygap wel;

// ------------------
//  声明外生变量
// ------------------
varexo eta;

// ------------------
//  声明参数
// ------------------
parameters bet vphi eps tht rhor phiy phipi rhoa siga;

// --------------
//  校准参数值
// --------------
param = paramfun;
bet      =      param(1);
vphi     =      param(2);
eps      =      param(3);
tht      =      param(4);
rhor     =      param(5);
phiy     =      param(6);
phipi    =      param(7);
rhoa     =      param(8);
siga     =      param(9);


// ---------------
//  声明模型
// ---------------
model;

// Welfare
wel = c - exp((1 + vphi) * n)/(1+vphi) + bet * wel(+1);

// equ(1): Euler Equation
exp(c - c(+1) - pi(+1) + r) = 1;

// equ(2): Labor supply
exp(c + vphi * n - w) = 1;

// equ(3): Production
exp(yi - n - a) = 1;

// equ(4): price
exp(p - s1 + s2) = 1;

// equ(5): s1
exp(s1) = (1-bet*tht)*exp(w - a + y - c) + bet*tht * exp(eps*pi(+1) + s1(+1));

// equ(6): s2
exp(s2) = (1-bet*tht)*exp(y - c) + bet*tht * exp((eps-1)*pi(+1) + s2(+1));

// equ(7): price index
tht*exp((eps-1)*pi) + (1-tht)*exp((1-eps)*p) = 1;

//equ(8): final goods market
exp(y) = exp(c);

//equ(9): intermediate goods market
exp(yi) = exp(y + v);

//equ(10): price disperse
exp(v) = tht * exp(eps*pi + v(-1)) + (1-tht)*exp(-eps*p);

// equ(11): monetary policy
@{policy[kk]}

// equ(12): output gap
exp(ygap) = exp(y - a);

// shocks
a = rhoa * a(-1) + eta;
end;

// ----------
// 声明初始值
// ----------

initval;
wel    =   - 1/((1-bet)*(1+vphi));
c      =   0;
n      =   0; 
y      =   0; 
yi     =   0; 
pi     =   0; 
w      =   0; 
r      =   0; 
p      =   0; 
s1     =   0; 
s2     =   0;
v      =   0; 
a      =   0;
ygap   =   0;
end;

steady;
check;

shocks;
var eta = siga^2;
end;

@#if SEC == 0
   stoch_simul(irf = 40, order = 1, nograph, noprint, periods = 100000, drop = 400);
@#else
   stoch_simul(irf = 40, order = 2, nograph, noprint, pruning, periods = 100000, drop = 400);
@#endif