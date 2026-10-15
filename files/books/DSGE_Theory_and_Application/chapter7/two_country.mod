%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 第7章：对称两国DSGE模型求解        
% 两种情形：                       
% （1）PCP；              
% （2）LCP；             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 作者: Deng, Guichuan             
% 邮箱: denggch@mail.sysu.edu.cn   
% 最新更新时间: 2026/10/15         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// 定义宏变量 case
// 如果 case = 1，表示PCP，默认值
// 如果 case = 2，表示LCP
@#ifndef case
    @#define case = 1
@#endif

// ------------
// 声明内生变量
// ------------
var 
// 弹性价格变量
yb ybst cb cbst pHb pFbst qb
// 黏性价格变量
c cst n nst y yst w wst r rst piH piFst pi pist pH pFst piHst piF
q s pHdot s1 s2 pFstdot x1st x2st vH vFst pHst pF  vHst vF de
@#if case == 2
   pHstdot s1st s2st pFdot x1 x2 cm
@#endif
a ast;

// -------------
//  声明外生变量
// -------------
varexo etaa etaast;

// --------
// 声明参数
// --------
parameters bet sig vph eps eta alp tht rhor phipi phiy
           rhoa rhoast siga sigast mu;
bet=0.99; 
sig=2.11; 
vph=3.00; 
eps=10.0;  
eta=2.50;
set_param_value('alp', alp);
tht=0.75; 
rhor=0.75;
phipi=1.50;
phiy=0.125;
rhoa=0.50;
siga=0.01;
rhoast=0.50;
sigast=0.01;
mu=(eps-1)/eps;

// --------
// 声明模型
// --------
model;
// ----------------
// 弹性价格系统
// ----------------
// equ(1-2): 消费-劳动权衡方程（劳动供给方程）
exp(sig*cb + vph*yb) = exp(pHb + (1+vph)*a);
exp(sig*cbst + vph*ybst) = exp(pFbst + (1+vph)*ast);

// equ(3-4): 产品市场均衡条件
exp(yb) = (1-alp) * exp(-eta*pHb + cb) + alp * exp(-eta*(pHb - qb) + cbst);
exp(ybst) = alp * exp(-eta*(pFbst + qb) + cb) + (1 - alp) * exp(-eta*pFbst + cbst);

// equ(5): 国际风险分担条件
exp(sig * (cb - cbst)) = exp(qb);

// equ(6-7): 价格指数方程
alp * exp((1 - eta) * (pFbst + qb)) + (1 - alp) * exp((1 - eta) * pHb) = 1;
alp *  exp((1 - eta) * (pHb - qb)) + (1 - alp) * exp((1 - eta) * pFbst) = 1;

// ---------------
//  黏性价格系统
// ---------------
// Equ(1-2): 消费-储蓄方程
exp(-sig * (c(+1)-c) + r - pi(+1)) = 1;
exp(-sig * (cst(+1)-cst) + rst - pist(+1)) = 1;

// Equ(3-4): 劳动供给方程
exp(w) = exp(sig*c + vph*n);
exp(wst) = exp(sig*cst + vph*nst);

// Equ(5): 国际风险分担条件
exp(-sig * (cst - c)) = exp(q);

// Equ(6-9): 本国生产价格设定
exp(pHdot - pH) = exp(s1 - s2);
exp(s1) = (1-bet*tht)*exp(-eta*pH + (1-sig)*c + w - a) + bet*tht*exp(eps*piH(+1) + s1(+1));
exp(s2) = (1-bet*tht)*exp((1 - eta)*pH + (1-sig)*c) + bet*tht * exp((eps-1)*piH(+1) + s2(+1));
tht*exp((eps-1)*piH) + (1-tht) * exp((1-eps)*(pHdot-pH)) = 1;

// Equ(10-13): 外国生产价格设定
exp(pFstdot - pFst) = exp(x1st - x2st);
exp(x1st) = (1-bet*tht)*exp(-eta*pFst + (1-sig)*cst+wst-ast) + bet*tht*exp(eps * piFst(+1)+x1st(+1));
exp(x2st) = (1-bet*tht)*exp((1-eta)*pFst + (1-sig)*cst) + bet*tht*exp((eps-1) * piFst(+1)+x2st(+1));
tht*exp((eps-1)*piFst) + (1-tht)*exp((1-eps)*(pFstdot-pFst)) = 1;

// Equ(14-17): 本国出口价格设定
@#if case == 1
   exp(pHst) = exp(pH - q);
@#endif
@#if case == 2
   exp(pHstdot - pHst) = exp(s1st - s2st);
   exp(s1st) = (1-bet*tht)*exp(-eta*pHst+cst-sig*c+w-a) + bet*tht*exp(eps*piHst(+1)+s1st(+1));
   exp(s2st) = (1-bet*tht)*exp((1-eta)*pHst+cst-sig*c+q) + bet*tht*exp((eps-1)*piHst(+1)+s2st(+1));
   tht*exp((eps-1)*piHst) + (1-tht)*exp((1-eps)*(pHstdot-pHst)) = 1;
@#endif

// Equ(18-21): 本国进口价格设定
@#if case == 1
   exp(pF) = exp(pFst + q);
@#endif
@#if case == 2
   exp(pFdot-pF) = exp(x1-x2);
   exp(x1) = (1-bet*tht)*exp(-eta*pF+c-sig*cst+wst-ast) + bet*tht*exp(eps*piF(+1)+x1(+1));
   exp(x2) =  (1-bet*tht)*exp((1-eta)*pF+c-sig*cst-q) + bet*tht*exp((eps-1)*piF(+1)+x2(+1));
   tht*exp((eps-1)*piF) + (1-tht)*exp((1-eps)*(pFdot-pF)) = 1;
@#endif

// Equ(22-23): 产品市场均衡条件
exp(y) = (1-alp)*exp(-eta*pH+c+vH) + alp*exp(-eta*pHst+cst+vHst);
exp(yst) = (1-alp)*exp(-eta*pFst+cst+vFst) + alp*exp(-eta*pF+c+vF);

// Equ(24-25) : 生产价格分散度
exp(vH) = tht*exp(eps*piH+vH(-1))+(1-tht)*exp(-eps*(pHdot-pH));
exp(vFst) = tht*exp(eps*piFst+vFst(-1)) + (1-tht)*exp(-eps*(pFstdot-pFst));

// Equ(26-27) : 贸易价格分散度
@#if case == 1
   exp(vHst) = exp(vH);
   exp(vF) = exp(vFst);
@#endif
@#if case == 2
   exp(vHst) = tht*exp(eps*piHst+vHst(-1))+(1-tht)*exp(-eps*(pHstdot-pHst));
   exp(vF) = tht*exp(eps*piF+vF(-1)) + (1-tht)*exp(-eps*(pFdot-pF));
@#endif

// Equ(28-29): 生产函数
exp(y) = exp(n + a);
exp(yst) = exp(ast + nst);

// Equ(30-31): 消费价格指数
alp*exp((1-eta)*pF) + (1-alp)*exp((1-eta)*pH) = 1;
alp*exp((1-eta)*pHst) + (1-alp)*exp((1-eta)*pFst) = 1;

// Equ(32-35): 实际价格定义方程
exp(piH) = exp(pi + pH - pH(-1));
exp(piFst) = exp(pist + pFst - pFst(-1));
exp(piHst) = exp(pist + pHst - pHst(-1));
exp(piF) = exp(pi + pF - pF(-1));

// Equ(36-37): 货币政策方程
r = rhor * r(-1) + (1-rhor) * (phiy * (y-yb) + phipi * piH);
rst = rhor * rst(-1) + (1-rhor) * (phiy * (yst - ybst) + phipi * piFst);

// Equ(38): 定义贸易条件
@#if case == 1
   exp(s) = exp(pFst - pH + q);
@#endif
@#if case == 2
   exp(s) = exp(pF - pH);
@#endif

// Equ(39): 定义货币错配
@#if case == 2
   exp(cm) = exp(q + pHst + pH);
@#endif

// Equ(40): 定义实际汇率
exp(q - q(-1)) = exp(de + pist - pi);

// Equ(41-42): 生产率冲击
a = rhoa * a(-1) + etaa;
ast = rhoast * ast(-1) + etaast;
end;

// ----------
// 声明初始值
// ----------
initval;
yb     = 0; 
ybst   = 0; 
cb     = 0; 
cbst   = 0; 
pHb    = 0; 
pFbst  = 0; 
qb     = 0; 
c      = 0; 
cst    = 0; 
n      = 0; 
nst    = 0; 
y      = 0; 
yst    = 0; 
w      = 0; 
wst    = 0;
r      = 0; 
rst    = 0; 
piH    = 0; 
piFst  = 0; 
pi     = 0; 
pist   = 0; 
pH     = 0; 
pFst   = 0; 
pHdot  = 0; 
s1     = 0; 
s2     = 0; 
pFstdot= 0; 
x1st   = 0; 
x2st   = 0; 
vH     = 0; 
vFst   = 0; 
pHst   = 0; 
pF     = 0; 
vHst   = 0; 
vF     = 0; 
s      = 0; 
q      = 0; 
a      = 0; 
ast    = 0; 
piHst  = 0; 
piF    = 0; 
de     = 0;
@#if case == 2
   pHstdot= 0; 
   s1st   = 0; 
   s2st   = 0; 
   pFdot  = 0; 
   x1     = 0; 
   x2     = 0; 
   cm     = 0;
@#endif
end;

steady;
check;

shocks;
var etaa = siga^2;
var etaast = sigast^2;
end;

stoch_simul(irf = 100, order = 1, noprint, nograph, periods = 0);