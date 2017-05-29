(* ::Package:: *)

(* ::Input:: *)
(*L = {{c11,c12},{c21, -c11}};*)
(*J = {{0,-1},{1,0}};*)
(*X = {{x/y,1/y},{-x^2/y-y,-x/y}};*)
(*HLie := Tr[L.X] + 3/2*Tr[J.X];*)
(*a = (u1-u2)/sqrt[3];*)
(*b = u0/3 - 2u1/3-2u2/3;*)
(*c= u0;*)
(*f1=(y*(b+2a*x-c*x^2+c*y^2))/(b+2*a*x-c*x^2-c*y^2);*)
(*f2 = (2*(a-c*x)*y^2)/(b+2*a*x-c*x^2-c*y^2);*)
(*Hh = v1*f1 + v2*f2;*)
(*H = Hlie + Hh;*)
(*D[H,x] ;*)
(*D[H,y]; *)
(*D[H,c11];*)
(*D[H,c12];*)
(*D[H,c21];*)
(*D[H,v1];*)
(*D[H,v2]; *)
(**)
(*D[H,x]/.{u0 -> 0, u1 -> 0, u2 -> 1};*)


(* ::InheritFromParent:: *)
(**)


(* ::Input:: *)
(*D[H,y]/.{u0 -> 0, u1 -> 1, u2 -> 0}*)
(*D[H,y]/.{u0 -> 1, u1 -> 0, u2 -> 0}*)


(* ::Input:: *)
(*D[H,c11]/.{u0 -> 0, u1 -> 0, u2 -> 1};*)
(*D[H,c11]/.{u0 -> 0, u1 -> 1, u2 -> 0};*)
(*D[H,c11]/.{u0 -> 1, u1 -> 0, u2 -> 0};*)
(*D[H,c12]/.{u0 -> 0, u1 -> 0, u2 -> 1};*)
(*D[H,c12]/.{u0 -> 0, u1 -> 1, u2 -> 0};*)
(*D[H,c12]/.{u0 -> 1, u1 -> 0, u2 -> 0};*)
(*D[H,c21]/.{u0 -> 0, u1 -> 0, u2 -> 1};*)
(*D[H,c21]/.{u0 -> 0, u1 -> 1, u2 -> 0};*)
(*D[H,c21]/.{u0 -> 1, u1 -> 0, u2 -> 0};*)
(*D[H,v1]/.{u0 -> 0, u1 -> 0, u2 -> 1};*)
(*D[H,v1]/.{u0 -> 0, u1 -> 1, u2 -> 0};*)
(*D[H,v1]/.{u0 -> 1, u1 -> 0, u2 -> 0};*)
(*D[H,v2]/.{u0 -> 0, u1 -> 0, u2 -> 1};*)
(*D[H,v2]/.{u0 -> 0, u1 -> 1, u2 -> 0};*)
(*D[H,v2]/.{u0 -> 1, u1 -> 0, u2 -> 0};*)


(* Line ? - ? are ported from https://github.com/flyspeck/publications-of-thomas-hales/blob/master/geometry/reinhardt-optimal-control/reinhardt_mathematica.m
 *)


Id[x_]:=x;
I2={{1,0},{0,1}};
A2={{a,b},{c,-a}};(*X*)frac[{{a_,b_},{c_,d_}},z_]:=(a z+b)/(c z+d);
Cx[{x_,y_}]:=x+I y;

Ctimes[{x1_,y1_},{x2_,y2_}]:={(x1 x2-y1 y2),(x1 y2+y1 x2)};
todisk[z_]:=(z-I)/(z+I);
Tx=Transpose;
Lie[X_,Y_]:=X.Y-Y.X;
R1={{1/2,Sqrt[3]/2},{-Sqrt[3]/2,1/2}};
R=Tx[R1];
powR[i_]:={{Cos[i Pi/3],-Sin[i Pi/3]},{Sin[i Pi/3],Cos[i Pi/3]}};
ExpJt={{Cos[t],-Sin[t]},{Sin[t],Cos[t]}};

v0={1,0};v1=R.v0;v2=R.v1;v3=R.v2;v4=R.v3;v5=R.v4;R.v5;

Clear[conj,real,imag,Reps];
test=Module[{z,conj,real,imag,u},conj[u_]:=u/;Head[u]==Symbol;
conj[u_]:=u/;Head[u]==Integer;
conj[Rational[a_,b_]]:=Rational[a,b];
conj[a_Real]:=a;
conj[Power[a_,b_]]:=Power[conj[a],b];
conj[Plus[a_,b_]]:=Plus[conj[a],conj[b]];
conj[Times[a_,b_]]:=Times[conj[a],conj[b]];
conj[Complex[a_,b_]]:=Complex[a,-b];
z=frac[{{a,b},{c,d}},Cx[{x,y}]];
real=(conj[z]+z)/2;
imag=(z-conj[z])/(2 I);
{real,imag}//Simplify//InputForm];

(*Reps[frac[{{a,b},{c,d}},Cx[{x,y}]]]//InputForm;*)

Repfrac[{{a_,b_},{c_,d_}},{x_,y_}]:={(b*(d+c*x)+a*(d*x+c*(x^2+y^2)))/(d^2+2*c*d*x+c^2*(x^2+y^2)),(-(b*c*y)+a*d*y)/(d^2+2*c*d*x+c^2*(x^2+y^2))};

polyapprox[f_,n_]:=Apply[Plus,Table[t^i SeriesCoefficient[f,{t,0,i}],{i,0,n}]];

adjoint[h_]:=-{D[h,x],D[h,y]}//Simplify;

(*hyperbolic plane,state eqn.*)

Xupper={{x/y,-x^2/y-y},{1/y,-x/y}};
(*upper half-plane*)
zofX[X_]:=Module[{x,y},y=1/X[[2,1]]//Simplify;
x=y*X[[1,1]]//Simplify;
{x,y}];

(*get a basis of sl2 from Xupper and its derivatives*)

Module[{stdbasis},stdbasis[X_]:={(X[[1,1]]-X[[2,2]])/2,X[[1,2]],X[[2,1]]};
Det[{stdbasis[Xupper],stdbasis[D[Xupper,x]],stdbasis[D[Xupper,y]]}]];

{Det[Xupper],Tr[Xupper.Xupper]}//Simplify;
UpperHalfODE=Module[{Yhh,Dhh,rhh,xh,yh,t,Dx,Dy,xsub,Dxsub},Yhh=Xupper/.{x->xh[t],y->yh[t]};
xsub={xh[t]->x,yh[t]->y};
Dxsub={xh'[t]->Dx,yh'[t]->Dy};
Dhh=(D[Yhh,t]/.Dxsub)/.xsub;
rhh=-2/Tr[A2.Xupper];
{Dx,Dy}/.Solve[Dhh==(Xupper.(rhh (A2)-Xupper)),{Dx,Dy}][[1]]/.xsub//Factor]

Clear[Hpayoff];
Hampayoff=lambda0 (3/2) (x^2+y^2+1)/y;
adjoint[Hampayoff];

(*general upper-half function*)

Hupper={lambda1,lambda2}.UpperHalfODE;
HamLie=Tr[{{aL,bL},{cL,-aL}}.Xupper];
Hamiltonian=Hupper+Hampayoff+HamLie//Simplify;

(*switching set,bang-bang*)

z001sub={a->-1/Sqrt[3],b->-2/3,c->0};
z010sub={a->1/Sqrt[3],b->-2/3,c->0};
z100sub={a->0,b->1/3,c->1};
z011sub={a->-m,b->-2/3,c->0};
Zsub={a->(u1-u2)/Sqrt[3],b->u0/3-2 u1/3-2 u2/3,c->u0};

u001sub={u0->0,u1->0,u2->1};
u010sub={u0->0,u1->1,u2->0};
u100sub={u0->1,u1->0,u2->0};

m001sub={m->1/Sqrt[3]};
m010sub={m->-1/Sqrt[3]};

Rm001=1/2 {{1,-1/m},{1/m,1}};
test=(Rm001/.m001sub)-R;

switch[i_,j_]:=Module[{zsub},zsub={z100sub,z010sub,z001sub};(Hupper/.zsub[[i]])-(Hupper/.zsub[[j]])//Simplify];

switch32=switch[3,2];
switch12=switch[1,2];
switch13=switch[1,3];

(*Bang-Bang 001,010*)

(*001,010 bang-bang sol of state eqn*)

Clear[alpha,u,c0,msub,hamc];
xy011sub={x->-m+c0 u,y->c0 alpha u};
(*x0y0sub=*)

c0alpha011sub=Solve[({x,y}/.(xy011sub/.u->1))=={x0,y0},{c0,alpha}][[1]];
usub={u->E^(alpha t)};
StateSolG=Module[{G2,XXh,eqnn,eqn,gsol,a,b,d,c,i,j},G2={{a[u],b[u]},{c[u],d[u]}};
XXh=(Xupper/(alpha u))/.xy011sub//Factor;
eqnn={D[G2,u],G2.XXh}//Simplify;
eqn[i_,j_]:=eqnn[[1,i,j]]==eqnn[[2,i,j]];
{{a[u],b[u]},{c[u],d[u]}}/.DSolve[{eqn[1,1],eqn[1,2],eqn[2,1],eqn[2,2],a[1]==1,b[1]==0,c[1]==0,d[1]==1},{a,b,c,d},u][[1]]];
(StateSolG-I2) alpha^2 c0 u/(u-1)//Factor;
Det[StateSolG]//Simplify

(*cost*)
cost011=(3/2) Integrate[(x^2+y^2+1)/y/.xy011sub/.usub,{t,0,t1}];

costTransform=Module[{Rt,x1,y1},Rt={{a,b},{-b,a}};
{x1,y1}=Repfrac[Rt,{x,y}];
(x1^2+y1^2+1)/y1//Simplify]

(*adjoint solution:Lie term 001,010*)

AdjointSolLie011=Module[{Ytt,Rs,Cs,eqns,a1,a2,a3},Ytt=Xupper/.xy011sub;
Rs={{a1[u],a2[u]},{a3[u],-a1[u]}};
Cs=(Rs.Ytt-Ytt.Rs)/(alpha u)//Simplify;
eqns[i_,j_]:=D[Rs,u][[i,j]]==Cs[[i,j]];
(*Rs/.*)Rs/.DSolve[{eqns[1,1],eqns[1,2],eqns[2,1]},{a1,a2,a3},u][[1]]];

(*How C[1],C[2],C[3] relate to initial {{aL,bL},{cL,-aL}}*)

AdjointSolLie011Init=AdjointSolLie011/.{u->1};

(*solve ODE 001,010,main term*)


HamLie011=Tr[AdjointSolLie011.Xupper]//Simplify;

UpperHalfODE/.{c->0}//Simplify;
Hupper011=lambda1 y+lambda2 y^2/(m+x);
Hamiltonian011=Hupper011+Hampayoff+HamLie011/.xy011sub//Simplify;

AdjointODE011=adjoint[HamLie011]+adjoint[Hampayoff]+adjoint[Hupper011]/.xy011sub//Simplify;



AdjointSol011=Module[{adj},adj=AdjointODE011/.{lambda1->lambda1[t],lambda2->lambda2[t]}/.usub;{lambda1[t],lambda2[t]}/.DSolve[{lambda1'[t]==adj[[1]],lambda2'[t]==adj[[2]]},{lambda1,lambda2},t][[1]]//Simplify];


ham011c=Hamiltonian011/.{lambda1->AdjointSol011[[1]],lambda2->AdjointSol011[[2]]}/.usub//Simplify//Factor;Clear[Rt,aLbar,bLbar,cLbar,junk,a,b,c,aL,bL,cL,zbar,fbar,blam,Dblam,DblamODE];

(*general rotation*)
Rt={{aR,bR},{-bR,aR}};
(*Lambda transform*)
{{aLbar,bLbar},{cLbar,junk}}=Rt.{{aL,bL},{cL,-aL}}.Inverse[Rt];
(*Z0 transform*)
{{aZ0bar,bZ0bar},{cZ0bar,junk}}=Rt.{{a,b},{c,-a}}.Inverse[Rt];

(*z={x,y},f={f1,f2} transform*)
{zbar,fbar}=Module[{z,zb,fb,sub},sub={x'[0]->f1,y'[0]->f2,x[0]->x,y[0]->y};
z={x[t],y[t]};
fb=D[Repfrac[Rt,z],t];
zb=Repfrac[Rt,z];
{zb,fb}/.t->0/.sub//Simplify];
(*{lambda1,lambda2} transform*)
blam=Module[{blam1,blam2},{blam1,blam2}/.Solve[({blam1,blam2}.fbar=={lam1,lam2}.{f1,f2})/.{f1->{1,0},f2->{0,1}},{blam1,blam2}]][[1]]//Simplify

test=blam.fbar//Simplify

Dblam=Module[{vv},vv=D[(blam/.{lam1->lam1[t],lam2->lam2[t],x->x[t],y->y[t]}),t];
vv/.{x'[t]->UpperHalfODE[[1]],y'[t]->UpperHalfODE[[2]],x[t]->x,y[t]->y,lam1'[t]->-D[Hamiltonian,x],lam2'[t]->-D[Hamiltonian,y],lam1[t]->lambda1,lam2[t]->lambda2}//Simplify];

DblamODE=Module[{Hambar,vv,xB,yB,aB,bB,cB,laB1,laB2,aLB,bLB,cLB},Hambar=Hamiltonian/.{x->xB,y->yB,lambda1->laB1,lambda2->laB2,aL->aLB,bL->bLB,cL->cLB,a->aB,b->bB,c->cB};
vv=-{D[Hambar,xB],D[Hambar,yB]}//Simplify;
vv/.{xB->zbar[[1]],yB->zbar[[2]],laB1->blam[[1]],laB2->blam[[2]],aB->aZ0bar,bB->bZ0bar,cB->cZ0bar,aLB->aLbar,bLB->bLbar,cLB->cLbar}/.{lam1->lambda1,lam2->lambda2}//Simplify];
Dblam-DblamODE//Simplify

(*Fillipov*)

(*Filippov*)
(*vertices distinct*)

test=(UpperHalfODE/.z001sub)-(UpperHalfODE/.z010sub)//Together//Simplify;
(*fixed sign*)
test=Module[{p1,p2,rssub},p1=(UpperHalfODE/.z001sub)//Simplify;
p2=(UpperHalfODE/.z010sub)//Simplify;
rssub=Solve[{{r,s}.p1==1,{r,s}.p2==1},{r,s}][[1]]//Simplify;
UpperHalfODE.{r,s}-1/.Zsub/.rssub/.{u0->1-u1-u2}//Simplify]

(*Smoothed Octagon 001 Parameters*)

(*smoothed octagon state eqns*)



Clear[a,b,c,d,tc,sc,cx,fj,s,ggt,t1,Ts];

(*This entire section is summarized by the constants:cmoctsub,alphaoctsub,tmaxoct*)

cmoctsub={c0->1/Sqrt[3],m->1/Sqrt[3]};
alphaoctsub={alpha->Sqrt[(2 Sqrt[2]-1)]};
tmaxoct=(1/(2 alpha)) Log[2]/.alphaoctsub;
umaxoct={u->Sqrt[2]};

(*curve in SL2,unnormalized t1*)

goct1=Module[{G2,u0,u2,a,b,c,d,gt,gamma0,gamma2,strel,dd,ssub,s0},s0=1-1/Sqrt[2];
G2={{a,b},{c,d}};
u0={1,0};u2={-1/2,Sqrt[3]/2};
gamma0:={1,0}+{0,t1};
gamma2:={-1/Sqrt[2],1/Sqrt[2]}+{-1/Sqrt[2],-1/Sqrt[2]} (s-s0);
dd=Det[{gamma0,gamma2}];
strel=(dd-(dd/.{s->0,t1->0}));
ssub=Solve[strel==0,s][[1]];
G2/.Solve[G2.Tx[{u0,u2}]==Tx[({gamma0,gamma2})/.ssub],{a,b,c,d}][[1]]//Factor];

smoothedoctdensity=8 Sqrt[3] (8-8 Sqrt[2]+Sqrt[2] Log[2])/(4 (-4+Sqrt[2]))/Sqrt[12];

test=Module[{Yt1,IYt1},Yt1=Inverse[goct1].D[goct1,t1]//Simplify;
{"Yt1",Tr[Yt1],Det[Yt1]}//Simplify;
(*payoff*)IYt1=(-3/2) Integrate[Yt1[[1,2]]-Yt1[[2,1]],{t1,0,s0}];
(*from paper*){4 IYt1/Sqrt[12],smoothedoctdensity}//N];


(*"Xoct" unit speed*)
goct=Module[{t1sub},t1sub={t1->1-E^(-alpha t)};
goct1/.t1sub//Simplify];



Xoct=Inverse[goct].D[goct,t]//Simplify;

(*tests*)
Xoct+Inverse[Xoct].D[Xoct,t]//Simplify;
({"Xoct",Tr[Xoct],Det[Xoct]}//Simplify)/.alphaoctsub;
{"X terminal",Module[{Xoct0,XoctfN},Xoct0=Xoct/.{t->0}/.alphaoctsub;
XoctfN=(Xoct/.{t->tmaxoct})/.alphaoctsub;
R1.Xoct0.Tx[R1]-XoctfN]//Simplify}//N;

Module[{IXoct},IXoct=(-3/2) Integrate[(Xoct[[1,2]]-Xoct[[2,1]])/.alphaoctsub,{t,0,tmaxoct}];
{4 IXoct/Sqrt[12],smoothedoctdensity}//Simplify]//N;


zoctupper=zofX[Xoct];
(*zoctupper1=Module[{xxt,yyt},yyt=1/Xoct[[2,1]]//Simplify;
xxt=yyt Xoct[[1,1]]//Simplify;
{xxt,yyt}];
test=zoctupper1-zoctupper//Simplify;*)

{"{x,y}",zoctupper};
4 Module[{Iz,xxt,yyt},{xxt,yyt}=zoctupper;
Iz=(3/2) NIntegrate[(xxt^2+yyt^2+1)/yyt/.alphaoctsub,{t,0,tmaxoct}]];


theta1[x_,y_] := 2*(x^2-y^2+1)/(x^4+2x^2*(y^2+1)+(y^2-1)^2);
theta2 [x_,y_]:= 4*x*y/(x^4+4*(y^2+1)+(y^2-1)^2);
theta1[zoctupper[1],zoctupper[2]] === octsol0[1];
theta2[zoctupper[1],zoctupper[2]] === octsol0[2];
u = (x^2+y^2-1)/(x^2+(y+1)^2);
v = -2*x/(x^2+(y+1)^2);
theta = ArcTan[v/u];
r = Sqrt[u^2+v^2];

alphaoctsub' = alphaoctsub[[All, 2]];
zoct1 = First[zoctupper]/.{alpha -> alphaoctsub'};
zoct2 = Last[zoctupper]/.{alpha -> alphaoctsub'};
octsolx = First[octsol0];
octsoly = Last[octsol0];
thetaoctx = D[theta, x]/.{x -> zoct1, y -> zoct2};
thetaocty = D[theta, y] /.{x -> zoct1, y-> zoct2};
roctx = D[r,x]/.{x -> zoct1, y -> zoct2};
rocty = D[r,y]/.{x -> zoct1, y-> zoct2};

foctx = (octsolx - thetaoctx)/rx;
focty = (octsoly - thetaocty)/ry;
foctx === focty;

foctx0 = (octsol0[[1]] - thetaoctx);
focty0 = (octsol0[[1]] - thetaocty);
c45sub' = c45sub[[All, 2]];
c123sub' = c123sub[[All, 2]];
c4 = c45sub'[[1]]/.{lambda0 -> -1};
c5 = c45sub'[[2]]/.{lambda0 -> -1};
c1 = c123sub'[[1]]/.{lambda0 -> -1, C[4] -> c4, C[5] -> c5};
c2 = c123sub'[[2]]/.{lambda0 -> -1, C[4] -> c4, C[5] -> c5};
c3 = c123sub'[[3]]/.{lambda0 -> -1, C[4] -> c4, C[5] -> c5};

foctx0' = foctx0/.{C[1] -> c1, C[2] -> c2, C[4] -> c4, lambda0 ->-1, t -> 0};
focty0' = focty0/.{C[1] -> c1, C[2] -> c2, C[4] -> c4, lambda0 ->-1, t -> 0};

N[foctx0', 4]
N[focty0', 4]

F(x_,y_) := foctx0'*x + focty0'*y;
