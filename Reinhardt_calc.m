(* ::Package:: *)

(* preliminary defintions *)

L = {{c11,c12},{c21, -c11}};
J = {{0,-1},{1,0}};
X = {{x/y,1/y},{-x^2/y-y,-x/y}};
HLie := Tr[L.X] - 3/2*Tr[J.X];
a = (u1-u2)/sqrt[3];
b = u0/3 - 2u1/3-2u2/3;
c= u0;
f1=(y*(b+2a*x-c*x^2+c*y^2))/(b+2*a*x-c*x^2-c*y^2);
f2 = (2*(a-c*x)*y^2)/(b+2*a*x-c*x^2-c*y^2);
Hh = v1*f1 + v2*f2;
H = HLie + Hh;

(* Compute partials of the Hamiltonian for arbitrary controls *)

D[H,x] ;
D[H,y];
D[H,v1];
D[H,v2];
D[H, c11]
D[H, c12]
D[H, c21]

(* partials at bang-bang controls *)

D[H,x]/.{u0 -> 0, u1 -> 0, u2 -> 1}
D[H,x]/.{u0 -> 0, u1 -> 1, u2 -> 0}
D[H,x]/.{u0 -> 1, u1 -> 0, u2 -> 0}
D[H,y]/.{u0 -> 0, u1 -> 0, u2 -> 1}
D[H,y]/.{u0 -> 0, u1 -> 1, u2 -> 0}
D[H,y]/.{u0 -> 1, u1 -> 0, u2 -> 0}
D[H,v1]/.{u0 -> 0, u1 -> 0, u2 -> 1}
D[H,v1]/.{u0 -> 0, u1 -> 1, u2 -> 0}
D[H,v1]/.{u0 -> 1, u1 -> 0, u2 -> 0}
D[H,v2]/.{u0 -> 0, u1 -> 0, u2 -> 1}
D[H,v2]/.{u0 -> 0, u1 -> 1, u2 -> 0}
D[H,v2]/.{u0 -> 1, u1 -> 0, u2 -> 0}

(* Hessian Computations *)

D[H, {{x,y}, 2}]
D[H, {{v1,v2}, 2}]
Hess := D[H, {{x,y,v1,v2,c11,c12,c21}, 2}];
Hess/.{u0 ->1 , u1 -> 0, u2 -> 0}
Hess/.{u0 ->0, u1 -> 1, u2 -> 0}
Hess/.{u0 ->0 , u1 -> 0, u2 -> 1}
