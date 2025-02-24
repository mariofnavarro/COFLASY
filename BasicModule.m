(* ::Package:: *)

(* ::Section:: *)
(*Basic module of COFLASY-M*)


(* ::Text:: *)
(*by Valerie Domcke (valerie.domcke@cern.ch), Miguel Escudero (miguel.escudero@cern.ch), Mario Fernandez Navarro (Mario.FernandezNavarro@glasgow.ac.uk), Stefan Sandner (stefan.sandner@lanl.gov)*)
(**)
(*If you use this code for scientific publications, please cite the paper :*)
(*  *)
(*"Lepton Flavour Asymmetries: from the early Universe to BBN",  Valerie Domcke, Miguel Escudero, Mario Fernandez Navarro and Stefan Sandner*)
(*[arXiv:2502.14960] (https://arxiv.org/abs/2502.14960), [INSPIRE] (https://inspirehep.net/literature/2893306) .*)
(*  *)


(* ::Text:: *)
(*This module contains all the definitions, constants and functions to solve the system of ODEs.*)
(**)
(*All equation numbering in this module refers to the technical notes in pdf accompanying this code, unless otherwise specified.*)


(* ::Subsubsection::Closed:: *)
(*Preliminaries*)


(*These are all constants to be used later*)

MeV=1;
GeV=10^3 MeV;
eV=10^-6 MeV;
MPr=2.435*10^18 GeV; (*reduced Planck mass*)
GF=1.166 10^-11 MeV;
MW=80.377 GeV;
SW2=0.223;(*squared sine of the weak mixing angle*)
me = 0.511MeV;
m\[Mu]=105.7MeV;

(*Mixing parameters used in the paper*)
\[CapitalDelta]msq21=7.53*10^-5 eV^2;
\[CapitalDelta]msq31=2.534*10^-3 eV^2; (*For inverted ordering, add minus sign*)
\[Theta]12=0.587;
\[Theta]23=0.831;
\[Theta]13=0.148;

\[Delta]CP=0.0;(*This may be modified by the user but typically has small impact on the final results, as described in the paper the largest difference (although still small) is with the variation of \[Theta]13*)

FAC\[CapitalDelta]n=1; (*Controls globally the overall normalization of asymmetries \[CapitalDelta]n=FAC\[CapitalDelta]n(n-nbar)/T^3, FAC\[CapitalDelta]n=1 matches the normalisation in our paper while e.g. FAC\[CapitalDelta]n=100*6 is used in some figures of [2110.11889]*)



(*These are useful functions to be used later*)

commutator[A_,B_]:=A . B-B . A;

unitmat=({
 {1, 0, 0},
 {0, 1, 0},
 {0, 0, 1}
});
gellMannMatrices={({
 {0, 1, 0},
 {1, 0, 0},
 {0, 0, 0}
}),({
 {0, -I, 0},
 {I, 0, 0},
 {0, 0, 0}
}),({
 {1, 0, 0},
 {0, -1, 0},
 {0, 0, 0}
}),({
 {0, 0, 1},
 {0, 0, 0},
 {1, 0, 0}
}),({
 {0, 0, -I},
 {0, 0, 0},
 {I, 0, 0}
}),({
 {0, 0, 0},
 {0, 0, 1},
 {0, 1, 0}
}),({
 {0, 0, 0},
 {0, 0, -I},
 {0, I, 0}
}),1/Sqrt[3] ({
 {1, 0, 0},
 {0, 1, 0},
 {0, 0, -2}
})}; 

(*Structure constants of SU(3), fijk = 1/(4I)Tr(Subscript[\[Lambda], i][Subscript[\[Lambda], j],Subscript[\[Lambda], k]]), with Subscript[\[Lambda], i] being Gell-Mann matrices, i,j,k=1,...,8*)
fijk[i_,j_,k_]:=1/(4I) Tr[gellMannMatrices[[i]] . commutator[gellMannMatrices[[j]],gellMannMatrices[[k]]]];(*Eq. (1)*)


(*This is  Eq. (32) of the paper, which can be used to decompose a general 3x3 complex matrix into its components in the basis (Identity, Gell-Mann)*)

(*gellMannComponents[X_]:=FullSimplify[{1/3Tr[X],Table[1/2Tr[X.GellMannMatrices[[i]]],{i,1,8}]}];*)

(*However, we can also use the following function for the decomposition*)

gellMannExpansion[a0_,a1_,a2_,a3_,a4_,a5_,a6_,a7_,a8_]:=a0*unitmat+a1*gellMannMatrices[[1]]+a2*gellMannMatrices[[2]]+a3*gellMannMatrices[[3]]+a4*gellMannMatrices[[4]]+a5*gellMannMatrices[[5]]+a6*gellMannMatrices[[6]]+a7*gellMannMatrices[[7]]+a8*gellMannMatrices[[8]];

gellMannComponents[X_]:=FullSimplify[Solve[gellMannExpansion[a0,a1,a2,a3,a4,a5,a6,a7,a8]==X,{a0,a1,a2,a3,a4,a5,a6,a7,a8}]]//First; 

(*We find the function above to deliver the same analytical results for the decomposition of all matrices except for the a0 component of the FD collision terms,
 where we find the two functions to deliver potentially different analytical functions which nevertheless agree numerically extremely well,
 with relative differences < 10^-12% in all cases when evaluating both functions with 200 random density matrices.
 However, the function above seems to deliver a more compact decomposition of the collision terms that is better handled by Mathematica,
 leading to a factor 2 faster evolution times in NDSolve with respect to using Eq. (32) of the paper for the decomposition*)


(* ::Subsubsection::Closed:: *)
(*Equations*)


Clear[r0,r1,r2,r3,r4,r5,r6,r7,r8,rbar0,rbar1,rbar2,rbar3,rbar4,rbar5,rbar6,rbar7,rbar8];

rvec={r1,r2,r3,r4,r5,r6,r7,r8};
rbarvec={rbar1,rbar2,rbar3,rbar4,rbar5,rbar6,rbar7,rbar8};

(*Density matrices decomposition, Eq. (15)*)

rmat=r0*unitmat+rvec . gellMannMatrices; 
rbarmat=rbar0*unitmat+rbarvec . gellMannMatrices;

(*stictly speaking the full density matrices are \[Rho]=fFD(\[Xi]=0)*r and Overscript[\[Rho], _]=fFD(\[Xi]=0)*Overscript[r, _]) given our ansatz in Eq. (2)*)



(*This is a reference temperature used to build the dimensionless variable x=(Tref/T). Tref=me is a common convention, but other values such as Tref=10 MeV may be used.*)
Tref=me; 

gs=2+4 7/8+7/8 Tr[rmat+rbarmat]; (*Effective number of degrees of freedom, Eq. (5)*)
hubble[x_]:=Sqrt[gs] \[Pi]/(Sqrt[90]MPr) (Tref/x)^2; (*Hubble parameter, Eq. (4)*)


(*Note that the following terms correspond to the various contributions of Eq. (3)*)

(*H0: vaccuum oscillations*)
uUnitary=\!\(\*
TagBox[
RowBox[{"(", GridBox[{
{"1", "0", "0"},
{"0", 
RowBox[{"Cos", "[", "\[Theta]23", "]"}], 
RowBox[{"Sin", "[", "\[Theta]23", "]"}]},
{"0", 
RowBox[{"-", 
RowBox[{"Sin", "[", "\[Theta]23", "]"}]}], 
RowBox[{"Cos", "[", "\[Theta]23", "]"}]}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\) . \!\(\*
TagBox[
RowBox[{"(", GridBox[{
{
RowBox[{"Cos", "[", "\[Theta]13", "]"}], "0", 
RowBox[{
RowBox[{"Sin", "[", "\[Theta]13", "]"}], 
RowBox[{"Exp", "[", 
RowBox[{
RowBox[{"-", "I"}], " ", "\[Delta]CP"}], "]"}]}]},
{"0", "1", "0"},
{
RowBox[{
RowBox[{"-", 
RowBox[{"Sin", "[", "\[Theta]13", "]"}]}], " ", 
RowBox[{"Exp", "[", 
RowBox[{"I", " ", "\[Delta]CP"}], "]"}]}], "0", 
RowBox[{"Cos", "[", "\[Theta]13", "]"}]}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\) . \!\(\*
TagBox[
RowBox[{"(", GridBox[{
{
RowBox[{"Cos", "[", "\[Theta]12", "]"}], 
RowBox[{"Sin", "[", "\[Theta]12", "]"}], "0"},
{
RowBox[{"-", 
RowBox[{"Sin", "[", "\[Theta]12", "]"}]}], 
RowBox[{"Cos", "[", "\[Theta]12", "]"}], "0"},
{"0", "0", "1"}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\); (*PMNS matrix*)

m2diag={{0,0,0},{0,\[CapitalDelta]msq21,0},{0,0,\[CapitalDelta]msq31}};

h0mat[x_]:=(x/Tref) \[Pi]^2/(36 Zeta[3])*FullSimplify[uUnitary . m2diag . ConjugateTranspose[uUnitary]]; (*Eq. (6)*)

(*Note that we project into the basis of Gell-Mann matrices and write the vectors in units of Hx, following the convention of the paper*)
h0vec[x_]:=1/(hubble[x]x) {a1,a2,a3,a4,a5,a6,a7,a8}/.gellMannComponents[h0mat[x]]; 

(*Vc: matter potential*)

(*Electron integrals computed via expansion for small Subscript[m, e]/T*)
epluspe[x_]:=(7\[Pi]^2)/45-1/6 me^2/T^2//.{T->Tref/x};(*Eq. (9)*)

(*Muon integrals computed in the Maxwell-Boltzmann approximation*)
eplusp\[Mu][x_]:=2/\[Pi]^2 (m\[Mu]/T)^3 BesselK[3,m\[Mu]/T]//.{T->Tref/x};(*Eq. (10)*)

eplusp\[Tau]=0.0; (*taus are assumed to be completely absent from the plasma*)

eplusPMat[x_]:=\!\(\*
TagBox[
RowBox[{"(", GridBox[{
{
RowBox[{"epluspe", "[", "x", "]"}], "0", "0"},
{"0", 
RowBox[{"eplusp\[Mu]", "[", "x", "]"}], "0"},
{"0", "0", "eplusp\[Tau]"}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\);

vcmat[x_]:=-2 Sqrt[2] (7\[Pi]^4)/(180 Zeta[3]) GF/MW^2 (Tref/x)^5 eplusPMat[x]; (*Eq. (7)*)


vcvec[x_]:=1/(hubble[x]x) {a1,a2,a3,a4,a5,a6,a7,a8}/.gellMannComponents[vcmat[x]]; 

(*Vs term*)

vs[x_]:=1/(hubble[x]x) 3/(2Sqrt[2]) GF Zeta[3]/\[Pi]^2 (Tref/x)^3; (*Eq. (8)*)

(*FD collisions implemented as in Eqs. (11-13)*)

eFD=0.049;(* This quantity is the expansion rarameter that controls the collsion term for FD statistics as compared for MB. It should be 0.049 for FD and 0 for MB. *)
Fac=0.995;(* This term multiplies the whole set of collision terms. It should be 0.995 for FD statistics and 1 for MB. *)

gL=SW2-1/2;gR=SW2;
GLmat=({
 {gL+1, 0, 0},
 {0, gL, 0},
 {0, 0, gL}
});GRmat=({
 {gR, 0, 0},
 {0, gR, 0},
 {0, 0, gR}
});

(*Eq. (11)*)
\[ScriptCapitalI]gL\[Nu]escatt=Fac 32/\[Pi]^3 GF^2 T^5 (1-eFD)(GLmat . rmat . GLmat . (unitmat-eFD rmat)+(unitmat-eFD rmat) . GLmat . rmat . GLmat-(rmat . GLmat . (unitmat-eFD rmat) . GLmat+GLmat . (unitmat-eFD rmat) . GLmat . rmat));
(*Eq. (12)*)
\[ScriptCapitalI]gR\[Nu]esann=Fac 8/\[Pi]^3 GF^2 T^5 ((GRmat . (unitmat-eFD rbarmat) . GRmat . (unitmat-eFD rmat)+(unitmat-eFD rmat) . GRmat . (unitmat-eFD rbarmat) . GRmat)-(1-eFD)^2 (rmat . GRmat . rbarmat . GRmat+GRmat . rbarmat . GRmat . rmat));
\[ScriptCapitalI]gL\[Nu]esann=Fac 8/\[Pi]^3 GF^2 T^5 ((GLmat . (unitmat-eFD rbarmat) . GLmat . (unitmat-eFD rmat)+(unitmat-eFD rmat) . GLmat . (unitmat-eFD rbarmat) . GLmat)-(1-eFD)^2 (rmat . GLmat . rbarmat . GLmat+GLmat . rbarmat . GLmat . rmat));
(*Eq. (13)*)
\[ScriptCapitalI]\[Nu]\[Nu]ann=Fac 1/4 8/\[Pi]^3 GF^2 T^5 (2Tr[rmat . rbarmat]*(unitmat-eFD(rmat+rbarmat))- (rmat . rbarmat+rbarmat . rmat)(Tr[unitmat]-eFD Tr[rmat+rbarmat]));

(*The same but for nubar*)
\[ScriptCapitalI]gL\[Nu]escattbar=Fac 32/\[Pi]^3 GF^2 T^5 (1-eFD)(GLmat . rbarmat . GLmat . (unitmat-eFD rbarmat)+(unitmat-eFD rbarmat) . GLmat . rbarmat . GLmat-(rbarmat . GLmat . (unitmat-eFD rbarmat) . GLmat+GLmat . (unitmat-eFD rbarmat) . GLmat . rbarmat));
\[ScriptCapitalI]gR\[Nu]esannbar=Fac 8/\[Pi]^3 GF^2 T^5 ((GRmat . (unitmat-eFD rmat) . GRmat . (unitmat-eFD rbarmat)+(unitmat-eFD rbarmat) . GRmat . (unitmat-eFD rmat) . GRmat)-(1-eFD)^2 (rbarmat . GRmat . rmat . GRmat+GRmat . rmat . GRmat . rbarmat));
\[ScriptCapitalI]gL\[Nu]esannbar=Fac 8/\[Pi]^3 GF^2 T^5 ((GLmat . (unitmat-eFD rmat) . GLmat . (unitmat-eFD rbarmat)+(unitmat-eFD rbarmat) . GLmat . (unitmat-eFD rmat) . GLmat)-(1-eFD)^2 (rbarmat . GLmat . rmat . GLmat+GLmat . rmat . GLmat . rbarmat));

\[ScriptCapitalI]\[Nu]\[Nu]annbar=Fac 1/4 8/\[Pi]^3 GF^2 T^5 (2Tr[rmat . rbarmat]*(unitmat-eFD(rmat+rbarmat))- (rmat . rbarmat+rbarmat . rmat)(Tr[unitmat]-eFD Tr[rmat+rbarmat]));


(* Total sum of the collision terms *)

\[ScriptCapitalI]totFDmat[x_]:=(\[ScriptCapitalI]gL\[Nu]escatt+\[ScriptCapitalI]gR\[Nu]esann+\[ScriptCapitalI]gL\[Nu]esann+\[ScriptCapitalI]\[Nu]\[Nu]ann)//.{T->Tref/x};

\[ScriptCapitalI]totFDbarmat[x_]:=(\[ScriptCapitalI]gL\[Nu]escattbar+\[ScriptCapitalI]gR\[Nu]esannbar+\[ScriptCapitalI]gL\[Nu]esannbar+\[ScriptCapitalI]\[Nu]\[Nu]annbar)//.{T->Tref/x};

ctotFDvec[x_]:=1/(hubble[x]x) {a1,a2,a3,a4,a5,a6,a7,a8}/.gellMannComponents[\[ScriptCapitalI]totFDmat[x]];
cbartotFDvec[x_]:=1/(hubble[x]x) {a1,a2,a3,a4,a5,a6,a7,a8}/.gellMannComponents[\[ScriptCapitalI]totFDbarmat[x]];

ctotFD0[x_]:=1/(hubble[x]x) a0/.gellMannComponents[\[ScriptCapitalI]totFDmat[x]];
cbartotFD0[x_]:=1/(hubble[x]x) a0/.gellMannComponents[\[ScriptCapitalI]totFDbarmat[x]];

(*Collision damping terms*)

\[ScriptCapitalI]DampingMat[x_]:=-(Tref/x)^5 (16 Fac GF^2 (1-eFD)(3+2 SW2^2))/\[Pi]^3 ({
 {0, r12, r13},
 {r21, 0, 0},
 {r31, 0, 0}
});  (*Eq. (14)*)
\[ScriptCapitalI]barDampingMat[x_]:=-(Tref/x)^5 (16  Fac (GF^2) (1-eFD)(3+2 SW2^2) )/\[Pi]^3 ({
 {0, rbar12, rbar13},
 {rbar21, 0, 0},
 {rbar31, 0, 0}
}); (*The same but for rbar*)

cDvec[x_]:=1/(hubble[x]x) {a1,a2,a3,a4,a5,a6,a7,a8}/.gellMannComponents[\[ScriptCapitalI]DampingMat[x]]/.{r12->rmat[[1,2]],r13->rmat[[1,3]],r21->rmat[[2,1]],r31->rmat[[3,1]]} ;
cbarDvec [x_]:=1/(hubble[x]x) {a1,a2,a3,a4,a5,a6,a7,a8}/.gellMannComponents[\[ScriptCapitalI]barDampingMat[x]]/.{rbar12->rbarmat[[1,2]],rbar13->rbarmat[[1,3]],rbar21->rbarmat[[2,1]],rbar31->rbarmat[[3,1]]} ;


(*Set of QKEs in vector form given by Eq. (16), allowing for both damping or FD collisions via input*)
(*Note that r0 and rbar0 are constant only when we consider damping collision terms, the system then contains 8*2 ODEs, otherwise it contains 9*2*)

drvec[x_,VsControl_,CollDamping_,CollFD_]:=FullSimplify[Insert[2Table[Sum[fijk[i,j,k]rvec[[k]](h0vec[x]+vcvec[x]-VsControl vs[x]rbarvec)[[j]],{j,1,8},{k,1,8}],{i,1,8}]+CollDamping cDvec[x]+ CollFD ctotFDvec[x],CollFD ctotFD0[x],1]];
drbarvec[x_,VsControl_,CollDamping_,CollFD_]:=FullSimplify[Insert[-2Table[Sum[fijk[i,j,k]rbarvec[[k]](h0vec[x]+vcvec[x]-VsControl vs[x]rvec)[[j]],{j,1,8},{k,1,8}],{i,1,8}]+CollDamping cbarDvec[x]+ CollFD cbartotFDvec[x],CollFD cbartotFD0[x],1]];

Clear[rhsvecFD,rhsvecbarFD,rhsvecD,rhsvecbarD];

rtorx={r0->r0[x],r1->r1[x],r2->r2[x],r3->r3[x],r4->r4[x],r5->r5[x],r6->r6[x],r7->r7[x],r8->r8[x],rbar0->rbar0[x],rbar1->rbar1[x],rbar2->rbar2[x],rbar3->rbar3[x],rbar4->rbar4[x],rbar5->rbar5[x],rbar6->rbar6[x],rbar7->rbar7[x],rbar8->rbar8[x]};

(*We store equations both in the FD and damping cases, to be used later*)
(*The function switch[x] is later used to switch between full system -> adiabatic during (or before) the evolution*)

rhsvecFD=drvec[x,switchVs[x],0.0,1.0]/.rtorx;
rhsvecbarFD=drbarvec[x,switchVs[x],0.0,1.0]/.rtorx;

rhsvecD=drvec[x,switchVs[x],1.0,0.0]/.rtorx;
rhsvecbarD=drbarvec[x,switchVs[x],1.0,0.0]/.rtorx;


(* ::Subsubsection::Closed:: *)
(*Initial conditions*)


(*Useful equations to connect (A,\[Phi]) to \[Xi]t=\[Xi]+\[Xi]^3/\[Pi]^2*)
(*Eq. (20)*)
\[Xi]tile[A_,\[Phi]_]:=A Sin[\[Theta]]Cos[\[Phi]]/.{\[Theta]->ArcTan[-Cos[\[Phi]]-Sin[\[Phi]],1]};
\[Xi]til\[Mu][A_,\[Phi]_]:=A Sin[\[Theta]]Sin[\[Phi]]/.{\[Theta]->ArcTan[-Cos[\[Phi]]-Sin[\[Phi]],1]};
\[Xi]til\[Tau][A_,\[Phi]_]:=A Cos[\[Theta]]/.{\[Theta]->ArcTan[-Cos[\[Phi]]-Sin[\[Phi]],1]};

(*Solve[\[Xi]t==\[Xi]+\[Xi]^3/\[Pi]^2,\[Xi]][[1,1]]//ToRadicals//Simplify*)
(*Inverting Eq. (22)*)
\[Xi]tildeto\[Xi][\[Xi]t_]:=N[((\[Pi]/6)^(2/3) (2 3^(1/3) \[Pi]^(2/3)-2^(1/3) (-9 \[Xi]t+Sqrt[12 \[Pi]^2+81 \[Xi]t^2])^(2/3)))/(-9 \[Xi]t+Sqrt[12 \[Pi]^2+81 \[Xi]t^2])^(1/3)];




(*Initial conditions given our momentum average ansatz*)

riniFD[\[Xi]_]:=(-4 PolyLog[3,-E^\[Xi]])/(3Zeta[3]);(*Eq. (18)*)

(*Added to preserve detailed balance given our ansatz*)

\[Epsilon]db[\[Xi]_]:=1/(-2+4 eFD) (riniFD[-\[Xi]]-Sqrt[4+4 eFD^2 (1+(riniFD[-\[Xi]]-riniFD[\[Xi]])^2)-4 eFD (2+(riniFD[-\[Xi]]-riniFD[\[Xi]])^2)+(riniFD[-\[Xi]]-riniFD[\[Xi]])^2]+riniFD[\[Xi]]-2 eFD (-1+riniFD[-\[Xi]]+riniFD[\[Xi]])); (*Eq. (19)*)

(*Eq. (17)*)
rinifMat[\[Xi]e_,\[Xi]\[Mu]_,\[Xi]\[Tau]_]:=({
 {riniFD[\[Xi]e]+\[Epsilon]db[\[Xi]e], 0, 0},
 {0, riniFD[\[Xi]\[Mu]]+\[Epsilon]db[\[Xi]\[Mu]], 0},
 {0, 0, riniFD[\[Xi]\[Tau]]+\[Epsilon]db[\[Xi]\[Tau]]}
});
rbarinifMat[\[Xi]e_,\[Xi]\[Mu]_,\[Xi]\[Tau]_]:=({
 {riniFD[-\[Xi]e]+\[Epsilon]db[\[Xi]e], 0, 0},
 {0, riniFD[-\[Xi]\[Mu]]+\[Epsilon]db[\[Xi]\[Mu]], 0},
 {0, 0, riniFD[-\[Xi]\[Tau]]+\[Epsilon]db[\[Xi]\[Tau]]}
});

rinif[\[Xi]e_,\[Xi]\[Mu]_,\[Xi]\[Tau]_]:={a1,a2,a3,a4,a5,a6,a7,a8}/.gellMannComponents[rinifMat[\[Xi]e,\[Xi]\[Mu],\[Xi]\[Tau]]];
rbarinif[\[Xi]e_,\[Xi]\[Mu]_,\[Xi]\[Tau]_]:={a1,a2,a3,a4,a5,a6,a7,a8}/.gellMannComponents[rbarinifMat[\[Xi]e,\[Xi]\[Mu],\[Xi]\[Tau]]];
r0inif[\[Xi]e_,\[Xi]\[Mu]_,\[Xi]\[Tau]_]:=a0/.gellMannComponents[rinifMat[\[Xi]e,\[Xi]\[Mu],\[Xi]\[Tau]]];
rbar0inif[\[Xi]e_,\[Xi]\[Mu]_,\[Xi]\[Tau]_]:=a0/.gellMannComponents[rbarinifMat[\[Xi]e,\[Xi]\[Mu],\[Xi]\[Tau]]];


(* ::Subsubsection::Closed:: *)
(*Functions to run, plot and perform oscillation average*)


(* ::Text:: *)
(*This section contains two functions to solve the system of ODEs : solveFD : Solves the full system of ODEs with FD collisions, solveDamping : Solves the full system of ODEs with damping collisions, Both functions may consider the adiabatic approximation or not depending on user input . A timeout set by the user also switches non-adiabiatic -> adiabatic during the evolution*)
(**)
(*These functions receive : *)
(*(\[Xi]e, \[Xi]\[Mu], \[Xi]\[Tau]) : initial values of the initial reduced chemical potentials for each flavour, *)
(*Filename : string containing the name for the . dat file were the output is saved, *)
(*Tini : the initial temperature in MeV, *)
(*Tave : the temperature at which the solver finishes and oscillation average is taken in MeV, *)
(*Tfinal : the results at Tave is extended asymptotically until this final temperature, *)
(*Accval : accuracy setting for NDSolve, *)
(*Precval : precision setting for NDSolve, *)
(*SolverMethod : method used by NDSolve to solve the ODEs, e . g . "BDF" or "StiffnessSwitching", *)
(*Adiabatic : True in order to use the adiabatic approximation from the beginning, False to solve the full system, *)
(*TimeoutAdiabatic : When the full system is being solved, this is the execution time in seconds after which the solver switches to the adiabatic approximation*)
(**)
(*Note that by definition one should respect Tini < Tave <= Tfinal*)
(**)
(*The output in all cases is a Mathematica plot and a saved . dat file with the results*)


solveFD[\[Xi]e_,\[Xi]\[Mu]_,\[Xi]\[Tau]_,Filename_:"output",Tini_:20MeV,Tave_:1.5MeV,Tfinal_:1MeV,Accval_:8,Precval_:8,SolverMethod_:"BDF",Adiabatic_:False,TimeoutAdiabatic_:300]:=(
(*Checking that user input of T values is consistent, otherwise default values are used*)
tini=AbsoluteTime[];
If[Tini>Tave>=Tfinal,Tinival=Tini;Taveval=Tave;Tfinalval=Tfinal,Tinival=20MeV;Taveval=1.5MeV;Tfinalval=1MeV;Print["Tini, Tave and Tfinal seem inconsistent, using instead Tini=20 MeV, Tave=1.5 MeV, Tfinal=1 MeV"]];

(* initial conditions *)
Clear[r0,rbar0,r1,r2,r3,r4,r5,r6,r7,r8,rbar1,rbar2,rbar3,rbar4,rbar5,rbar6,rbar7,rbar8];
r0ini=r0inif[\[Xi]e,\[Xi]\[Mu],\[Xi]\[Tau]];
rbar0ini=rbar0inif[\[Xi]e,\[Xi]\[Mu],\[Xi]\[Tau]];
rini=rinif[\[Xi]e,\[Xi]\[Mu],\[Xi]\[Tau]];
rbarini=rbarinif[\[Xi]e,\[Xi]\[Mu],\[Xi]\[Tau]];

xini=Tref/Tinival;
xave=Tref/Taveval;

(*Solving the system*)
ControlAdiabatic=Adiabatic;
If[ControlAdiabatic==True,Print["Starting evolution of the system with \!\(\*SubscriptBox[\(V\), \(s\)]\)=0"];VsControl=0.0,Print["Starting evolution of the full system"];VsControl=1.0];
tin=AbsoluteTime[];
Monitor[sol=First[NDSolve[{ 
{r0'[x],r1'[x],r2'[x],r3'[x],r4'[x],r5'[x],r6'[x],r7'[x],r8'[x]}==rhsvecFD,{rbar0'[x],rbar1'[x],rbar2'[x],rbar3'[x],rbar4'[x],rbar5'[x],rbar6'[x],rbar7'[x],rbar8'[x]}==rhsvecbarFD,
r0[xini]==r0ini,r1[xini]==rini[[1]],r2[xini]==rini[[2]],r3[xini]==rini[[3]],r4[xini]==rini[[4]],r5[xini]==rini[[5]],r6[xini]==rini[[6]],r7[xini]==rini[[7]],r8[xini]==rini[[8]],rbar0[xini]==rbar0ini,rbar1[xini]==rbarini[[1]],rbar2[xini]==rbarini[[2]],rbar3[xini]==rbarini[[3]],rbar4[xini]==rbarini[[4]],rbar5[xini]==rbarini[[5]],rbar6[xini]==rbarini[[6]],rbar7[xini]==rbarini[[7]],rbar8[xini]==rbarini[[8]], switchVs[xini]==VsControl,WhenEvent[(AbsoluteTime[]-tin>TimeoutAdiabatic)&&(ControlAdiabatic==False),{Print["Timeout of "<>ToString[TimeoutAdiabatic]<>"s reached at T= "<>ToString[Tref/x]<>" MeV \[Rule] Switching to adiabatic \!\(\*SubscriptBox[\(V\), \(s\)]\)=0"],switchVs[x]->0.0,ControlAdiabatic=True}]},{r0,r1,r2,r3,r4,r5,r6,r7,r8,rbar0,rbar1,rbar2,rbar3,rbar4,rbar5,rbar6,rbar7,rbar8,switchVs},{x,xini,xave},PrecisionGoal->Precval,AccuracyGoal->Accval,MaxSteps->10^8,MaxStepSize->0.1,Method->{SolverMethod},DiscreteVariables->{switchVs},EvaluationMonitor:>{(monitor=Row[{"T= ", CForm[Tref/x]," MeV"}]),(time=x)}]],{monitor,ProgressIndicator[time,{xini,xave}]}];

Print["System evolved until \!\(\*SubscriptBox[\(T\), \(ave\)]\)="<>ToString[Taveval]," MeV"];

(*Asymmetries are obtained via Eq. (22)*)
Clear[\[CapitalDelta]ne,\[CapitalDelta]n\[Mu],\[CapitalDelta]n\[Tau]];

\[CapitalDelta]ne[T_]:=FAC\[CapitalDelta]n 3/4 Zeta[3]/\[Pi]^2 Re[rmat[[1,1]]-rbarmat[[1,1]]]/.rtorx/.{x->Tref/T}/.sol;
\[CapitalDelta]n\[Mu][T_]:=FAC\[CapitalDelta]n 3/4 Zeta[3]/\[Pi]^2 Re[rmat[[2,2]]-rbarmat[[2,2]]]/.rtorx/.{x->Tref/T}/.sol;
\[CapitalDelta]n\[Tau][T_]:=FAC\[CapitalDelta]n 3/4 Zeta[3]/\[Pi]^2 Re[rmat[[3,3]]-rbarmat[[3,3]]]/.rtorx/.{x->Tref/T}/.sol;

(*Oscillation average*)

(*First I store the density matrix at Tave*)

rTave=rmat/.rtorx/.sol/.x->Tref/Taveval;
rbarTave=rbarmat/.rtorx/.sol/.x->Tref/Taveval;

(*Then I apply Eq. (21) to extract the asymptotic final asymetries*)

rFinal=uUnitary . DiagonalMatrix[Diagonal[ConjugateTranspose[uUnitary] . rTave . uUnitary]] . ConjugateTranspose[uUnitary];
rbarFinal=uUnitary . DiagonalMatrix[Diagonal[ConjugateTranspose[uUnitary] . rbarTave . uUnitary]] . ConjugateTranspose[uUnitary];

\[CapitalDelta]nefinal=FAC\[CapitalDelta]n 3/4 Zeta[3]/\[Pi]^2 Re[rFinal[[1,1]]-rbarFinal[[1,1]]];
\[CapitalDelta]n\[Mu]final=FAC\[CapitalDelta]n 3/4 Zeta[3]/\[Pi]^2 Re[rFinal[[2,2]]-rbarFinal[[2,2]]];
\[CapitalDelta]n\[Tau]final=FAC\[CapitalDelta]n 3/4 Zeta[3]/\[Pi]^2 Re[rFinal[[3,3]]-rbarFinal[[3,3]]];

(*Inverting Eq. (22)*)
\[Xi]\[Alpha]final={\[Xi]tildeto\[Xi][6\[CapitalDelta]nefinal],\[Xi]tildeto\[Xi][6\[CapitalDelta]n\[Mu]final],\[Xi]tildeto\[Xi][6\[CapitalDelta]n\[Tau]final]};
\[Xi]efinal=\[Xi]\[Alpha]final[[1]];

(*Relevant output*)

Print["Final \[Xi]e: ",\[Xi]efinal];

(* Plot *)

plotBeforeTave=LogLinearPlot[{\[CapitalDelta]n\[Tau][T],\[CapitalDelta]n\[Mu][T],\[CapitalDelta]ne[T]},{T,Tinival,Taveval},Frame->True,
FrameLabel->{Style["T [MeV]",FontSize->20],Style["\!\(\*FractionBox[SubscriptBox[\((n - \*OverscriptBox[\(n\), \(_\)])\), \(\[Alpha]\[Alpha]\)], SuperscriptBox[\(T\), \(3\)]]\)",FontSize->30]},PlotStyle->{{Darker@Blue,Thickness[0.005]},{Darker@Red,Thickness[0.005]},{Black,Thickness[0.005]}},ScalingFunctions->{"Reverse",Identity},BaseStyle->Times,FrameStyle->Directive[Black,Thick],RotateLabel->False,ImageSize->575,LabelStyle->{FontSize->20},PlotLegends->{SwatchLegend[{Black,Darker@Red,Darker@Blue},{"e","\[Mu]","\[Tau]"},LegendLayout->{"Column",1},LegendMarkers->{"Line"},LegendMarkerSize->25]}];

If[Taveval==Tfinalval,plotFinal=plotBeforeTave,plotAfterTave=LogLinearPlot[{\[CapitalDelta]n\[Tau]final,\[CapitalDelta]n\[Mu]final,\[CapitalDelta]nefinal},{T,Taveval,Tfinalval},Frame->True,
FrameLabel->{Style["T [MeV]",FontSize->20],Style["\!\(\*FractionBox[SubscriptBox[\((n - \*OverscriptBox[\(n\), \(_\)])\), \(\[Alpha]\[Alpha]\)], SuperscriptBox[\(T\), \(3\)]]\)",FontSize->30]},PlotStyle->{{Darker@Blue,Thickness[0.005]},{Darker@Red,Thickness[0.005]},{Black,Thickness[0.005]}},ScalingFunctions->{"Reverse",Identity},BaseStyle->Times,FrameStyle->Directive[Black,Thick],RotateLabel->False,ImageSize->575,LabelStyle->{FontSize->20}];plotFinal=Show[plotBeforeTave,plotAfterTave,PlotRange->Full];];

(*Export results to output .dat file*)

SetDirectory[NotebookDirectory[]];
Tvec=Reverse[Table[N[10^i],{i,Log10[Taveval],Log10[Tinival],(Log10[Tinival]-Log10[Taveval])/1000}]];
Matvec=Table[0,{i,1,Length[Tvec]+2},{j,1,4}];

Do[
Matvec[[i,1]]=Tvec[[i]];
Matvec[[i,2]]=\[CapitalDelta]ne[Tvec[[i]]];
Matvec[[i,3]]=\[CapitalDelta]n\[Mu][Tvec[[i]]];
Matvec[[i,4]]=\[CapitalDelta]n\[Tau][Tvec[[i]]];
,{i,1,Length[Tvec],1}];

Matvec[[Length[Tvec]+1,1]]=Taveval;
Matvec[[Length[Tvec]+1,2]]=\[CapitalDelta]nefinal;
Matvec[[Length[Tvec]+1,3]]=\[CapitalDelta]n\[Mu]final;
Matvec[[Length[Tvec]+1,4]]=\[CapitalDelta]n\[Tau]final;
Matvec[[Length[Tvec]+2,1]]=Tfinalval;
Matvec[[Length[Tvec]+2,2]]=\[CapitalDelta]nefinal;
Matvec[[Length[Tvec]+2,3]]=\[CapitalDelta]n\[Mu]final;
Matvec[[Length[Tvec]+2,4]]=\[CapitalDelta]n\[Tau]final;

Matvec=SetPrecision[Matvec[[All,All]],7];
Print["Thermodynamics are output to: output/"<> Filename<>".dat"];
Export["output/"<> Filename<>".dat",Matvec,"Table","FieldSeparators"->"      ","TableHeadings"-> {"# T[MeV]" ,"nnbare","nnbarmu","nnbartau"}];


tfinal=AbsoluteTime[];
Print["Total execution time: ",If[tfinal-tini<60,N[tfinal-tini,4]"s",N[(tfinal-tini)/60,4]"min"]];

Return[plotFinal];
);



solveDamping[\[Xi]e_,\[Xi]\[Mu]_,\[Xi]\[Tau]_,Filename_:"output",Tini_:20MeV,Tave_:1.5MeV,Tfinal_:1MeV,Accval_:8,Precval_:8,SolverMethod_:"BDF",Adiabatic_:False,TimeoutAdiabatic_:300]:=(
(*Checking that user input of T values is consistent, otherwise default values are used*)
tini=AbsoluteTime[];
If[Tini>Tave>=Tfinal,Tinival=Tini;Taveval=Tave;Tfinalval=Tfinal,Tinival=20MeV;Taveval=1.5MeV;Tfinalval=1MeV;Print["Tini, Tave and Tfinal seem inconsistent, using instead Tini=20 MeV, Tave=1.5 MeV, Tfinal=1 MeV"]];

(* initial conditions *)
Clear[r0,rbar0,r1,r2,r3,r4,r5,r6,r7,r8,rbar1,rbar2,rbar3,rbar4,rbar5,rbar6,rbar7,rbar8];
r0ini=r0inif[\[Xi]e,\[Xi]\[Mu],\[Xi]\[Tau]];
rbar0ini=rbar0inif[\[Xi]e,\[Xi]\[Mu],\[Xi]\[Tau]];
rini=rinif[\[Xi]e,\[Xi]\[Mu],\[Xi]\[Tau]];
rbarini=rbarinif[\[Xi]e,\[Xi]\[Mu],\[Xi]\[Tau]];

xini=Tref/Tinival;
xave=Tref/Taveval;

(*Solving the system*)
ControlAdiabatic=Adiabatic;
If[ControlAdiabatic==True,Print["Starting evolution of the system with \!\(\*SubscriptBox[\(V\), \(s\)]\)=0"];VsControl=0.0,Print["Starting evolution of the full system"];VsControl=1.0];
tin=AbsoluteTime[];
Monitor[sol=First[NDSolve[{ 
{r0'[x],r1'[x],r2'[x],r3'[x],r4'[x],r5'[x],r6'[x],r7'[x],r8'[x]}==rhsvecD,{rbar0'[x],rbar1'[x],rbar2'[x],rbar3'[x],rbar4'[x],rbar5'[x],rbar6'[x],rbar7'[x],rbar8'[x]}==rhsvecbarD,
r0[xini]==r0ini,r1[xini]==rini[[1]],r2[xini]==rini[[2]],r3[xini]==rini[[3]],r4[xini]==rini[[4]],r5[xini]==rini[[5]],r6[xini]==rini[[6]],r7[xini]==rini[[7]],r8[xini]==rini[[8]],rbar0[xini]==rbar0ini,rbar1[xini]==rbarini[[1]],rbar2[xini]==rbarini[[2]],rbar3[xini]==rbarini[[3]],rbar4[xini]==rbarini[[4]],rbar5[xini]==rbarini[[5]],rbar6[xini]==rbarini[[6]],rbar7[xini]==rbarini[[7]],rbar8[xini]==rbarini[[8]], switchVs[xini]==VsControl,WhenEvent[(AbsoluteTime[]-tin>TimeoutAdiabatic)&&(ControlAdiabatic==False),{Print["Timeout of "<>ToString[TimeoutAdiabatic]<>"s reached at T= "<>ToString[Tref/x]<>" MeV \[Rule] Switching to adiabatic \!\(\*SubscriptBox[\(V\), \(s\)]\)=0"],switchVs[x]->0.0,ControlAdiabatic=True}]},{r0,r1,r2,r3,r4,r5,r6,r7,r8,rbar0,rbar1,rbar2,rbar3,rbar4,rbar5,rbar6,rbar7,rbar8,switchVs},{x,xini,xave},PrecisionGoal->Precval,AccuracyGoal->Accval,MaxSteps->10^8,MaxStepSize->0.1,Method->{SolverMethod},DiscreteVariables->{switchVs},EvaluationMonitor:>{(monitor=Row[{"T= ", CForm[Tref/x]," MeV"}]),(time=x)}]],{monitor,ProgressIndicator[time,{xini,xave}]}];

Print["System evolved until \!\(\*SubscriptBox[\(T\), \(ave\)]\)="<>ToString[Taveval]," MeV"];

(*Asymmetries are obtained via Eq. (22)*)
Clear[\[CapitalDelta]ne,\[CapitalDelta]n\[Mu],\[CapitalDelta]n\[Tau]];

\[CapitalDelta]ne[T_]:=FAC\[CapitalDelta]n 3/4 Zeta[3]/\[Pi]^2 Re[rmat[[1,1]]-rbarmat[[1,1]]]/.rtorx/.{x->Tref/T}/.sol;
\[CapitalDelta]n\[Mu][T_]:=FAC\[CapitalDelta]n 3/4 Zeta[3]/\[Pi]^2 Re[rmat[[2,2]]-rbarmat[[2,2]]]/.rtorx/.{x->Tref/T}/.sol;
\[CapitalDelta]n\[Tau][T_]:=FAC\[CapitalDelta]n 3/4 Zeta[3]/\[Pi]^2 Re[rmat[[3,3]]-rbarmat[[3,3]]]/.rtorx/.{x->Tref/T}/.sol;

(*Oscillation average*)

(*First I store the density matrix at Tave*)

rTave=rmat/.rtorx/.sol/.x->Tref/Taveval;
rbarTave=rbarmat/.rtorx/.sol/.x->Tref/Taveval;

(*Then I apply Eq. (21) to extract the asymptotic final asymetries*)

rFinal=uUnitary . DiagonalMatrix[Diagonal[ConjugateTranspose[uUnitary] . rTave . uUnitary]] . ConjugateTranspose[uUnitary];
rbarFinal=uUnitary . DiagonalMatrix[Diagonal[ConjugateTranspose[uUnitary] . rbarTave . uUnitary]] . ConjugateTranspose[uUnitary];

\[CapitalDelta]nefinal=FAC\[CapitalDelta]n 3/4 Zeta[3]/\[Pi]^2 Re[rFinal[[1,1]]-rbarFinal[[1,1]]];
\[CapitalDelta]n\[Mu]final=FAC\[CapitalDelta]n 3/4 Zeta[3]/\[Pi]^2 Re[rFinal[[2,2]]-rbarFinal[[2,2]]];
\[CapitalDelta]n\[Tau]final=FAC\[CapitalDelta]n 3/4 Zeta[3]/\[Pi]^2 Re[rFinal[[3,3]]-rbarFinal[[3,3]]];

(*Inverting Eq. (22)*)
\[Xi]\[Alpha]final={\[Xi]tildeto\[Xi][6\[CapitalDelta]nefinal],\[Xi]tildeto\[Xi][6\[CapitalDelta]n\[Mu]final],\[Xi]tildeto\[Xi][6\[CapitalDelta]n\[Tau]final]};
\[Xi]efinal=\[Xi]\[Alpha]final[[1]];

(*Relevant output*)

Print["Final \[Xi]e: ",\[Xi]efinal];

(* Plot *)

plotBeforeTave=LogLinearPlot[{\[CapitalDelta]n\[Tau][T],\[CapitalDelta]n\[Mu][T],\[CapitalDelta]ne[T]},{T,Tinival,Taveval},Frame->True,
FrameLabel->{Style["T [MeV]",FontSize->20],Style["\!\(\*FractionBox[SubscriptBox[\((n - \*OverscriptBox[\(n\), \(_\)])\), \(\[Alpha]\[Alpha]\)], SuperscriptBox[\(T\), \(3\)]]\)",FontSize->30]},PlotStyle->{{Darker@Blue,Thickness[0.005]},{Darker@Red,Thickness[0.005]},{Black,Thickness[0.005]}},ScalingFunctions->{"Reverse",Identity},BaseStyle->Times,FrameStyle->Directive[Black,Thick],RotateLabel->False,ImageSize->575,LabelStyle->{FontSize->20},PlotLegends->{SwatchLegend[{Black,Darker@Red,Darker@Blue},{"e","\[Mu]","\[Tau]"},LegendLayout->{"Column",1},LegendMarkers->{"Line"},LegendMarkerSize->25]}];

If[Taveval==Tfinalval,plotFinal=plotBeforeTave,plotAfterTave=LogLinearPlot[{\[CapitalDelta]n\[Tau]final,\[CapitalDelta]n\[Mu]final,\[CapitalDelta]nefinal},{T,Taveval,Tfinalval},Frame->True,
FrameLabel->{Style["T [MeV]",FontSize->20],Style["\!\(\*FractionBox[SubscriptBox[\((n - \*OverscriptBox[\(n\), \(_\)])\), \(\[Alpha]\[Alpha]\)], SuperscriptBox[\(T\), \(3\)]]\)",FontSize->30]},PlotStyle->{{Darker@Blue,Thickness[0.005]},{Darker@Red,Thickness[0.005]},{Black,Thickness[0.005]}},ScalingFunctions->{"Reverse",Identity},BaseStyle->Times,FrameStyle->Directive[Black,Thick],RotateLabel->False,ImageSize->575,LabelStyle->{FontSize->20}];plotFinal=Show[plotBeforeTave,plotAfterTave,PlotRange->Full];];

(*Export results to output .dat file*)

SetDirectory[NotebookDirectory[]];
Tvec=Reverse[Table[N[10^i],{i,Log10[Taveval],Log10[Tinival],(Log10[Tinival]-Log10[Taveval])/1000}]];
Matvec=Table[0,{i,1,Length[Tvec]+2},{j,1,4}];

Do[
Matvec[[i,1]]=Tvec[[i]];
Matvec[[i,2]]=\[CapitalDelta]ne[Tvec[[i]]];
Matvec[[i,3]]=\[CapitalDelta]n\[Mu][Tvec[[i]]];
Matvec[[i,4]]=\[CapitalDelta]n\[Tau][Tvec[[i]]];
,{i,1,Length[Tvec],1}];

Matvec[[Length[Tvec]+1,1]]=Taveval;
Matvec[[Length[Tvec]+1,2]]=\[CapitalDelta]nefinal;
Matvec[[Length[Tvec]+1,3]]=\[CapitalDelta]n\[Mu]final;
Matvec[[Length[Tvec]+1,4]]=\[CapitalDelta]n\[Tau]final;
Matvec[[Length[Tvec]+2,1]]=Tfinalval;
Matvec[[Length[Tvec]+2,2]]=\[CapitalDelta]nefinal;
Matvec[[Length[Tvec]+2,3]]=\[CapitalDelta]n\[Mu]final;
Matvec[[Length[Tvec]+2,4]]=\[CapitalDelta]n\[Tau]final;

Matvec=SetPrecision[Matvec[[All,All]],7];
Print["Thermodynamics are output to: output/"<> Filename<>".dat"];
Export["output/"<> Filename<>".dat",Matvec,"Table","FieldSeparators"->"      ","TableHeadings"-> {"# T[MeV]" ,"nnbare","nnbarmu","nnbartau", "Neff_trace"}];


tfinal=AbsoluteTime[];
Print["Total execution time: ",If[tfinal-tini<60,N[tfinal-tini,4]"s",N[(tfinal-tini)/60,4]"min"]];

Return[plotFinal];
);
