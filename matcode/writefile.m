(* ::Package:: *)

Close[str];
Clear[str];


HessianH[f_, x_List?VectorQ] := D[f, {x, 2}];

Unprotect[Power];
Format[Power[a_,n_Integer?Positive], CForm] := Distribute[
    ConstantArray[Hold[a],n],
    Hold, List, HoldForm, Times
]
Protect[Power];

ClearAll[cform];
cform[(-1)]:= "-0.1e1";
cform[(1)]:= "0.1e1";

cform[Times[a_, b_]] := "("<>cform[a] <> ")" <> " * " <> "(" <> cform[b] <> ")"; 
cform[Rational[a_, b_]]:= cform[a]<>"/"<>cform[b]<>".";
cform[Power[a_, 3]]:= cform[a]<>"*"<>cform[a]<>"*"<>cform[a];
cform[Power[a_, 2]]:= cform[a]<>"*"<>cform[a];
cform[Power[a_, -1]]:= "0.1e1"<>"/"<>cform[a];
cform[Power[a_, -2]]:= "0.1e1"<>"/"<>cform[a]<>"*"<>cform[a]
cform[Power[a_, -3]]:= "0.1e1"<>"/"<>cform[a]<>"*"<>cform[a]<>"*"<>cform[a]
cform[Power[a_, 1/2]]:= ToString[Sqrt[cform[a]]];
cform[Power[a_, -1/2]]:= "0.1e1"<>"/"<>ToString[Sqrt[cform[a]]];

cform[Plus[a_, b_]] := cform[a] <> " + " <> cform[b];

cform[Complex[a_, b_]]:=  Which[a == 0 && b != 0, cform[b]<>"*"<>ToString[A], a != 0 && b==0, cform[a] ,a!=0 && b!=0, cform[a]<>"+"<>cform[b]<>"*"<>ToString[A]];
(*cform[a_?AtomQ] := ToString[a];*)
cform[h_[args__]] := cform[h] <> "(" <> StringJoin[Riffle[cform /@ {args}, ","]] <> ")";
(*cform[I] := ToString[II];*)

cform[Cg] := "C_g";
cform[Cgs] := "C_gs";

cform[CMassElectron]:= "C_MassElectron";
cform[CMassMu]:= "C_MassMu";
cform[CMassTau]:= "C_MassTau";
cform[CMassUp]:= "C_MassUp";
cform[CMassDown]:= "C_MassDown";
cform[CMassTop] := "C_MassTop";
cform[CMassBottom]:= "C_MassBottom";
cform[CMassStrange]:= "C_MassStrange";
cform[CMassCharm]:= "C_MassCharm";
cform[Conjugate[a_]]:= "conj"<>"("<>ToString[a]<>")";
cform[allOther_] := ToString[allOther, CForm];


(* Higgs fields configuration *)
\[CapitalPhi]1 = (1/Sqrt[2]){{\[Rho]1 + I \[Eta]1}, {\[Xi]1  + I \[Psi]1}};
\[CapitalPhi]2 = (1/Sqrt[2]){{\[Rho]2  + I \[Eta]2}, {\[Xi]2 + I \[Psi]2 }};
\[CapitalPhi]1T = (ConjugateTranspose[\[CapitalPhi]1]/.{Conjugate[\[Eta]1]->\[Eta]1, Conjugate[\[Rho]1]->\[Rho]1, Conjugate[\[Xi]1]->\[Xi]1, Conjugate[\[Psi]1]->\[Psi]1, Conjugate[\[Xi]1 + \[Omega]1]->\[Xi]1 + \[Omega]1})//Flatten//Simplify;
\[CapitalPhi]2T =  (ConjugateTranspose[\[CapitalPhi]2]/.{Conjugate[\[Eta]2]->\[Eta]2, Conjugate[\[Rho]2 ]->\[Rho]2 , Conjugate[\[Xi]2]->\[Xi]2, Conjugate[\[Psi]2]->\[Psi]2, Conjugate[\[Xi]2 + \[Omega]2]->\[Xi]2 + \[Omega]2, Conjugate[\[Psi]2 + \[Omega]CP]->\[Psi]2 + \[Omega]CP})//Flatten//Simplify;
\[CapitalPhi]S = \[Xi]S;

\[CapitalPhi]1c = Simplify[I*sigma2.Conjugate[\[CapitalPhi]1],\[CapitalPhi]real]; (* to generate up quark mass *)
\[CapitalPhi]real = {\[Rho]1 \[Element] Reals, \[Eta]1 \[Element] Reals, \[Xi]1 \[Element] Reals, \[Psi]1 \[Element] Reals, \[Rho]2 \[Element] Reals, \[Eta]2 \[Element] Reals, \[Xi]2 \[Element] Reals, \[Psi]2 \[Element] Reals};

(* Gauge fields configuration *)
Gaugebasis = {W1, W2, W3, B0};
Gaugereal = {W1 \[Element] Reals, W2 \[Element] Reals, W3 \[Element] Reals, B0 \[Element] Reals};

sigma0 = {{1,0},{0,1}};
sigma1 = {{0,1}, {1,0}};
sigma2 = {{0, -I},{I, 0}};
sigma3 = {{1,0},{0,-1}};
sigma = {sigma1, sigma2, sigma3, sigma0};
Dmu = -I (Cg/2) Sum[Gaugebasis[[i]] sigma[[i]], {i,1,Length@Gaugebasis-1}]- I (Cgs/2) Gaugebasis[[4]]sigma[[4]];
DmuT = I (Cg/2) ConjugateTranspose[Sum[Gaugebasis[[i]] sigma[[i]], {i,1,Length@Gaugebasis-1}]] + I (Cgs/2) ConjugateTranspose[Gaugebasis[[4]]sigma[[4]]];

(* Lepton fields configuration *)
PiLep = {{ye, 0,0},{0, ymu, 0},{0,0, ytau}};

NuL = {veL, vmuL, vtauL};
ER = {eR, muR, tauR};
EL = {eL, muL, tauL};
LepBase = {eL, eR, muL, muR, tauL, tauR, veL, vmuL, vtauL};

(* Quark field configuration *)
UL = {{uL}, {cL}, {tL}};
DL = {{dL}, {sL}, {bL}};
UR = {{uR}, {cR}, {tR}};
DR = {{dR}, {sR}, {bR}};

VCKM = {{V11, V12, V13}, {V21,V22, V23}, {V31, V32, V33}};
CKML = {V11->Vud, V12->Vus, V13->Vub, V21->Vcd, V22->Vcs, V23->Vcb, V31->Vtd, V32->Vts, V33->Vtb};

Dd = {{yd, 0, 0}, {0, ys, 0}, {0, 0, yb}};
DU = {{yu, 0, 0}, {0, yc, 0}, {0, 0, yt}};
QuarkBase = {uR, cR, tR, dR, sR, bR, uL, cL, tL, dL, sL, bL};


(* Tree level configurations of fields and vev-input *)
\[Phi] = {\[Rho]1, \[Eta]1, \[Rho]2, \[Eta]2, \[Psi]1, \[Psi]2, \[Xi]1,  \[Xi]2, \[Xi]S};
\[Phi]vev = {\[Rho]1->0, \[Eta]1->0, \[Rho]2->\[Omega]CB, \[Eta]2->0, \[Xi]1->\[Omega]1, \[Psi]1->0, \[Xi]2->\[Omega]2, \[Psi]2->\[Omega]CP, \[Xi]S->\[Omega]S};
(*rotate = {\[Omega]1\[Rule]v, \[Omega]2\[Rule]0, \[Omega]CB\[Rule]0, \[Omega]CP\[Rule]0, \[Omega]S\[Rule]0}*)
rotate = {\[Omega]1->vh, \[Omega]2->0, \[Omega]CB->0, \[Omega]CP->0, \[Omega]S->0};
rotate0 = {\[Omega]1->0, \[Omega]2->0, \[Omega]CB->0, \[Omega]CP->0, \[Omega]S->0};


(* Higgs potential *)
VTree  = ms11 \[CapitalPhi]1T.\[CapitalPhi]1 + ms22 \[CapitalPhi]2T.\[CapitalPhi]2  + (\[Lambda]1/2)(\[CapitalPhi]1T.\[CapitalPhi]1)^2 + (\[Lambda]2/2)(\[CapitalPhi]2T.\[CapitalPhi]2)^2  + \[Lambda]3 (\[CapitalPhi]1T.\[CapitalPhi]1) (\[CapitalPhi]2T.\[CapitalPhi]2) + \[Lambda]4  (\[CapitalPhi]1T.\[CapitalPhi]2) (\[CapitalPhi]2T.\[CapitalPhi]1) + (Areal + I Aimag)\[CapitalPhi]1T.\[CapitalPhi]2 \[CapitalPhi]S + (Areal - I Aimag) \[CapitalPhi]2T.\[CapitalPhi]1 \[CapitalPhi]S + (\[Lambda]5/2)( (\[CapitalPhi]1T.\[CapitalPhi]2)^2 +  (\[CapitalPhi]2T.\[CapitalPhi]1)^2) +mss \[CapitalPhi]S^2/2 + (\[Lambda]6/4)\[CapitalPhi]S^4 + (\[Lambda]7/2)(\[CapitalPhi]1T.\[CapitalPhi]1)\[CapitalPhi]S^2 +(\[Lambda]8/2)(\[CapitalPhi]2T.\[CapitalPhi]2)\[CapitalPhi]S^2;

(* counter-terms Higgs potential *)
Vtemp = \[Delta]T1 (\[Xi]1 + \[Omega]1) + \[Delta]T2 (\[Xi]2 + \[Omega]2) + \[Delta]TCP (\[Psi]2 + \[Omega]CP) + \[Delta]TCB (\[Rho]2 + \[Omega]CB) + \[Delta]TS (\[Xi]S + \[Omega]S);
VCT  = \[Delta]ms11 \[CapitalPhi]1T.\[CapitalPhi]1 + \[Delta]ms22 \[CapitalPhi]2T.\[CapitalPhi]2  + (\[Delta]\[Lambda]1/2)(\[CapitalPhi]1T.\[CapitalPhi]1)^2 + (\[Delta]\[Lambda]2/2)(\[CapitalPhi]2T.\[CapitalPhi]2)^2  + \[Delta]\[Lambda]3 (\[CapitalPhi]1T.\[CapitalPhi]1) (\[CapitalPhi]2T.\[CapitalPhi]2) + \[Delta]\[Lambda]4  (\[CapitalPhi]1T.\[CapitalPhi]2) (\[CapitalPhi]2T.\[CapitalPhi]1) + (\[Delta]Areal + I \[Delta]Aimag)\[CapitalPhi]1T.\[CapitalPhi]2 \[CapitalPhi]S + (\[Delta]Areal - I \[Delta]Aimag) \[CapitalPhi]2T.\[CapitalPhi]1 \[CapitalPhi]S + (\[Delta]\[Lambda]5/2)( (\[CapitalPhi]1T.\[CapitalPhi]2)^2 +  (\[CapitalPhi]2T.\[CapitalPhi]1)^2) +\[Delta]mss \[CapitalPhi]S^2/2 + (\[Delta]\[Lambda]6/4)\[CapitalPhi]S^4 + (\[Delta]\[Lambda]7/2)(\[CapitalPhi]1T.\[CapitalPhi]1)\[CapitalPhi]S^2 +(\[Delta]\[Lambda]8/2)(\[CapitalPhi]2T.\[CapitalPhi]2)\[CapitalPhi]S^2 + Vtemp//Simplify;

(* Gauge potential *)
Vgauge = FullSimplify[(\[CapitalPhi]1T.DmuT).(Dmu.\[CapitalPhi]1) + (\[CapitalPhi]2T.DmuT).(Dmu.\[CapitalPhi]2), Gaugereal][[1]];

(* Lepton potential *)
VFLep = (NuL.PiLep.ER)\[CapitalPhi]1[[1]] + (EL.PiLep.ER)\[CapitalPhi]1[[2]]//Simplify;

(* Quark potential *)
VF= Simplify[((Transpose[UL].VCKM.Dd.DR) \[CapitalPhi]1[[1]] + Transpose[DL].Dd.DR \[CapitalPhi]1[[2]] + Transpose[UL].DU.UR \[CapitalPhi]1c[[1]] + Transpose[DL].ConjugateTranspose[VCKM].DU.UR \[CapitalPhi]1c[[2]]//Simplify)[[1]][[1]], \[CapitalPhi]real];


(* relation among parameters expressing via independent parameters from minization of tree level potential *)
(* cannot automize now *)
minimize = ConstantArray[0, Length@\[Phi]];
Do[
minimize[[i]] = D[VTree[[1]], \[Phi][[i]]]/.\[Phi]vev/.rotate
,{i,1, Length@\[Phi]}];

eqs=DeleteCases[minimize, 0, Infinity]
linearrep = Solve[eqs[[1]]==0, {ms11}]

(*VTree = (VTree/.linearrep)[[1]] (* remove all dependent parameters from potential *);*)


(* Relation between Lepton Yukawa couplings and its masses *)
(* cannot automize now *)
MassLep = HessianH[VFLep[[1]], LepBase]/.\[Phi]vev/.rotate;
MLep = Eigenvalues[MassLep];
RepLepMass = Solve[MLep[[5]] == CMassElectron/vh && MLep[[7]] == CMassMu/vh && MLep[[9]] == CMassTau/vh, {ye, ymu, ytau}, Reals]

(* Relation between Quark Yukawa couplings and its masses *)
MQuark = HessianH[VF, QuarkBase]/.\[Phi]vev/.rotate;
mQ = Eigenvalues[MQuark]
RepQMass = Solve[mQ[[2]]==CMassBottom/vh && mQ[[4]]==CMassCharm/vh && mQ[[6]] == CMassDown/vh && mQ[[8]] == CMassStrange/vh && mQ[[10]] == CMassTop/vh && mQ[[12]]==CMassUp/vh, {yb, yc, yd, ys, yt, yu}];



Template = CPDark;
TEMPLATE = CPDARK;

(* all counter terms parameters and parameters in the potential are listed here *)
countertermspara = {"ms11CT", "ms22CT", "mssCT", "L1CT", "L2CT", "L3CT", "L4CT", "L5CT", "L6CT", "L7CT", "L8CT", "ArealCT", "AimagCT", "TCBCT", "TCPCT", "T1CT", "T2CT", "TSCT"};
inputparas = {"ms11", "ms22", "mss", "L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "Areal", "Aimag", "vh"};
Cform = {I->II, \[Delta]\[Lambda]1->L1CT, \[Delta]\[Lambda]2->L2CT, \[Delta]\[Lambda]3->L3CT, \[Delta]\[Lambda]4->L4CT, \[Delta]\[Lambda]5->L5CT, \[Delta]\[Lambda]6->L6CT, \[Delta]\[Lambda]7->L7CT, \[Delta]\[Lambda]8->L8CT, \[Delta]Areal->ArealCT, \[Delta]Aimag->AimagCT, 
		\[Delta]ms11->ms11CT, \[Delta]ms22->ms22CT, \[Delta]mss->mssCT, \[Delta]TCB->TCBCT, \[Delta]TCP->TCPCT, \[Delta]T1->T1CT, \[Delta]T2->T2CT, \[Delta]TS->TSCT,
		(*  label for Vtree *)
		\[Lambda]1->L1, \[Lambda]2->L2, \[Lambda]3->L3, \[Lambda]4->L4, \[Lambda]5->L5, \[Lambda]6->L6, \[Lambda]7->L7, \[Lambda]8->L8};

Higgsphysicals = {"G^+","G^-", "H^+", "H^-", "G^0", "A", "h_SM", "h_1", "h_H"};
NNeutralHiggs = 5; 
NChargedHiggs = 4;
NLepton = Length@LepBase;
NQuarks = Length@QuarkBase;
NGauge = Length@Gaugebasis;
NHiggs = Length@Higgsphysicals;

(* \[Phi] = {\[Rho]1, \[Eta]1, \[Rho]2, \[Eta]2, \[Psi]1, \[Psi]2, \[Xi]1,  \[Xi]2, \[Xi]S};
\[Phi]vev = {\[Rho]1->0, \[Eta]1->0, \[Rho]2->\[Omega]CB, \[Eta]2->0, \[Xi]1->\[Omega]1, \[Psi]1->0, \[Xi]2->\[Omega]2, \[Psi]2->\[Omega]CP, \[Xi]S->\[Omega]S};
(*rotate = {\[Omega]1\[Rule]v, \[Omega]2\[Rule]0, \[Omega]CB\[Rule]0, \[Omega]CP\[Rule]0, \[Omega]S\[Rule]0}*)
rotate = {\[Omega]1->vh, \[Omega]2->0, \[Omega]CB->0, \[Omega]CP->0, \[Omega]S->0}; *)

vevsT = {"omega1(T_c)", "omega2(T_c)", "omegaS(T_c)", "omegaCB(T_c)", "omegaCP(T_c)"};
vevs = {"omega1", "omega2", "omegaS", "omegaCB", "omegaCP"};
vevsvalues = {"vh", 0, 0, 0, 0};
VevOrder = {7, 8, 9, 3, 6};

nVEV = Length@vevs; 
nPar = Length@inputparas;
nParCT = Length@countertermspara;

UseVTreeSimplified = false;
UseVCounterSimplified = false;


SetDirectory[NotebookDirectory[]];
str = OpenWrite["ClassTemplate.cpp"];


WriteString[str,
"#include <ext/alloc_traits.h> \n",
"#include <stddef.h> \n",
"#include <algorithm> \n",
"#include <iostream> \n",
"#include <memory> \n",
"#include <BSMPT/models/SMparam.h> \n",
"#include "<>ToString["Eigen/Dense", InputForm]<>" \n",
"#include "<>ToString["Eigen/Eigenvalues", InputForm]<>" \n",
"#include "<>ToString["Eigen/IterativeLinearSolvers", InputForm]<>" \n\n",
"#include <BSMPT/models/"<>ToString[Template]<>".h> \n",
"#include <BSMPT/models/IncludeAllModels.h> \n",
"#include <BSMPT/utility.h> \n",
"using namespace std; \n"
"using namespace Eigen; \n",
"namespace BSMPT{ \n",
"namespace Models{ \n",
"Class_"<>ToString[Template]<>"::Class_"<>ToString[Template]<>" () \n",
"{ \n",
"    Model = ModelID::ModelIDs::"<>ToString[TEMPLATE]<>"; \n",
"    NNeutralHiggs ="<>ToString[NNeutralHiggs]<>"; \n",
"    NChargedHiggs = "<>ToString[NChargedHiggs]<>"; \n\n",
"    NLepton = "<>ToString[NLepton]<>";\n",
"    NQuarks = "<>ToString[NQuarks]<>";\n",
"    NGauge = "<>ToString[NGauge]<>";\n",
"    NHiggs = "<>ToString[NHiggs]<>";\n",

"    nPar = "<>ToString[nPar]<>"; \n",
"    nParCT = "<>ToString[nParCT]<>"; \n\n",

"    nVEV = "<>ToString[nVEV]<>"; \n",
"    VevOrder.resize("<>ToString[nVEV]<>"); \n"];
Do[
WriteString[str,
"    VevOrder["<>ToString[i-1]<>"] = "<>ToString[VevOrder[[i]]]<>"; \n"]
, {i,1,Length@VevOrder}];

WriteString[str,
"    UseVTreeSimplified = "<>ToString[UseVTreeSimplified]<>"; \n",
"    UseVCounterSimplified = "<>ToString[UseVCounterSimplified]<>"; \n",
"} \n\n"];


WriteString[str,
(* starting labels Legends for counter-terms *)
"Class_"<>ToString[Template]<>"::~Class_"<>ToString[Template]<>" (){} \n",
"std::vector<std::string> Class_"<>ToString[Template]<>"::addLegendCT() const{\n",
"    std::vector<std::string> labels;\n"
];

Do[
WriteString[str,
"    labels.push_back("<>ToString[countertermspara[[i]], InputForm]<>"); \n"]
, {i,1, Length@countertermspara}];
WriteString[str, 
"    return labels;\n",
"} \n\n"];

(* labels for Tc, vc, xi_c, vevs at Temperature *)
WriteString[str, 
"std::vector<std::string> Class_"<>ToString[Template]<>"::addLegendTemp() const{\n",
"    std::vector<std::string> labels; \n",
"    labels.push_back("<>ToString["T_c", InputForm]<>"); \n",
"    labels.push_back("<>ToString["v_c", InputForm]<>"); \n",
"    labels.push_back("<>ToString["v_c/T_c", InputForm]<>"); \n"]
Do[
WriteString[str,
"    labels.push_back("<>ToString[vevsT[[i]], InputForm]<>"); \n"];
,{i,1 Length@vevsT}]
WriteString[str,
"    return labels; \n",
"}\n\n",
(* Tripple H coupling -- lables for physical Higgs fields after mixing *)
"std::vector<std::string> Class_"<>ToString[Template]<>"::addLegendTripleCouplings() const{ \n",
"    std::vector<std::string> labels; \n",
"    std::vector<std::string> particles; \n",
"    particles.resize(NHiggs); \n"];

Do[
WriteString[str,
"    particles.push_back("<>ToString[Higgsphysicals[[i]], InputForm]<>");\n"];
,{i,1, Length@Higgsphysicals}];

WriteString[str,
"    std::string out ="<>ToString["Tree_", InputForm]<>";\n",
"    for(std::size_t i=0; i<NHiggs;i++)\n",
"    {\n",
"        for(std::size_t j=i;j<NHiggs;j++)\n",
"        {\n",
"            for(std::size_t k=j; k<NHiggs;j++) \n",
"            {\n",
"                labels.push_back("<>ToString["Tree_", InputForm]<>"+particles.at(i)+particles.at(j)+particles.at(k));\n",
"                labels.push_back("<>ToString["CT_", InputForm]<>"+particles.at(i)+particles.at(j)+particles.at(k));\n",
"                labels.push_back("<>ToString["CW_", InputForm]<>"+particles.at(i)+particles.at(j)+particles.at(k));\n",
"            }\n",
"        }\n",
"    }\n",
"    return labels;\n",
"}\n\n",

(* vevs at T=0, vevs labels *)
"std::vector<std::string> Class_"<>ToString[Template]<>"::addLegendVEV() const{\n",
"    std::vector<std::string> labels; \n"];

Do[
WriteString[str,
"    labels.push_back("<>ToString[vevs[[i]], InputForm]<>");\n"
];
, {i,1, Length@vevs}];

WriteString[str,
"    return labels; \n",
"}\n\n",

(* read input from files into sstream *)
"void Class_"<>ToString[Template]<>"::ReadAndSet(const std::string& linestr, std::vector<double>& par )\n",
"{\n",
"    std::stringstream ss(linestr); \n",
"    double tmp;\n",
"    if (UseIndexCol){\n",
"        ss >> tmp;",
"    }\n",
"    for(int k=1; k<="<> ToString[Length@countertermspara]<>";k++)\n",
"    {\n",
"        ss >> tmp;\n",
"        if(k==1) par[0] = tmp; \n"];

Do[
WriteString[str,
"        else if(k=="<>ToString[i]<>") par["<>ToString[i-1]<>"] = tmp; \n"
];
, {i, 2, Length@inputparas}]
WriteString[str,
"    }\n"];
(*Do[
WriteString[str,
"    par["<>ToString[i-1]<>"] = "<>ToString[inputparas[[i]]]<>";\n"
]
,{i, 1, Length@inputparas}];
*)
WriteString[str,
"    set_gen(par);\n",
"    return;\n",
"}\n\n"];

(* set input parameters, including dependent and independent parameters *)
WriteString[str,
"void Class_"<>ToString[Template]<>"::set_gen(const std::vector<double>& par) {\n"];

Do[
WriteString[str, 
"    "<>ToString[inputparas[[i]]]<>" = par["<>ToString[i-1]<>"];\n"
]
,{i,1, Length@inputparas}];
WriteString[str,
"    scale = C_vev0;\n",
"    vevTreeMin.resize(nVEV);\n"];
(* other dependent parameters from minimize potential are here *)
WriteString[str, "    L1 = C_MassSMHiggs*C_MassSMHiggs/(C_vev0*C_vev0);\n", 
"    ms11 = - L1*vh*vh/0.2e1;\n"];

Do[
WriteString[str, 
"    vevTreeMin["<>ToString[i-1]<>"] =  "<>ToString[vevsvalues[[i]]]<>";\n"]
,{i,1, Length@vevsvalues}];

WriteString[str,
"    vevTree.resize(NHiggs);\n",
"    vevTree = MinimizeOrderVEV(vevTreeMin);\n",
"    if(!SetCurvatureDone) SetCurvatureArrays();\n",
"}\n\n"];


(* setting counter-terms --\[Rule] input counterterm or fixed label for counter-terms in the calulations *)
WriteString[str,
"void Class_"<>ToString[Template]<>"::set_CT_Pot_Par(const std::vector<double>& par){\n"];
Do[
WriteString[str,
"    "<>ToString[countertermspara[[i]]]<>" = par["<>ToString[i-1]<>"]; \n"]
,{i,1, Length@countertermspara}]


(* need to adjust the cform to get better C++ form for numerical evaluation *)
(* need to imrove cform*)
Do[
temp = (D[VCT[[1]], \[Phi][[i]]]/.\[Phi]vev/.rotate0//Simplify)/.Cform;
If[temp === 0, Continue, WriteString[str, "    Curvature_Higgs_CT_L1["<>ToString[i-1, InputForm]<>"]="<>ToString[temp//CForm]<>";\n"]];
, {i,1 , NHiggs}];

WriteString[str, "\n"];

Do[
temp =(D[D[VCT[[1]], \[Phi][[i]]], \[Phi][[j]]]/.\[Phi]vev/.rotate0//Simplify)/.Cform;
If[temp === 0, Continue, WriteString[str, "    Curvature_Higgs_CT_L2["<>ToString[i-1]<>"]["<>ToString[j-1]<>"]=",ToString[temp//CForm]<>";\n"]];
, {i,1, NHiggs}, {j,i, NHiggs}];

WriteString[str, "\n"];

Do[
temp =  (D[D[D[VCT[[1]], \[Phi][[i]]], \[Phi][[j]]], \[Phi][[k]]]/.\[Phi]vev/.rotate0//Simplify)/.Cform;
If[temp ===0, Continue, WriteString[str, "    Curvature_Higgs_CT_L3["<>ToString[i-1]<>"]["<>ToString[j-1]<>"]["<>ToString[k-1]<>"] ="<> ToString[temp//CForm]<>";\n"]];
, {i,1, NHiggs}, {j,i, NHiggs}, {k,j, NHiggs}];

WriteString[str, "\n"];

Do[
temp = (D[D[D[D[VCT[[1]], \[Phi][[i]]], \[Phi][[j]]], \[Phi][[k]]], \[Phi][[l]]]/.\[Phi]vev/.rotate0//Simplify)/.Cform;
If[temp === 0, Continue, WriteString[str,"    Curvature_Higgs_CT_L4["<>ToString[i-1, InputForm]<>"]["<>ToString[j-1, InputForm]<>"]["<>ToString[k-1]<>"]["<>ToString[l-1, InputForm]<>"]="<>ToString[ temp//CForm]<>";\n"]];
Clear[temp]
, {i,1, NHiggs}, {j,i, NHiggs}, {k,j, NHiggs}, {l,k, NHiggs}];

WriteString[str, "\n"];
WriteString[str, 
"    for(std::size_t k1=0;k1<NHiggs;k1++)\n",
"    {\n",
"        for(std::size_t k2=k1;k2<NHiggs;k2++)\n",
"        {\n",
"            Curvature_Higgs_CT_L2[k2][k1] = Curvature_Higgs_CT_L2[k1][k2];\n",
"            for(std::size_t k3=k2;k3<NHiggs;k3++)\n",
"            {\n",
"                Curvature_Higgs_CT_L3[k1][k3][k2] = Curvature_Higgs_CT_L3[k1][k2][k3];\n",
"                Curvature_Higgs_CT_L3[k2][k1][k3]= Curvature_Higgs_CT_L3[k1][k2][k3];\n",
"                Curvature_Higgs_CT_L3[k2][k3][k1] = Curvature_Higgs_CT_L3[k1][k2][k3];\n",
"                Curvature_Higgs_CT_L3[k3][k1][k2] = Curvature_Higgs_CT_L3[k1][k2][k3];\n",
"                Curvature_Higgs_CT_L3[k3][k2][k1] = Curvature_Higgs_CT_L3[k1][k2][k3];\n",
"                for(std::size_t k4=k3;k4<NHiggs;k4++)\n",
"                {\n",
"                    Curvature_Higgs_CT_L4[k2][k3][k4][k1] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_CT_L4[k3][k4][k1][k2] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_CT_L4[k4][k1][k2][k3] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_CT_L4[k2][k1][k3][k4] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_CT_L4[k4][k2][k1][k3] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_CT_L4[k3][k4][k2][k1] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_CT_L4[k1][k3][k4][k2] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_CT_L4[k3][k2][k1][k4] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_CT_L4[k4][k3][k2][k1] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_CT_L4[k1][k4][k3][k2] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_CT_L4[k2][k1][k4][k3] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_CT_L4[k4][k2][k3][k1] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_CT_L4[k1][k4][k2][k3] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_CT_L4[k3][k1][k4][k2] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_CT_L4[k2][k3][k1][k4] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_CT_L4[k1][k3][k2][k4] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_CT_L4[k4][k1][k3][k2] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_CT_L4[k2][k4][k1][k3] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_CT_L4[k3][k2][k4][k1] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_CT_L4[k1][k2][k4][k3] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_CT_L4[k3][k1][k2][k4] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_CT_L4[k4][k3][k1][k2] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_CT_L4[k2][k4][k3][k1] = Curvature_Higgs_CT_L4[k1][k2][k3][k4];\n",
"                }\n",
"            }\n",
"        }\n",
"    }\n",
"    return;\n"
];
WriteString[str, "}\n\n"];


WriteString[str, "void Class_"<>ToString[Template]<>"::write() const { \n",
"    typedef std::numeric_limits<double> dbl;\n",
"    std::cout.precision(dbl::max_digits10);\n",
"    std::cout << "<>ToString["Model = ", InputForm]<>" << Model << "<>ToString[" \n ", InputForm]<>";\n",
"    std::cout << "<>ToString["The parameters are : \n", InputForm]<>";\n",
"    std::cout <<"<>ToString["Renorm Scale = ", InputForm]<>" << scale << "<>ToString["\n", InputForm]<>";\n"];

Do[
WriteString[str,
"    std::cout <<"<>ToString[StringJoin[inputparas[[i]], " = "], InputForm]<>" << "<>ToString[inputparas[[i]]]<>" << std::endl;\n"]
,{i,1, Length@inputparas}];

WriteString[str,
"    std::cout << "<>ToString["The counterterm parameters are :  \n", InputForm]<>";\n"];
Do[
WriteString[str,
"    std::cout << "<>ToString[StringJoin[countertermspara[[i]], " = "], InputForm]<>" << "<>countertermspara[[i]]<> " << std::endl;\n"]
, {i,1, Length@countertermspara}]

(* not going further for matrix rotation --> use scanMS project *)

WriteString[str, "}\n\n\n"]


<< counterterms`


WriteString[str, "std::vector<double>Class_"<>ToString[Template]<>"::calc_CT() const {\n",
    "    std::vector<double> parCT;\n",
    "    if(!SetCurvatureDone){\n",
    "        std::string retmes = __func__;\n",
    "        retmes += "<>ToString["was called before SetCurvatureArrays()!\n"]<>"\n";
    "        throw std::runtime_error(retmes);\n",
    "    }\n",
    "    if(!CalcCouplingsdone){\n",
    "        std::string retmes = __func__;\n",
    "        throw std::runtime_error(retmes);\n",
    "    }\n",
    "    std::vector<double> WeinbergNabla, WeinbergHesse;\n",
    "    WeinbergNabla = WeinbergFirstDerivative();\n",
    "    WeinbergHesse = WeinbergSecondDerivative();\n",
    "    VectorXd NablaWeinberg(NHiggs);\n",
    "    MatrixXd HesseWeinberg(NHiggs, NHiggs), HiggsRot(NHiggs, NHiggs);\n",
    "    for(std::size_t i=0;i<NHiggs;i++){\n",
    "        NablaWeinberg[i] = WeinbergNabla[i];\n",
    "        for(std::size_t j=0;j<NHiggs;j++) HesseWeinberg(i,j) = WeinbergHesse.at(j*NHiggs+i);\n",
    "    }\n\n"
];
Do[
WriteString[str, "    parCT.push_back((double)("<>ToString[solutions[[1]][[i]][[2]]//CForm]<>"));\n"]
,{i, 1, Length@solutions[[1]]}]
WriteString[str, "    return parCT;\n",
"}\n"];


(* need the triplehiggs for EWSB running, cannot make without it *)
WriteString[str, "void Class_"<>ToString[Template]<>"::TripleHiggsCouplings()\n",
"{\n",
"	if(!SetCurvatureDone)SetCurvatureArrays();\n",
"	if(!CalcCouplingsdone)CalculatePhysicalCouplings();\n",
"	std::vector<double> HiggsOrder(NHiggs);\n",
"	// Here you have to set the vector HiggsOrder. By telling e.g. HiggsOrder[0] = 5 you always want your 6th lightest\n",
"	// particle to be the first particle in the vector (which has the index 5 because they are sorted by mass)\n",
"	// example for keeping the mass order\n",
"	for(std::size_t i=0;i<NHiggs;i++) {\n",
"		HiggsOrder[i]=i;\n",
"	}\n",
"	std::vector<double> TripleDeriv;\n",
"    TripleDeriv=WeinbergThirdDerivative();\n"
"	std::vector<std::vector<std::vector<double>>> GaugeBasis(NHiggs, std::vector<std::vector<double>>(NHiggs,std::vector<double>(NHiggs)));\n",
"	for(std::size_t i=0;i<NHiggs;i++)\n",
"	  {\n",
"		for(std::size_t j=0;j<NHiggs;j++)\n",
"		{\n",
"		  for(std::size_t k=0;k<NHiggs;k++)\n",
"			{\n",
"			  GaugeBasis[i][j][k] = TripleDeriv.at(i+j*NHiggs+k*NHiggs*NHiggs);\n",
"			}\n",
"		}\n",
"	  }\n",
"	MatrixXd HiggsRot(NHiggs,NHiggs);\n",
"	for(std::size_t i=0;i<NHiggs;i++)\n",
"	{\n",
"		for(std::size_t j=0;j<NHiggs;j++)\n",
"		{\n",
"			HiggsRot(i,j) = HiggsRotationMatrix[i][j];\n",
"		}\n",
"	}\n",
"	MatrixXd HiggsRotSort(NHiggs,NHiggs);\n",
"	for(std::size_t i=0;i<NHiggs;i++)\n",
"	{\n",
"		HiggsRotSort.row(i) = HiggsRot.row(HiggsOrder[i]);\n",
"	}\n",
"	TripleHiggsCorrectionsCWPhysical.resize(NHiggs);\n",
"	TripleHiggsCorrectionsTreePhysical.resize(NHiggs);\n",
"	TripleHiggsCorrectionsCTPhysical.resize(NHiggs);\n",
"	for(std::size_t i=0;i<NHiggs;i++) {\n",
"		TripleHiggsCorrectionsTreePhysical[i].resize(NHiggs);\n",
"		TripleHiggsCorrectionsCWPhysical[i].resize(NHiggs);\n",
"		TripleHiggsCorrectionsCTPhysical[i].resize(NHiggs);\n",
"		for(std::size_t j=0;j<NHiggs;j++) {\n",
"			TripleHiggsCorrectionsCWPhysical[i][j].resize(NHiggs);\n",
"			TripleHiggsCorrectionsTreePhysical[i][j].resize(NHiggs);\n",
"			TripleHiggsCorrectionsCTPhysical[i][j].resize(NHiggs);\n",
"		}\n",
"	}\n",
"	for(std::size_t i=0;i<NHiggs;i++)\n",
"	{\n",
"		for(std::size_t j=0;j<NHiggs;j++)\n",
"		{\n",
"			for(std::size_t k=0;k<NHiggs;k++)\n",
"			{\n",
"			    TripleHiggsCorrectionsCWPhysical[i][j][k] = 0;\n",
"			    TripleHiggsCorrectionsTreePhysical[i][j][k] = 0;\n",
"			    TripleHiggsCorrectionsCTPhysical[i][j][k] = 0;\n",
"			    for(std::size_t l=0;l<NHiggs;l++)\n",
"			    {\n",
"				    for(std::size_t m=0;m<NHiggs;m++)\n",
"				    {\n",
"					    for(std::size_t n=0;n<NHiggs;n++)\n",
"					    {\n",
"						     double RotFac = HiggsRotSort(i,l)*HiggsRotSort(j,m)*HiggsRotSort(k,n);\n",
"						     TripleHiggsCorrectionsCWPhysical[i][j][k] += RotFac*GaugeBasis[l][m][n];\n",
"						     TripleHiggsCorrectionsTreePhysical[i][j][k] += RotFac*LambdaHiggs_3[l][m][n];\n",
"						     TripleHiggsCorrectionsCTPhysical[i][j][k] += RotFac*LambdaHiggs_3_CT[l][m][n];\n",
"					    }\n",
"				    }\n",
"			    }\n",
"			}\n",
"		}\n",
"	}\n",
"}\n"];





(* no need to use triple Higgs coupling function now *)
(* The most important function, curvaturearray *)
WriteString[str, 
"void Class_"<>ToString[Template]<>"::SetCurvatureArrays(){\n",
"    initVectors();\n",
"    SetCurvatureDone=true;\n",
"    for(std::size_t i=0; i<NHiggs;i++) HiggsVev[i] = vevTree[i];\n"
];
(* curvature_Higgs *)
Do[
temp = (D[VTree[[1]], \[Phi][[i]]]/.\[Phi]vev/.rotate0//Simplify)/.Cform;
If[temp === 0, Continue, WriteString[str, "    Curvature_Higgs_L1["<>ToString[i-1, InputForm]<>"]="<>ToString[temp//CForm]<>";\n"]]
, {i,1, NHiggs}];

WriteString[str, "\n"];

Do[
temp = (D[VTree[[1]], \[Phi][[i]], \[Phi][[j]]]/.\[Phi]vev/.rotate0//Simplify)/.Cform;
If[temp === 0, Continue, WriteString[str, "    Curvature_Higgs_L2["<>ToString[i-1]<>"]["<>ToString[j-1]<>"]=",ToString[temp//CForm]<>";\n"]]
,{i,1, NHiggs}, {j,i, NHiggs}];

WriteString[str, "\n"];

Do[
temp = (D[VTree[[1]], \[Phi][[i]], \[Phi][[j]], \[Phi][[k]]]/.\[Phi]vev/.rotate0//Simplify)/.Cform;
If[temp ===0, Continue, WriteString[str, "    Curvature_Higgs_L3["<>ToString[i-1]<>"]["<>ToString[j-1]<>"]["<>ToString[k-1]<>"] ="<> ToString[temp//CForm]<>";\n"]]
, {i, 1, NHiggs}, {j,i, NHiggs}, {k,j, NHiggs}];

WriteString[str, "\n"];

Do[
temp = (D[VTree[[1]], \[Phi][[i]], \[Phi][[j]], \[Phi][[k]], \[Phi][[l]]]/.\[Phi]vev/.rotate0//Simplify)/.Cform;
If[temp === 0, Continue, WriteString[str,"    Curvature_Higgs_L4["<>ToString[i-1, InputForm]<>"]["<>ToString[j-1, InputForm]<>"]["<>ToString[k-1]<>"]["<>ToString[l-1, InputForm]<>"]="<>ToString[ temp//CForm]<>";\n"]]
,{i, 1, NHiggs}, {j,i, NHiggs}, {k,j, NHiggs}, {l,k, NHiggs}];
Clear[temp];

WriteString[str, "\n"];
WriteString[str, 
"    for(std::size_t k1=0;k1<NHiggs;k1++)\n",
"    {\n",
"        for(std::size_t k2=k1;k2<NHiggs;k2++)\n",
"        {\n",
"            Curvature_Higgs_L2[k2][k1] = Curvature_Higgs_L2[k1][k2];\n",
"            for(std::size_t k3=k2;k3<NHiggs;k3++)\n",
"            {\n",
"                Curvature_Higgs_L3[k1][k3][k2] = Curvature_Higgs_L3[k1][k2][k3];\n",
"                Curvature_Higgs_L3[k2][k1][k3]= Curvature_Higgs_L3[k1][k2][k3];\n",
"                Curvature_Higgs_L3[k2][k3][k1] = Curvature_Higgs_L3[k1][k2][k3];\n",
"                Curvature_Higgs_L3[k3][k1][k2] = Curvature_Higgs_L3[k1][k2][k3];\n",
"                Curvature_Higgs_L3[k3][k2][k1] = Curvature_Higgs_L3[k1][k2][k3];\n",
"                for(std::size_t k4=k3;k4<NHiggs;k4++)\n",
"                {\n",
"                    Curvature_Higgs_L4[k2][k3][k4][k1] = Curvature_Higgs_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_L4[k3][k4][k1][k2] = Curvature_Higgs_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_L4[k4][k1][k2][k3] = Curvature_Higgs_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_L4[k2][k1][k3][k4] = Curvature_Higgs_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_L4[k4][k2][k1][k3] = Curvature_Higgs_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_L4[k3][k4][k2][k1] = Curvature_Higgs_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_L4[k1][k3][k4][k2] = Curvature_Higgs_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_L4[k3][k2][k1][k4] = Curvature_Higgs_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_L4[k4][k3][k2][k1] = Curvature_Higgs_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_L4[k1][k4][k3][k2] = Curvature_Higgs_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_L4[k2][k1][k4][k3] = Curvature_Higgs_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_L4[k4][k2][k3][k1] = Curvature_Higgs_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_L4[k1][k4][k2][k3] = Curvature_Higgs_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_L4[k3][k1][k4][k2] = Curvature_Higgs_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_L4[k2][k3][k1][k4] = Curvature_Higgs_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_L4[k1][k3][k2][k4] = Curvature_Higgs_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_L4[k4][k1][k3][k2] = Curvature_Higgs_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_L4[k2][k4][k1][k3] = Curvature_Higgs_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_L4[k3][k2][k4][k1] = Curvature_Higgs_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_L4[k1][k2][k4][k3] = Curvature_Higgs_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_L4[k3][k1][k2][k4] = Curvature_Higgs_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_L4[k4][k3][k1][k2] = Curvature_Higgs_L4[k1][k2][k3][k4];\n",
"                    Curvature_Higgs_L4[k2][k4][k3][k1] = Curvature_Higgs_L4[k1][k2][k3][k4];\n",
"                }\n",
"            }\n",
"         }\n",
"    }\n"
]
WriteString[str, "\n\n"];

(* Curvature_Gauge *)
Do[
	temp = D[Vgauge, Gaugebasis[[a]], Gaugebasis[[b]], \[Phi][[i]], \[Phi][[j]]]/.Cform;
	If[temp === 0, Continue, WriteString[str, "    Curvature_Gauge_G2H2["<>ToString[a-1]<>"]["<>ToString[b-1]<>"]["<>ToString[i-1]<>"]["<>ToString[j-1]<>"]="<>cform[temp]<>";\n"]]
,{a,1, NGauge}, {b, 1, NGauge}, {i, 1, NHiggs}, {j,1,NHiggs}]
Clear[temp];
WriteString[str, "\n\n"];

(* Curvature_Lepton *)
WriteString[str, "    std::complex<double> A(0,1);\n"]
Do[
	temp = D[VFLep[[1]], LepBase[[a]], LepBase[[b]], \[Phi][[i]]]/.RepLepMass;
	If[temp[[1]] === 0, Continue, WriteString[str, "    Curvature_Lepton_F2H1["<>ToString[a-1]<>"]["<>ToString[b-1]<>"]["<>ToString[i-1]<>"]="<>cform[temp[[1]]]<>";\n"]]
,{a, 1, NLepton}, {b, 1, NLepton}, {i, 1, NHiggs}]

Clear[temp];
WriteString[str, "\n\n"];

(* Curvature_Quark *)
WriteString[str,
"    std::complex<double> V11, V12, V13, V21, V22, V23, V31, V32, V33;\n",
"    V11 = C_Vud;\n",
"    V12 = C_Vus;\n",
"    V13 = C_Vub;\n",
"    V21 = C_Vcd;\n",
"    V22 = C_Vcs;\n",
"    V23 = C_Vcb;\n",
"    V31 = C_Vtd;\n",
"    V32 = C_Vts;\n",
"    V33 = C_Vtb;\n"
];
Do[
	temp = D[VF, QuarkBase[[a]], QuarkBase[[b]], \[Phi][[i]]]/.RepQMass;
	If[temp[[1]] === 0, Continue, WriteString[str, "    Curvature_Quark_F2H1["<>ToString[a-1]<>"]["<>ToString[b-1]<>"]["<>ToString[i-1]<>"]="<>cform[temp[[1]]]<>";\n"]]
,{a, 1, NQuarks}, {b, 1, NQuarks}, {i, 1, NHiggs}]

Clear[temp];

WriteString[str, 
"    SetCurvatureDone = true;\n",
"}\n\n"]


WriteString[str, "bool Class_"<>ToString[Template]<>"::CalculateDebyeSimplified(){\n",
"    return false;\n",
"}\n\n"]


WriteString[str, "bool Class_"<>ToString[Template]<>"::CalculateDebyeGaugeSimplified(){\n",
"    return false;\n",
"}\n\n"]


WriteString[str, "double Class_"<>ToString[Template]<>"::VTreeSimplified(const std::vector<double>& v) const\n",
"{\n",
"    (void) v;\n",
"    return 0;\n",
"}\n\n"
]


WriteString[str, "double Class_"<>ToString[Template]<>"::VCounterSimplified(const std::vector<double>& v) const\n",
"{\n",
"    (void) v;\n",
"    return 0;\n",
"}\n\n"
]


WriteString[str, "void Class_"<>ToString[Template]<>"::Debugging(const std::vector<double>& input, std::vector<double>& output) const{\n",
"    (void) input;\n",
"    (void) output;\n",
"    std::cout<<"<>ToString["Start", InputForm]<>";\n",
"    bool Debug = true;\n",
"    if(not Debug) return;\n",
"    write();\n",
"    if(!SetCurvatureDone){\n",
"        std::string retmes = __func__;\n",
"        retmes += "<>ToString[" was called before SetCurvatureArrays()!\n", InputForm]<>";\n",
"        throw std::runtime_error(retmes);\n",
"    }\n",
"    if(!CalcCouplingsdone){\n",
"        std::string retmes = __func__;\n",
"        retmes += "<>ToString["was called before CalculatePhysicalCouplings()!\n", InputForm]<>";\n",
"        throw std::runtime_error(retmes);\n",
"    }",
"        std::vector<double> WeinbergNabla,WeinbergHesse;\n",
"    WeinbergNabla=WeinbergFirstDerivative();\n",
"    WeinbergHesse=WeinbergSecondDerivative();\n",
"        VectorXd NablaWeinberg(NHiggs);\n",
"        MatrixXd HesseWeinberg(NHiggs,NHiggs),HiggsRot(NHiggs,NHiggs);\n",
"        for(std::size_t i=0;i<NHiggs;i++)\n",
"        {\n",
"                NablaWeinberg[i] = WeinbergNabla[i];\n",
"                for(std::size_t j=0;j<NHiggs;j++) HesseWeinberg(i,j) = WeinbergHesse.at(j*NHiggs+i);\n",
"        }\n",
"        std::cout << "<>ToString["Start ", InputForm]<>" << std::endl;\n",
"\n",
"        std::vector<double> LOMasses(NHiggs), NLOMasses(NHiggs);\n",
"    LOMasses=HiggsMassesSquared(vevTree,0);\n",
"        MatrixXd TreeMassMatrix(NHiggs,NHiggs), CTMassMatrix(NHiggs,NHiggs), CWMassMatrix(NHiggs,NHiggs);\n",
"        for(std::size_t i=0;i<NHiggs;i++)\n",
"        {\n",
"                for(std::size_t j=0;j<NHiggs;j++){\n",
"                        TreeMassMatrix(i,j) = 0;\n",
"                        CTMassMatrix(i,j) = 0;\n",
"                        CWMassMatrix(i,j) = 0;\n",
"                }\n",
"        }\n",
"        for(std::size_t i=0;i<NHiggs;i++)\n",
"        {\n",
"                for(std::size_t j=0;j<NHiggs;j++)\n",
"                {\n",
"                        TreeMassMatrix(i,j) = Curvature_Higgs_L2[i][j];\n",
"                        for(std::size_t k=0;k<NHiggs;k++)\n",
"                        {\n",
"                                TreeMassMatrix(i,j) += Curvature_Higgs_L3[i][j][k]*vevTree[k];\n",
"                                for(std::size_t l=0;l<NHiggs;l++) TreeMassMatrix(i,j) += 0.5*Curvature_Higgs_L4[i][j][k][l] * vevTree[k]*vevTree[l];\n",
"                        }\n",
"                        CTMassMatrix(i,j) = Curvature_Higgs_CT_L2[i][j];\n",
"                        for(std::size_t k=0;k<NHiggs;k++)\n",
"                        {\n",
"                                CTMassMatrix(i,j) += Curvature_Higgs_CT_L3[i][j][k]*vevTree[k];\n",
"                                for(std::size_t l=0;l<NHiggs;l++) CTMassMatrix(i,j) += 0.5*Curvature_Higgs_CT_L4[i][j][k][l] * vevTree[k]*vevTree[l];\n",
"                        }\n",
"                }\n",
"        }\n",
"        CWMassMatrix = HesseWeinberg;\n",
"        SelfAdjointEigenSolver<MatrixXd> es(TreeMassMatrix+CTMassMatrix+CWMassMatrix,EigenvaluesOnly);\n",
"        for(std::size_t i=0;i<NHiggs;i++) {\n",
"                NLOMasses[i] = es.eigenvalues()[i];\n",
"        }\n",
"        std::cout << "<>ToString["Print LO | NLO Mass^2 ", InputForm]<>" << std::endl;\n",
"        for(std::size_t i=0;i<NHiggs;i++){\n",
"                std::cout << LOMasses[i] << "<>ToString[" | ", InputForm]<>" << NLOMasses[i] << std::endl;\n",
"        }\n",
(* added consistency condition here to check numerical issue, cannot be automize yet *)
"}\n"
]



WriteString[str, "}\n",
"}\n"
];



Close[str];
Clear[str];



NotebookDirectory[]
