(* ::Package:: *)

pu


(* ::Input:: *)
(*(* Renormalization sheme *)*)
(*(* Using NLO mass = LO mass in the minimzation of effective CW_potential with vev configuration at T=0 *)*)
(*(* Need RREFpivot modules - Gaus partial pivoting elemination method *)*)


(* ::Input:: *)
(**)


(* field configurations-input*)
\[CapitalPhi]1 = (1/Sqrt[2]){{\[Rho]1 + I \[Eta]1}, {\[Xi]1  + I \[Psi]1}};
\[CapitalPhi]2 = (1/Sqrt[2]){{\[Rho]2  + I \[Eta]2}, {\[Xi]2 + I \[Psi]2 }};
\[CapitalPhi]1T = (ConjugateTranspose[\[CapitalPhi]1]/.{Conjugate[\[Eta]1]->\[Eta]1, Conjugate[\[Rho]1]->\[Rho]1, Conjugate[\[Xi]1]->\[Xi]1, Conjugate[\[Psi]1]->\[Psi]1, Conjugate[\[Xi]1 + \[Omega]1]->\[Xi]1 + \[Omega]1})//Flatten//Simplify;
\[CapitalPhi]2T =  (ConjugateTranspose[\[CapitalPhi]2]/.{Conjugate[\[Eta]2]->\[Eta]2, Conjugate[\[Rho]2 ]->\[Rho]2 , Conjugate[\[Xi]2]->\[Xi]2, Conjugate[\[Psi]2]->\[Psi]2, Conjugate[\[Xi]2 + \[Omega]2]->\[Xi]2 + \[Omega]2, Conjugate[\[Psi]2 + \[Omega]CP]->\[Psi]2 + \[Omega]CP})//Flatten//Simplify;
\[CapitalPhi]S = \[Xi]S;

\[CapitalPhi]real = {\[Rho]1 \[Element] Reals, \[Eta]1 \[Element] Reals, \[Xi]1 \[Element] Reals, \[Psi]1 \[Element] Reals, \[Rho]2 \[Element] Reals, \[Eta]2 \[Element] Reals, \[Xi]2 \[Element] Reals, \[Psi]2 \[Element] Reals};


(* set of propest counter-terms input *)
countertermspara = {\[Delta]ms11, \[Delta]ms22, \[Delta]mss, \[Delta]Areal ,\[Delta]Aimag, \[Delta]\[Lambda]1,\[Delta]\[Lambda]2,\[Delta]\[Lambda]3, \[Delta]\[Lambda]4, \[Delta]\[Lambda]5,\[Delta]\[Lambda]6, \[Delta]\[Lambda]7, \[Delta]\[Lambda]8, \[Delta]T1, \[Delta]T2, \[Delta]TCP, \[Delta]TCB, \[Delta]TS};


Vtemp = \[Delta]T1 (\[Xi]1 + \[Omega]1) + \[Delta]T2 (\[Xi]2 + \[Omega]2) + \[Delta]TCP (\[Psi]2 + \[Omega]CP) + \[Delta]TCB (\[Rho]2 + \[Omega]CB) + \[Delta]TS (\[Xi]S + \[Omega]S);
VCT  = \[Delta]ms11 \[CapitalPhi]1T.\[CapitalPhi]1 + \[Delta]ms22 \[CapitalPhi]2T.\[CapitalPhi]2  + (\[Delta]\[Lambda]1/2)(\[CapitalPhi]1T.\[CapitalPhi]1)^2 + (\[Delta]\[Lambda]2/2)(\[CapitalPhi]2T.\[CapitalPhi]2)^2  + \[Delta]\[Lambda]3 (\[CapitalPhi]1T.\[CapitalPhi]1) (\[CapitalPhi]2T.\[CapitalPhi]2) + \[Delta]\[Lambda]4  (\[CapitalPhi]1T.\[CapitalPhi]2) (\[CapitalPhi]2T.\[CapitalPhi]1) + (\[Delta]Areal + I \[Delta]Aimag)\[CapitalPhi]1T.\[CapitalPhi]2 \[CapitalPhi]S + (\[Delta]Areal - I \[Delta]Aimag) \[CapitalPhi]2T.\[CapitalPhi]1 \[CapitalPhi]S + (\[Delta]\[Lambda]5/2)( (\[CapitalPhi]1T.\[CapitalPhi]2)^2 +  (\[CapitalPhi]2T.\[CapitalPhi]1)^2) +\[Delta]mss \[CapitalPhi]S^2/2 + (\[Delta]\[Lambda]6/4)\[CapitalPhi]S^4 + (\[Delta]\[Lambda]7/2)(\[CapitalPhi]1T.\[CapitalPhi]1)\[CapitalPhi]S^2 +(\[Delta]\[Lambda]8/2)(\[CapitalPhi]2T.\[CapitalPhi]2)\[CapitalPhi]S^2 + Vtemp//Simplify


(* Tree level configurations of fields and vev-input *)
\[Phi] = {\[Rho]1, \[Eta]1, \[Rho]2, \[Eta]2, \[Psi]1, \[Psi]2, \[Xi]1,  \[Xi]2, \[Xi]S};
\[Phi]vev = {\[Rho]1->0, \[Eta]1->0, \[Rho]2->\[Omega]CB, \[Eta]2->0, \[Xi]1->\[Omega]1, \[Psi]1->0, \[Xi]2->\[Omega]2, \[Psi]2->\[Omega]CP, \[Xi]S->\[Omega]S};
(*rotate = {\[Omega]1\[Rule]v, \[Omega]2\[Rule]0, \[Omega]CB\[Rule]0, \[Omega]CP\[Rule]0, \[Omega]S\[Rule]0}*)
rotate = {\[Omega]1->vh, \[Omega]2->0, \[Omega]CB->0, \[Omega]CP->0, \[Omega]S->0};
rotate0 = {\[Omega]1->0, \[Omega]2->0, \[Omega]CB->0, \[Omega]CP->0, \[Omega]S->0};




firstderivatives  = ConstantArray[0, Length[\[Phi]]];
secondderivatives = ConstantArray[0, {Length[\[Phi]], Length[\[Phi]]}];
thirdderivatives = ConstantArray[0, {Length[\[Phi]], Length[\[Phi]], Length[\[Phi]]}];
fourthderivatives = ConstantArray[0, {Length[\[Phi]], Length[\[Phi]], Length[\[Phi]], Length[\[Phi]]}];
tempfirst = ConstantArray[0, Length[\[Phi]]];
tempsecond = ConstantArray[0, {Length[\[Phi]], Length[\[Phi]]}];
tempthird = ConstantArray[0, {Length[\[Phi]], Length[\[Phi]], Length[\[Phi]]}];
tempfourth = ConstantArray[0, {Length[\[Phi]], Length[\[Phi]], Length[\[Phi]], Length[\[Phi]]}];

Do[
firstderivatives[[i]] = D[VCT[[1]], \[Phi][[i]]]/.\[Phi]vev/.rotate//Simplify;
tempfirst[[i]] =  Subscript[NW, {i}];

secondderivatives[[i]][[j]] = D[D[VCT[[1]], \[Phi][[i]]], \[Phi][[j]]]/.\[Phi]vev/.rotate//Simplify;
tempsecond[[i]][[j]] = Subscript[HW, {i, j}];

thirdderivatives[[i]][[j]][[k]] = D[D[D[VCT[[1]], \[Phi][[i]]], \[Phi][[j]]], \[Phi][[k]]]/.\[Phi]vev/.rotate//Simplify;
tempthird[[i]][[j]][[k]] = Subscript[TW, {i, j, k}];

fourthderivatives[[i]][[j]][[k]][[l]] = D[D[D[D[VCT[[1]], \[Phi][[i]]], \[Phi][[j]]], \[Phi][[k]]], \[Phi][[l]]]/.\[Phi]vev/.rotate//Simplify;
tempfourth[[i]][[j]][[k]][[l]] = Subscript[FW, {i,j,k,l}]

,{i, 1, Length@\[Phi]}, {j,1, Length@\[Phi]}, {k, 1, Length@\[Phi]}, {l, 1, Length@\[Phi]}
]


pos1 = Position[firstderivatives, 0];
pos1a = Position[tempfirst, Subscript[NW, {_}]];
realpos = Delete[pos1a, pos1];
Nterms = Extract[tempfirst, realpos];
firstderivativeseq= Extract[firstderivatives,realpos] + Nterms;


pos2 = Position[secondderivatives, 0];
pos2a = Position[tempsecond,  Subscript[HW, {_, _}]];
realpos = Complement[pos2a, pos2];
Hterms = Extract[tempsecond, realpos];
secondderivativeseq=Extract[secondderivatives, realpos] + Hterms;



pos3 = Position[thirdderivatives,0];
pos3a = Position[tempthird,  Subscript[TW, {_, _,_}]];
realpos = Complement[pos3a, pos3];
Tterms = Extract[tempthird, realpos];
thirdderivativeseq = Extract[thirdderivatives, realpos] + Tterms;


pos4 = Position[fourthderivatives, 0];
pos4a = Position[tempfourth,  Subscript[FW, {_,_,_,_}]];
realpos = Complement[pos4a, pos4];
Fterms = Extract[tempfourth, realpos];
fourthderivativeseq=Extract[fourthderivatives, realpos] + Fterms;


{constterms, coeffmatrix} = CoefficientArrays[{firstderivativeseq, secondderivativeseq, thirdderivativeseq, fourthderivativeseq}//Flatten, countertermspara]//Normal;


(* Delete duplicate row and grouping the consistent conditions *)
rownames = Array[constterms[[#]]&, Length[coeffmatrix]];
rownamesnumbers = Array[ToExpression[ToString[#]]&, Length[coeffmatrix]];
groups = GroupBy[Thread[rownames->coeffmatrix], Last->First];
coeffnames=Values[groups];
Length@coeffnames


coeffnames0 = ConstantArray[0, Length[coeffnames]];
Do[coeffnames0[[i]] = coeffnames[[i]][[1]], {i,1, Length@coeffnames}];
MatrixForm[coeffmatrix0 =DeleteDuplicates[coeffmatrix]]


(* This step cannot be done automatically, have to be treated seperately, case by case*)
(* WARNING: without careful handle the 0 row or colums inside the coeffmatrix0, the program will crash *)
SetDirectory[NotebookDirectory[]]
<<pivotgauss`
rref = RREFpivot[coeffmatrix0, -coeffnames0, 0];


(* need an imporvement to the CForm to export file in C++ for solutions *)
(* adapt the notation in class template for N and H into WeinbergNabla(i) and WeinbergHesse(i,j) *)


{mult, eqs, const} = rref//Simplify;
values = const[[1;;Length@countertermspara]];
joinlist = Join[{countertermspara}, {values}]//Transpose;
solutions = Rule@@@joinlist//MatrixForm;
coeffnames;


solutions;
consistents=ConstantArray[0, Length@const[[Length@countertermspara+1;;-1]]];
joinlist =Join[{const[[Length@countertermspara+1;;-1]]}, {consistents}]//Transpose;
consistentsa = Rule@@@joinlist//MatrixForm;


solutions/.{Subscript[HW, {a_, b_}]->HW[a-1, b-1], Subscript[NW, {a_}]->NW[a - 1], Subscript[TW, {a_, b_, c_}]->TW[a-1, b-1, c-1], Subscript[FW, {a_, b_,c_,d_}]->FW[a-1, b-1, c-1, d-1]}


testa = solutions/.{Subscript[HW, {a_, b_}]->HW[a-1,b-1], Subscript[NW, {a_}]->NW[a - 1], Subscript[TW, {a_, b_, c_}]->TW[a-1,b-1,c-1], Subscript[FW, {a_, b_,c_,d_}]->FW[a-1,b-1,c-1,d-1]};
testa = testa/.{\[Delta]ms11 -> ms11CT, \[Delta]ms22->ms22CT, \[Delta]mss->mssCT, \[Delta]Areal->ArealCT, \[Delta]Aimag->AimagCT, \[Delta]\[Lambda]1->L1CT, \[Delta]\[Lambda]2->L2CT, \[Delta]\[Lambda]3->L3CT, \[Delta]\[Lambda]4->L4CT, \[Delta]\[Lambda]5->L5CT, \[Delta]\[Lambda]6->L6CT, \[Delta]\[Lambda]7->L7CT, \[Delta]\[Lambda]8->L8CT, \[Delta]T1->T1CT, \[Delta]T2->T2CT, \[Delta]TCP->TCPCT, \[Delta]TCB->TCBCT, \[Delta]TS->TSCT}
Length@testa[[1]]


str = OpenWrite["equations.cpp"];

Do[
	WriteString[str, "\tparCT.push_back((double)("<>ToString[CForm[testa[[1]][[i]][[1]]], InputForm]<>"="<>ToString[CForm[testa[[1]][[i]][[2]]], InputForm]<>");\n"];
,{i,1,Length@testa[[1]]}]

Close[str];
Clear[str];


counta = 0;
countb = 0;
count1 = 0;
str = OpenWrite["consistents.cpp"];
Print["xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"];
Print["                                      grouped consistency conditions: "];
Print["xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"];
Do[
test = coeffnames[[j]]/.{Subscript[HW, {a_, b_}]->HW[a-1][b-1], Subscript[NW, {a_}]->NW[a - 1], Subscript[TW, {a_, b_, c_}]->TW[a-1][b-1][c-1], Subscript[FW, {a_, b_,c_,d_}]->FW[a-1][b-1][c-1][d-1]};
stra = {};
	Do[
		If[Length@test > 1 && i < Length@test, 
			stra = StringJoin[ToString[test[[i]]], " - ", ToString[test[[i+1]]], "= 0"];
			WriteString[str,  "\tparcons.push_back(abs((double)("<>ToString[test[[i]], InputForm]<>"-"<>ToString[test[[i+1]], InputForm]<>")));\n"];
			counta += 1;
			Print[stra],
			
			count1 += 1;		
		]
	,{i,1,Length@test}]
, {j,1,Length@coeffnames}]

consistents1 = consistentsa/.{Subscript[HW, {a_, b_}]->HW[a-1, b-1], Subscript[NW, {a_}]->NW[a - 1], 
				Subscript[TW, {a_, b_, c_}]->TW[a-1, b-1, c-1], 
				Subscript[FW, {a_, b_,c_,d_}]->FW[a-1, b-1, c-1, d-1]};
				
				
Print["xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"];
Print["                             Consistency conditions from redundancy degree of freedom in system of equations: "];	
Print["xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"];			

WriteString[str, " "];
Do[
stra = StringJoin[ToString[consistents1[[1]][[i]][[1]], InputForm], " = 0"];
WriteString[str,  "\tparcons.push_back(abs((double)("<>ToString[CForm[consistents1[[1]][[i]][[1]]], InputForm]<>")));\n"];
countb += 1;
Print[stra]
,{i,1,Length@consistents1[[1]]}]

Print["xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"];							
Print["                     number of consistency conditions: ", counta + countb]
Print["                     numbmer of individual not grouped consistency conditions :", count1]
Print["xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"];


Close[str];
Clear[str];	



(*solutions= solutions/.Cform/.{Subsuperscript[HesseWeinberg, a_, b_]->HesseWeinberg[a, b], Subscript[NablaWeinberg, a_]->NablaWeinberg[a]};*)



(*Do[
WriteString[str, "parCT.push_back((double)("<>ToString[solutions[[1]][[i]][[2]]//CForm]<>"));\n"]
,{i, 1, Length@solutions[[1]]}];*)




