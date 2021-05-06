(* ::Package:: *)

(* ::Input:: *)
(**)


(* improved from non-pivoting Gauss-elimination from *)
(* https://mathematica.stackexchange.com/questions/139606/find-elementary-matrices-that-produce-rref?rq=1*)
(* to pivoting Gauss-elimination *)
(* written by N. N. Tien Dat, August, 8th, 2020 *)


Format[matWithDiv[n_,opts:OptionsPattern[Grid]][m_?MatrixQ]]:=MatrixForm[{{Grid[m,opts,Dividers->{n->{Red,Dashed}}]}}];

switchrows[mat_, {i_Integer, j_Integer}]:= Module[{lmat = mat, len=Length[mat]}, 
	lmat[[{i, j}]] = lmat[[{j, i}]];
	lmat
];  


RREFpivot[mat_?(MatrixQ[#] &), b_?(VectorQ[#] &), verbose_]:= Module[{i, j, multiplier, pivot, augmat, npara = Length@mat[[1]], nsol = Length@mat, inverse, eqs, sol},
	If[Length@b != Length@mat, Return["Size of b vector not the same as number of row in A matrix"]];

	augmat = Transpose[Join[Transpose[mat], {b}]];
	Print["rank of (mat, augmat) is (", MatrixRank[mat],",", MatrixRank[augmat],")"];
	Print[">>>>>>>>>>>>>>>>Starting forward Gaussian elimination phase using", augmat[[All, 1;; npara+ 1]]//matWithDiv[npara+1, Background->LightOrange],">>>>>>>>>>>>>>>>>"];
	(*augmat = ArrayFlatten[{{augmat, IdentityMatrix[nsol]}}]*);
	augmat = Transpose[Join[Transpose[augmat], IdentityMatrix[nsol]]];
	Do[
	    If[verbose === 0, Continue,
		Print["pivot now is (",pivot,",",pivot ")"]];
		Do[
			If[augmat[[pivot, pivot]] == 0,
				(* find new row to swap *)
				posrow = Position[augmat[[;;, pivot]], 0];
				posrow = List[Complement[Range@Length@augmat[[;;, pivot]], posrow//Flatten]]//Transpose;
				posrow  = List[Select[posrow//Flatten, #>pivot &]]//Transpose;
				If[Length[posrow] != 0,
					posrow = posrow[[1]][[1]];
					If[verbose === 0, Continue,
					Print["element (",pivot,",", pivot,") is (",augmat[[pivot, pivot]],") thus swaping row ", pivot," to row ", posrow]];
					(* swap pivot to posrow *)
					augmat = switchrows[augmat, {posrow, pivot}];
					If[verbose == 0,
						Continue,
						Print[augmat[[All, 1;; npara + 1]]//matWithDiv[npara+1, Background->LightOrange], MatrixForm[augmat[[All, npara+2 ;;]]]]
					];
					
					If[verbose === 0, Continue,
					Print["elements below current pivot (",pivot,",", pivot,") are zero, move to the next pivot"]];
					Goto[begin];
				], augmat
			];
		
			multiplier = augmat[[j, pivot]]/augmat[[pivot, pivot]];
			If[verbose === 0, Continue,
			Print["will now zero out element (",j,",",pivot,") by subtracting ",multiplier," times row ",pivot," from row ",j]];
			augmat[[j, pivot ;;]] = augmat[[j, pivot ;;]] - multiplier*augmat[[pivot, pivot;;]];
			Label[begin];
			If[verbose == 1,
				Print[augmat[[All, 1;; npara + 1]]//matWithDiv[npara+1, Background->LightOrange], MatrixForm[augmat[[All, npara+2 ;;]]]],
				Continue
				]
				, {j, pivot + 1, nsol}
			]
			, {pivot, 1, npara}
	];
	If[verbose == 2, Print[augmat[[All, 1;;npara+1]]//matWithDiv[npara+1, Background->LightOrange], MatrixForm[augmat[[All, npara +2;;]]]], Continue];
	
	Print[">>>>>>>>>>>>>>>>>>>Starting backward elimination phase>>>>>>>>>>>>>>>>>>>>>>"];
	Do[
		Do[
			multiplier = augmat[[j, pivot]]/augmat[[pivot, pivot]];
			If[verbose === 0, Continue,	
			Print["will now zero out element (",j,"," ,pivot,") by subtracting ",  multiplier, " times row ", pivot," from row ", j]];
			augmat[[j, pivot;;]] = augmat[[j, pivot ;;]] - multiplier*augmat[[pivot, pivot;;]];
			If[verbose == 1, 
				Print[augmat[[All, 1;;npara+1]]//matWithDiv[npara+1, Background->LightOrange], MatrixForm[augmat[[All, npara +2;;]]]],
				Continue
				]
			,{j,1, pivot-1}
		],
		{pivot, 2, npara}
	];
	If[verbose == 2, Print[augmat[[All, 1;;npara+1]]//matWithDiv[npara+1, Background->LightOrange], MatrixForm[augmat[[All, npara +2;;]]]], Continue];

	Print[">>>>>>>>>>>>>>>>>>>Starting Final phase, convert reduced echelon to identity matrix>>>>>>>>>>>>>>>"]; Do[augmat[[j,;;]]=augmat[[j,;;]]/augmat[[j,j]],{j,1,npara}];
	Print[augmat[[All,1;;npara+1]]//matWithDiv[npara+1,Background->LightOrange],MatrixForm[augmat[[All,npara+2;;]]]];
	
	Print["Inverse Matrix is ",MatrixForm[augmat[[All,npara+2;;]]]];
	Print["Solution  vector is ",MatrixForm[augmat[[All,npara+1]]]];
	inverse = augmat[[All, npara+2;;]];
	sol = augmat[[All, npara+1]];
	eqs = augmat[[;;, ;;npara]];
	Return[{inverse, eqs, sol}]
]


(* test *)
(*ClearAll[p, q, r, t]
mat={{5,2,18},{0,0,2},{4,0,12},{2,3,8}, {1, -1 , 0}};
X = {x, y, z};
b={p,q,r,t, w};

rref = RREFpivot[mat, b, 0];*)


(*(* solving system of equations with reduce *)
eqs = mat.X
(* solving equations *)
Reduce[eqs[[1]]==b[[1]] && eqs[[2]]==b[[2]] && eqs[[3]] == b[[3]] && eqs[[4]] == b[[4]] && eqs[[5]] == b[[5]], X, Reals]
(* solving constraints and solution from RREFpivot *)
ClearAll[p, q, r, t]
constrain = {rref[[3]][[4]], rref[[3]][[5]]};
sol = {rref[[3]][[1]], rref[[3]][[2]], rref[[3]][[3]]}
consistents = Reduce[constrain[[1]] ==0 && constrain[[2]] == 0, {q,p}, Reals]
sol/.ToRules[consistents]//Simplify*)


(* ::Input:: *)
(**)
