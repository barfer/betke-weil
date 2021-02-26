#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Title:: *)
(*Analysis of the Betke-Weil Inequality*)


(* ::Text:: *)
(*by Ferenc A. Bartha*)
(*barfer@math.u-szeged.hu*)
(**)
(*Bolyai Institute, University of Szeged*)


(* ::Section:: *)
(*Map f*)


(* ::Subsection:: *)
(*Definition of the map f*)


a1 = 2;


f2 := Function[{a2, a3}, (a1 + a2 + a3) * (-a1 + a2 + a3) * (a1 - a2 + a3) * (a1 + a2 - a3)];


f1 := Function[{a2, a3, t1, t2, t3},((a1 + a2 + a3)^2 + f2[a1, a2, a3] * (t1 / a1 + t2 / a2 + t3 / a3)^2)];


f := Function[{a2, a3, t1, t2, t3}, f1[a2, a3, t1, t2, t3] - 3 * Sqrt[3] * Sqrt[f2[a2, a3]] * (1 + t1 * t2 + t2 * t3 + t3 * t1)];


(* ::Subsection:: *)
(*f(2, 2, t, t, t) = 0 for all t*)


f[2, 2, t, t, t] // Simplify


(* ::Subsection:: *)
(*Hessian of the map f*)


Hf = D[f[a2, a3, t1, t2, t3], {{a2, a3, t1, t2, t3}, 2}];


DiagX = DiagonalMatrix[{x, x, x, x, x}];
DiagX //MatrixForm


(* ::Subsubsection:: *)
(*Hf(2, 2, 0, 0, 0)*)


HfEval = Hf /. {a2 -> 2, a3 -> 2, t1 -> 0, t2 -> 0, t3 -> 0};
"Hf(2, 2, 0, 0, 0)"
HfEval //MatrixForm
"Characteristic polynomial"
CharPoly = Det[HfEval - DiagX] //Expand
"Derivative of the charactersitic polynomial"
DCharPoly = D[CharPoly, {x}] //Expand
Plot[{ CharPoly}, {x,-10.,47.}, ImageSize -> Scaled[0.5], PlotLabel -> "Characteristic polynomial"]
Plot[{DCharPoly}, {x,-10.,47.}, ImageSize->Scaled[0.5], PlotLabel -> "Derivative of the characteristic polynomial"]


(* ::Subsubsection:: *)
(*Hf(2 + r1, 2 + r2, 0 + r3, 0 + r4, 0 + r5)  for rk \in [0, 1/6] random*)


"r1, r2, r3, r4, r5"
diffs = RandomReal[{0, 1/6}, 5]
HfEval = Hf /. {a2 -> 2 + diffs[[1]], a3 -> 2 + diffs[[2]], t1 -> 0 + diffs[[3]], t2 -> 0 + diffs[[4]], t3 -> 0 + diffs[[5]]};
"Hf(2 + r1, 2 + r2, r3, r4, r5)"
HfEval //MatrixForm
"Characteristic polynomial"
CharPoly = Det[HfEval - DiagX] //Expand
"Derivative of the charactersitic polynomial"
DCharPoly = D[CharPoly, {x}] //Expand
Plot[ {CharPoly},{x,-10.,47.}, ImageSize -> Scaled[0.5], PlotLabel -> "Characteristic polynomial"]
Plot[{DCharPoly},{x,-10.,47.}, ImageSize -> Scaled[0.5], PlotLabel -> "Derivative of the characteristic polynomial"]


(* ::Section:: *)
(*Map f~*)


(* ::Subsection:: *)
(*Coordinate transform: S*)


e1 = {1, 0,         0,           0,           0};
e2 = {0, 1,         0,           0,           0};
e3 = {0, 0, 1/Sqrt[2], - 1/Sqrt[2],           0};
e4 = {0, 0, 1/Sqrt[6],   1/Sqrt[6], - 2/Sqrt[6]};
e5 = {0, 0, 1/Sqrt[3],   1/Sqrt[3],   1/Sqrt[3]};
S = Transpose[{e1, e2, e3, e4, e5}];
S // MatrixForm
S //N //MatrixForm


(* ::Subsection:: *)
(*Definition of the map f~*)


ftilde = Function[{y1, y2, y3, y4, y5}, 
	y = Transpose[{y1, y2, y3, y4, y5}];
	Sy = S . y;
	f[Sy[[1]], Sy[[2]], Sy[[3]], Sy[[4]], Sy[[5]]]
];


(* ::Subsection:: *)
(*Hessian of the map f~*)


Hftilde = D[ftilde[y1, y2, y3, y4, y5], {{y1, y2, y3, y4, y5}, 2}];


(* ::Subsubsection:: *)
(*Hf~(2, 2, 0, 0, 0)*)


HftildeEval = Hftilde /. {y1 -> 2, y2 -> 2, y3 -> 0, y4 -> 0, y5 -> 0};
"Hf~(2, 2, 0, 0, 0)"
HftildeEval //MatrixForm


(* ::Subsubsection:: *)
(*Principal minor of Hf~(2, 2, 0, 0, 0)*)


DiagXFour = DiagonalMatrix[{x, x, x, x}];
DiagXFour //MatrixForm


HftildeEvalPrincipal = HftildeEval[[1 ;; 4, 1 ;; 4]];
HftildeEvalPrincipal //MatrixForm


"Characteristic polynomial"
CharPoly = Det[HftildeEvalPrincipal - DiagXFour] //Expand
"Derivative of the charactersitic polynomial"
DCharPoly = D[CharPoly, {x}] //Expand
Plot[ {CharPoly},{x,-10.,47.}, ImageSize -> Scaled[0.5], PlotLabel -> "Characteristic polynomial"]
Plot[{DCharPoly},{x,-10.,47.}, ImageSize -> Scaled[0.5], PlotLabel -> "Derivative of the characteristic polynomial"]


(* ::Subsubsection:: *)
(*Principal minor of Hf~(2 + r1, 2 + r2, 0 + r3, 0 + r4, 0 + r5) for r1, r2 \in [0, 1/6];    r3, r4 \in [-0.14, 0.14];     r5 \in [0, 0.3]*)


"r1, r2, r3, r4, r5"
diffs = RandomReal[{0, 1/6}, 5];
diffs[[3]] = RandomReal[{-0.14, 0.14}, 1][[1]];
diffs[[4]] = RandomReal[{-0.14, 0.14}, 1][[1]];
diffs[[5]] = RandomReal[{ 0   , 0.3 }, 1][[1]];
diffs
HftildeEval = Hftilde /. {y1 -> 2 + diffs[[1]], y2 -> 2 + diffs[[2]], y3 -> 0 + diffs[[3]], y4 -> 0 + diffs[[4]], y5 -> 0 + diffs[[5]]};
"Hf~(2 + r1, 2 + r2, r3, r4, r5)"
HftildeEval //MatrixForm
"Principal minor of Hf~(2 + r1, 2 + r2, r3, r4, r5)"
HftildeEvalPrincipal = HftildeEval[[1 ;; 4, 1 ;; 4]];
HftildeEvalPrincipal //MatrixForm
"Characteristic polynomial"
CharPoly = Det[HftildeEvalPrincipal - DiagXFour] //Expand
"Derivative of the charactersitic polynomial"
DCharPoly = D[CharPoly, {x}] //Expand
Plot[ {CharPoly},{x,-2.,47.}, ImageSize -> Scaled[0.5], PlotLabel -> "Characteristic polynomial"]
Plot[{DCharPoly},{x,-2.,47.}, ImageSize -> Scaled[0.5], PlotLabel -> "Derivative of the characteristic polynomial"]



