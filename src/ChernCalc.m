(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)

(* :Title: ChernCalc *)
(* :Context: ChernCalc` *)
(* :Author: c8888 *)
(* :Date: 2017-11-23 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: *)
(* :Copyright: (c) 2017 c8888 *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["ChernCalc`"]
(* Exported symbols added here with SymbolName::usage *)

Begin["`Private`"]

link[ckBZ_, i1_, j1_, i2_, j2_] :=
    Exp[I Arg@complexDotProduct[ckBZ[[i2, j2]], ckBZ[[i1, j1]] ] ]

FxyT[ckBZ_] := Table[Log[
  link[ckBZ, i, j, i + 1, j]*
      link[ckBZ, i + 1, j, i + 1, j + 1]*
      link[ckBZ, i + 1, j + 1, i, j + 1]*
      link[ckBZ, i, j + 1, i, j]
], {i, 1, Dimensions[ckBZ][[1]] - 1}, {j, 1, Dimensions[ckBZ][[2]] - 1}]


End[] (* `Private` *)

EndPackage[]