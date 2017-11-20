(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)

(* :Title: WaveFunctions *)
(* :Context: WaveFunctions` *)
(* :Author: c8888 *)
(* :Date: 2017-11-19 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: *)
(* :Copyright: (c) 2017 c8888 *)
(* :Keywords: *)
(* :Discussion:
      * All wave functions implemented here should take an ElementaryCell2D argument and return this table with wave function values in all its points.
      * Support is a table with which we multiply the wave function to get its final space distribution. Support is constructed in "Space2D.m", not here!
      * Each wave function should have periodic boundary conditions on the boundary of the elementary cell!
*)

BeginPackage["WaveFunctions`"]
(* Exported symbols added here with SymbolName::usage *)

Needs["Space2D`"];

(* ---- Helping functions ---- *)

wannier::usage =
    "wannier[{x,y}, {xi, yi}, \[Sigma]w, \[Beta]] returns value of wannier function localized at {xj, yj}. \[Sigma]w is the Gauss peak width and \[Beta] is the normalisation factor."

wannierNormalisationFactor::usage =
    "wannierNormalisationFactor[\[Sigma]w_, \[Delta]x_, \[Delta]y_, space2DProbingPoints_] returns a normalisation factor for Wannier function localized at {0,0}"

interferenceProfileSin::usage =
    "interferenceProfileSin[k0x, k0y, kx, ky, n, m, a] returns value (positive of negative) for wave function in momentum space (definition: page 18 in lab. notebook)"

complexDotProduct::usage =
    "complexDotProduct[x, y] takes vectors x and y and returns hermitian_conjugate(y) . x"

overlapArray::usage =
    "overlapArray[array1, array2, dx, dy] calculates the squared overlap magnitude of two 2D matrices normalized with respect to spacings dx, dy between probing points"

supportHelper::usage =
    "supportHelper[elementaryCellXYTable_] returns wave function with value one in the whole elementary cell"

(* ---- Wave functions of specific models ---- *)

waveFunctionAB::usage =
    "waveFunctionAB[elementaryCellXYTable_, k0_, a_, \[Sigma]w_, hamiltonianThetaPhi___] takes elementaryCellXYTable compatible with a (that is, points with coordinates {x,y}
where -a/2 <= x <= 3/2a , -a/2 <= y <= a/2 ), state signature k0 = {k0x, k0y}, \[Sigma]w -- wannier function width and hamiltonianThetaPhi[k0] that should return
proper thetak0 and phik0 angles depending on the properties of specific system. A table with structure of elementaryCellXYTable with calculated wave function values is returned."



Begin["`Private`"]


(* ---- Helping functions ---- *)

wannierNormalisationFactor[\[Sigma]w_, \[Delta]x_, \[Delta]y_, space2DProbingPoints_] :=
    1/Sqrt[Total@Total[
    Chop@Map[wannier[#, {0, 0}, \[Sigma]w, 1] &, space2DProbingPoints, {2}]^2
  ] * \[Delta]x * \[Delta]y
]

wannier[r_,ri_,\[Sigma]w_, wannierNormalisationFactor_]:=
     wannierNormalisationFactor*Exp[(-(N@r[[1]]-N@ri[[1]])^2 - (N@r[[2]]-N@ri[[2]])^2)/2./N@\[Sigma]w^2]

interferenceProfileSin[k0x_, k0y_, kx_, ky_, n_, m_, a_] :=
    Which[Mod[a (kx - k0x), Pi] != 0,
      Sin[(2. n + 1.) a (kx - k0x)]/Sin[a (kx - k0x)]*
          Which[Mod[a (ky - k0y), Pi] != 0,
            Sin[(m + 0.5) a (ky - k0y)]/Sin[0.5 a (ky - k0y)],
            Mod[a (ky - k0y), Pi] == 0, (2 m + 1.)],
      Mod[kx, k0x] == 0, (2 n + 1.)*
        Which[Mod[a (ky - k0y), Pi] != 0,
          Sin[(m + 0.5) a (ky - k0y)]/Sin[0.5 a (ky - k0y)],
          Mod[ky, k0y] == 0, (2 m + 1.)]
    ]

complexDotProduct[x_, y_] := Chop[Dot[x, Conjugate[y]]]

overlapArray[array1_, array2_, dx_, dy_]:=
    Return[
      Total[Total[Abs[Conjugate[array1] * array2]]]^2 * dx^2 * dy^2
    ]


(* ---- Wave functions of specific models ---- *)

waveFunctionAB[elementaryCellXYTable_, k0_, a_, \[Sigma]w_, hamiltonianThetaPhi___]:=
    Module[
      {
        ret,
        \[Theta]k0,
        \[Phi]k0,
        c\[Theta]k0,
        s\[Theta]k0,
        expi\[Phi]k0,
        nodesA = {{0,0}, {0,a}, {0,-a}, {2a,a}, {2a,0}, {2a,-a}},
        nodesB = {{-a, a}, {-a, 0}, {-a, -a}, {a,a}, {a,0}, {a,-a}},
        elCellDims = Dimensions@elementaryCellXYTable
      },
      {\[Theta]k0, \[Phi]k0} = hamiltonianThetaPhi[k0];
      {c\[Theta]k02, s\[Theta]k02, expi\[Phi]k0} = {Cos[\[Theta]k0/2], Sin[\[Theta]k0/2], Exp[I \[Phi]k0]}; (* We do not want to calculate it many times. *)

      (* Site A is located at {0,0}, site B is located at {a, 0}.
         Periodicity in x is 2a, periodicity in y is a.
       *)
      (* We take into account only input from next-nearest neighbours defined in nodesA, nodesB (see p. 23 in lab notebook).
*)
      (* Let's use convention presented in Lewenstein PRL 113 045303 (2014): choice of spinors.*)
      ret = Map[
            s\[Theta]k02 *Exp[I k0.#] * Total@Function[r, Map[ wannier[r, #, \[Sigma]w, 1.]&, nodesA, {1}]][#]
            - c\[Theta]k02 * expi\[Phi]k0 * Exp[I k0.(#-{a,0})] * Total@Function[r, Map[ wannier[r, #, \[Sigma]w, 1.]&, nodesB, {1}]][#] & ,
        elementaryCellXYTable, {2}
      ];
      Return[ret] (**)
    ]

supportHelper[elementaryCellXYTable_]:=
    Map[1 &, elementaryCellXYTable, {2}]



End[] (* `Private` *)

EndPackage[]