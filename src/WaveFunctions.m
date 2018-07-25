(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)

(* :Title: WaveFunctions *)
(* :Context: WaveFunctions` *)
(* :Author: c8888 *)
(* :Date: 2017-11-19 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 11.1*)
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
    "wannierNormalisationFactor[\[Sigma]w_, \[Delta]x_, \[Delta]y_, rangeRectangleSizeX_, rangeRectangleSizeY_] returns a normalisation factor for Wannier function localized at {0,0}"

interferenceProfileSin::usage =
    "interferenceProfileSin[k0_, k_, numCellsX_, numCellsY_, a_] returns value (positive or negative) for wave
    function in momentum space (definition: page 18 in lab. notebook)"

complexDotProduct::usage =
    "complexDotProduct[x, y] takes vectors x and y and returns hermitian_conjugate(y) . x"

overlapArray::usage =
    "overlapArray[array1, array2, \[Delta]x_, \[Delta]y_] calculates the squared overlap magnitude of two 2D matrices. We asssume
    dx and dy do not change anywhere during execution of the program!"

supportHelper::usage =
    "supportHelper[elementaryCellXYTable_] returns wave function with value one in the whole elementary cell"

wannierRectangleTable::usage =
    "wannierRectangleTable[fullSpace2DNodes_, nodesNeighbourhood_, wannierNormalisationFactorValue_, \[Sigma]w_, \[Delta]x_, \[Delta]y_]
    returns table of wannier function values in the neighbourhood of every node to significantly speed up
    calculations of overlap."

blochNodePhaseConjugateTable::usage =
    "blochNodePhaseTable[k0x_, k0y_, nodesXYTable_] :=
    Exp[I {k0x, k0y}.#] & /@ nodesXYTable returns values of conjugate bloch phases at nodes positions"

wannierProject::usage =
    "wannierProject[array2D_, nodesNeighbourhood_, wannierRectangleTableValues_, \[Delta]x_, \[Delta]y_] returns a wave function
    vector in wannier function basis"

kProject::usage =
    "kProject[ckWannier_, blochNodePhaseConjugateTable_, q_, ckSupportMemberTable_] takes table from wannierProject and cancels out the phase
    resulting from bloch wave."

overlapWannier::usage =
    "overlapWannier[ck1_, ck2_] returns Abs[complexDotProduct[ck1, ck2]]^2"

wannierFTProfileTable::usage =
    "wannierFTProfileTableSq[\[Delta]kx_, \[Delta]ky_, \[Sigma]w_, dimx_, dimy_] returns squared table of
    wannier envelope"

interferenceProfileSin00Table::usage =
    "interferenceProfileSin00Table[numCellsX_, numCellsY_, \[Delta]kx_, \[Delta]ky_, ax_, ay_, \[Sigma]w_, dimx_, dimy_, dimBZx_, dimBZy_ ] returns a 2D table with dimenstions (dimx + dimBZx)x(dimy + dimBZy) that contains interference profile from a single fermion located at k = (0,0) (corresponds to center of the table), with highest peak with height 1. Further usage: see fermionIntProfile function."

fermionIntProfile::usage =
    "fermionIntProfile[nx_, ny_, intProfSin00Table_, gx_, gy_, dimx_, dimy_, dimBZx_, dimBZy_] returns intProfSin00Table[[-nx gx + dimBZx/2 + 1 ;; -nx gx + dimBZx/2 +
   dimx, -ny gy + dimBZy/2 + 1 ;; -ny gy + dimBZy/2 + dimy]] which is an interference profile of fermion located at
   (nx \[Delta]k0x, ny \[Delta]k0y)."

findHeight::usage =
    "findHeight[nx_, ny_, allFermionsAbsSq_, gx_, gy_, dimx_, dimy_] returns height of peak coming from a fermion with
    indices (nx, ny)."

hamiltonianHarperQ::usage =
    "hamiltonianHarperQ[k0x_, k0y_, J_, J1_, q_] returns a q*q matrix of Harper hamiltonian"

myES::usage =
    "myES[hamiltonian_] returns an eigensystem of the hamiltonian matrix sorted by energy."

addNoise::usage =
    "addNoise[wfKSpaceAbsSq_, SNR_, FBZnColnRow_] returns noised square of magnitude. signal is mean over FBZ."

hamiltonianHaldane::usage =
    "hamiltonianHaldane[kx_, ky_, J_,
  J1_, \[Theta]_, \[CapitalDelta]_] is the Haldane hamiltonian on a brick wall"


stateVectorHaldaneABAB::usage =
    "stateVectorHaldaneABAB[stateVectorHaldaneAB_] :=
    Normalize@Flatten[{stateVectorHaldaneAB, stateVectorHaldaneAB}] takes AB vector and makes normalized ABAB vector"

stateVectorHaldaneAB::usage =
    "stateVectorHaldaneAB[stateVectorHaldaneABAB_] takes ABAB vector and returns normalized AB vector"

(* ---- Wave functions of specific models ---- *)

blochWave::usage =
    "blochWave[fullSpaceXYTable_, k0_] returns bloch wave values in the whole space. Wave function has to be multiplied by it."

waveFunctionAB::usage =
    "waveFunctionAB[elementaryCellXYTable_, k0_, a_, \[Sigma]w_, hamiltonianThetaPhi___] takes elementaryCellXYTable compatible with a (that is, points with coordinates {x,y}
where -a/2 <= x <= 3/2a , -a/2 <= y <= a/2 ), state signature k0 = {k0x, k0y}, \[Sigma]w -- wannier function width and hamiltonianThetaPhi[k0] that should return
proper thetak0 and phik0 angles depending on the properties of specific system. A table with structure of elementaryCellXYTable with calculated wave function values is returned."

stateVectorABGround::usage =
    "stateVectorABGround[k0_, a_, J_, J1_] returns a normalized vector {uA, uB} for a ground state in harper model with 2 sites"

stateVectorABExcited::usage =
    "stateVectorABExcited[k0_, a_, J_, J1_] returns a normalized vector {uA, uB} for a ground state in harper model with 2 sites"

waveFunctionABkSpace::usage =
    "waveFunctionABkSpace[fullSpaceKxKyTable_, k0_, a_, \[Sigma]w_, numCellsX_, numCellsY_, J_, J1_, n_] returns one particle wave function in K-space"

waveFunctionABkSpaceFermionsAbsSquared::usage =
    "waveFunctionABkSpaceFermionsAbsSquared[fullSpaceKxKyTable_, k0Table_, a_, \[Sigma]w_, numCellsX_, numCellsY_, J_, J1_, n_] returns squared modulus of wave function in K-space for many fermions"

waveFunctionQ::usage =
    "waveFunctionQ[elementaryCellXYTable_, k0x_, k0y_, \[Sigma]w_,
  stateVectorQ_, nodesExactPositions_] returns table of wave function values in elementary cell"

fastFullSpaceWfQRSpace::usage =
    "fastFullSpaceWfQRSpace[k0x_, k0y_, stateVectorQ_, \[Sigma]w_, numCellsX_,
  numCellsY_, nodesExactPositions_, elementaryCellXYTable_,
  fullSpaceXYTable_, \[Delta]x_, \[Delta]y_, dimx_, dimy_, support_, ax_, ay_] return full space table generated from
  bloch wave
  and wave
   function in elementary cell."
fastFullSpaceWfQRSpaceAbs::usage =
    "fastFullSpaceWfQRSpaceAbs[stateVectorQ_, \[Sigma]w_, numCellsX_,
  numCellsY_, nodesExactPositions_, elementaryCellXYTable_,
  fullSpaceXYTable_, \[Delta]x_, \[Delta]y_, dimx_, dimy_, support_, ax_, ay_] works identical to
  fastFullSpaceWfQRSpace except it
  returns only the absolute value. Because of that it is faster since no bloch wave has to be constructed. Gives wf
  not normalized properly if the support is not rectangular."

Begin["`Private`"]


(* ---- Helping functions ---- *)

wannierNormalisationFactor[\[Sigma]w_, \[Delta]x_, \[Delta]y_, rangeRectangleSizeX_, rangeRectangleSizeY_] :=
    Module[
      {
        nx = Round[rangeRectangleSizeX / 2 / \[Delta]x],
        ny = Round[rangeRectangleSizeX / 2 / \[Delta]y]
      },
      Return@(1 / Sqrt[Total@Total[
        Chop@Map[wannier[#, {0, 0}, \[Sigma]w, 1] &,
          Table[{x, y}, {x, -nx * \[Delta]x, rangeRectangleSizeX / 2, \[Delta]x}, {y, -ny * \[Delta]y, rangeRectangleSizeY / 2, \[Delta]y}],
          {2}]^2
      (*We explicitly make sure that we cross the zero point during normalisation!!!! We translate the rectangular net if necessary. *)
      ] * \[Delta]x * \[Delta]y
      ])
    ]

wannier[r_, ri_, \[Sigma]w_, wannierNormalisationFactor_] :=
    wannierNormalisationFactor * Exp[(-(N@r[[1]] - N@ri[[1]])^2 - (N@r[[2]] - N@ri[[2]])^2) / 2. / N@\[Sigma]w^2]

wannierFTProfileTable[\[Delta]kx_, \[Delta]ky_, \[Sigma]w_, dimx_, dimy_]:=
    Table[Exp[-((i * \[Delta]kx)^2 + (j * \[Delta]ky)^2) * \[Sigma]w^2/2], {i, -dimx/2, dimx/2-1},{j, -dimy/2,
      dimy/2-1}];

interferenceProfileSin00Table[numCellsX_, numCellsY_, \[Delta]kx_, \[Delta]ky_, ax_, ay_, \[Sigma]w_, dimx_, dimy_, dimBZx_, dimBZy_ ] :=
    Module[{ret = {}},
      ret = Table[(*Exp[-((i * \[Delta]kx)^2 + (j * \[Delta]ky)^2) / 2 * \[Sigma]w^2]**)
          If[Abs[Mod[1 / 2 ax i * \[Delta]kx, Pi]] < 16 $MachineEpsilon ||
              Abs[Mod[1 / 2 ax i * \[Delta]kx, Pi] - Pi] < 16 $MachineEpsilon ,
            1, Sin[(numCellsX + 1 / 2) ax i \[Delta]kx] /
              Sin[ax * 1 / 2 * (i \[Delta]kx)] / (2 numCellsX + 1)] *
          If[Abs[Mod[1 / 2 ay j * \[Delta]ky, Pi]] < 16 $MachineEpsilon ||
              Abs[Mod[1 / 2 ay j * \[Delta]ky, Pi] - Pi] < 16 $MachineEpsilon,
            1, Sin[(numCellsY + 1 / 2) ay j \[Delta]ky] /
              Sin[ay * 1 / 2 * (j \[Delta]ky)] / (2 numCellsY + 1)] ,
        {i, -(dimx + dimBZx) / 2, (dimx + dimBZx) / 2 - 1, 1},
        {j, -(dimy + dimBZy) / 2, (dimy + dimBZy) / 2 - 1, 1}];
      Return[ret];
    ]


fermionIntProfile[nx_, ny_, intProfSin00Table_, gx_, gy_, dimx_, dimy_, dimBZx_, dimBZy_] :=
    Return[
      intProfSin00Table[[-nx gx + dimBZx / 2 + 1 ;; -nx gx + dimBZx / 2 +
          dimx, -ny gy + dimBZy / 2 + 1 ;; -ny gy + dimBZy / 2 + dimy]]
    ]

findHeight[nx_, ny_, allFermionsAbsSq_, gx_, gy_, dimx_, dimy_]:=
    Return[allFermionsAbsSq[[dimx/2 + 1 + nx gx, dimy/2 + 1 + ny gy]]]

complexDotProduct[x_, y_] := Chop[Dot[x, Conjugate[y]]]

overlapArray[array1_, array2_, \[Delta]x_, \[Delta]y_] :=
    Return[
      Abs[Total[Total[Conjugate[array1] * array2]] \[Delta]x * \[Delta]y]^2
    ]

supportHelper[elementaryCellXYTable_] :=
    Map[1 &, elementaryCellXYTable, {2}]

wannierRectangleTable[fullSpace2DNodes_, nodesNeighbourhood_, wannierNormalisationFactorValue_, \[Sigma]w_, \[Delta]x_, \[Delta]y_] :=
    Module[ (* All data calculated using only first node *)
      {
        rectangleDims = First@nodesNeighbourhood,
        firstNode = First@fullSpace2DNodes
      },
      Return@Map[wannierNormalisationFactorValue * (* See page 25 in lab notebook *)
          Exp[-Total[((firstNode[[1]] - #) * {\[Delta]x, \[Delta]y} - firstNode[[2]])^2] / 2. / N@\[Sigma]w^2] &,
        Table[{i, j}, {i, rectangleDims[[1, 1]], rectangleDims[[2, 1]]}, {j, rectangleDims[[1, 2]], rectangleDims[[2, 2]]}], {2}]
    ]

wannierProject[array2D_, nodesNeighbourhood_, wannierRectangleTableValues_, \[Delta]x_, \[Delta]y_] :=
    Return@(Map[Total@Total[wannierRectangleTableValues * array2D[[ #[[1, 1]] ;; #[[2, 1]] , #[[1, 2]] ;; #[[2, 2]] ]]]&,
      nodesNeighbourhood, {1}] * \[Delta]x * \[Delta]y)

blochNodePhaseConjugateTable[k0x_, k0y_, nodesXYTable_] :=
    Exp[I {k0x, k0y}.#] & /@ nodesXYTable

kProject[ckWannier_, blochNodePhaseConjugateTable_, q_, ckSupportMemberTable_] := (*requires so that the nodes are in
    good order in the nodes
    list.*)
    Normalize@Table[
      Total[(ckWannier*blochNodePhaseConjugateTable*ckSupportMemberTable)[[
          1 + i*Length[ckWannier]/q ;; (i + 1)*Length[ckWannier]/q]]]/
          Total[ckSupportMemberTable[[
              1 + i*Length[ckWannier]/q ;; (i + 1)*Length[ckWannier]/q]]], {i, 0,
        q - 1}]

overlapWannier[ck1_, ck2_] :=
    Abs[complexDotProduct[ck1, ck2]]^2

Ay[m_, q_] := 2 Pi m * 1/q

hamiltonianHarperQ[k0x_, k0y_, J_, J1_, q_] := SparseArray[{
  {i_, i_} -> -2. J1 Cos[k0y - Ay[i, q]], (*diag*)
  {i_, j_} /;
      i - j == 1 -> -J Exp[-I k0x], (*pozadiag*)
  {i_, j_} /;
      i - j == -1 -> -J Exp[I k0x]},(*pozadiag*)
  {q, q}] +
    SparseArray[{
      {i_, j_} /; {i, j} == {1, q} -> -J Exp[-I k0x],(*róg p-
    górny*)
      {i_, j_} /; {i, j} == {q, 1} -> -J Exp[
        I k0x]  (*róg l-dolny*)
    },
      {q, q}]

myES[hamiltonian_] :=
    Module[{es =
        Transpose@
            Sort[Transpose[Eigensystem[hamiltonian]], #1[[1]] < #2[[1]] &]},
      es[[2]] = Map[Exp[-I Arg[#[[1]]]]*# &, es[[2]]]; Return[es] (*es[[2,
  i,1]] is always real*)]

addNoise[wfKSpaceAbsSq_, SNR_, FBZnColnRow_] := Module[
  {
    measuredSq = mirror2DSpace[wfKSpaceAbsSq],
    A, \[Sigma]n
  },
  A = Mean[
    Flatten[measuredSq[[FBZnColnRow[[1, 1]] ;; FBZnColnRow[[1, 2]],
        FBZnColnRow[[2, 1]] ;;
            FBZnColnRow[[2, 2]] ]]]];(*signal strength*)
  \[Sigma]n = A/SNR;
  Return[wfKSpaceAbsSq +
      RandomReal[{0, \[Sigma]n}, {Dimensions[wfKSpaceAbsSq][[1]],
        Dimensions[wfKSpaceAbsSq][[2]]}]]
]

hamiltonianHaldane[kx_, ky_, J_,
  J1_, \[Theta]_, \[CapitalDelta]_] := {
  {
    \[CapitalDelta] - 2 J1 (Cos[\[Theta] + 2 kx] + Cos[\[Theta] - kx - ky] +
        Cos[\[Theta] - kx + ky]),
    -J (2 Cos[kx] + Exp[-I ky])
  }, {
    -J (2 Cos[kx] + Exp[I ky]),
    -\[CapitalDelta] - 2 J1 (Cos[\[Theta] - 2 kx] + Cos[\[Theta] + kx + ky] +
        Cos[\[Theta] + kx - ky])
  }
}

stateVectorHaldaneABAB[stateVectorHaldaneAB_] :=
    Normalize@Flatten[{stateVectorHaldaneAB, stateVectorHaldaneAB}]

stateVectorHaldaneAB[stateVectorHaldaneABAB_] :=
    Normalize[{Mean[stateVectorHaldaneABAB[[{1, 3}]]],
      Mean[stateVectorHaldaneABAB[[{2, 4}]]]}]

(* ---- Wave functions of specific models ---- *)

blochWave[fullSpaceXYTable_, k0_] :=
    Return@Map[Exp[I k0.#]&, fullSpaceXYTable, {2}]

stateVectorABGround[k0_, a_, J_, J1_] :=
    Return@Normalize[{(-4 * J * Cos[k0[[2]]] - Sqrt[2] * Sqrt[9 * J^2 + 6 * J * J1 + 5 * J1^2 + (3 * J^2 + 10 * J * J1
        + 3 * J1^2) * Cos[2 * k0[[1]]] + 4 * J^2 * Cos[2 * k0[[2]]]]), (-2 * ((2 * J + 2 * J1) * Cos[k0[[1]]] + I(-J + J1) * Sin[k0[[1]]]))}]

stateVectorABExcited[k0_, a_, J_, J1_] :=
    Return@Normalize[{(-4 * J * Cos[k0[[2]]] + Sqrt[2] * Sqrt[9 * J^2 + 6 * J * J1 + 5 * J1^2 + (3 * J^2 + 10 * J * J1
        + 3 * J1^2) * Cos[2 * k0[[1]]] + 4 * J^2 * Cos[2 * k0[[2]]]]), (-2 * ((2 * J + 2 * J1) * Cos[k0[[1]]] + I(-J + J1) * Sin[k0[[1]]]))}]

waveFunctionAB[elementaryCellXYTable_, k0_, a_, \[Sigma]w_, J_, J1_, n_] :=
    Module[
      {
        ret,
        uA, uB,
      (*nodesA = {{0,0}, {0,a}, {0,-a}, {2a,a}, {2a,0}, {2a,-a}},
        nodesB = {{-a, a}, {-a, 0}, {-a, -a}, {a,a}, {a,0}, {a,-a}},*)
        nodesA = {{-a / 2, a}, {-a / 2, 0}, {-a / 2, -a}, {3a / 2, a}, {3a / 2, 0}, {3a / 2, -a}},
        nodesB = {{-3a / 2, a}, {-3a / 2, 0}, {-3a / 2, -a}, {a / 2, a}, {a / 2, 0}, {a / 2, -a}},
        elCellDims = Dimensions@elementaryCellXYTable
      },
      {uA, uB} = Which[n == 0, stateVectorABGround[k0, a, J, J1], n == 1, stateVectorABExcited[k0, a, J, J1] ];

      (* Site A is located at {-a/2,0}, site B is located at {a/2, 0}.
         Periodicity in x is 2a, periodicity in y is a.
       *)
      (* We take into account only the input from next-nearest neighbours defined in nodesA, nodesB (see p. 23 in lab notebook).
*)
      (* Let's use convention presented in Lewenstein PRL 113 045303 (2014): choice of spinors.*)
      ret = Map[
        uA (*Exp[I k0.#]*) * Total@Function[r, Map[ wannier[r, #, \[Sigma]w, 1.]&, nodesA, {1}]][#]
            + uB (*Exp[I k0.#]*) (** Exp[-I k0.{a,0}] *) * Total@Function[r, Map[ wannier[r, #, \[Sigma]w, 1.]&, nodesB, {1}]][#] & ,
        elementaryCellXYTable, {2}
      ];
      Return[ret] (*not normalized yet, no bloch wave*)
    ]

waveFunctionABkSpace[fullSpaceKxKyTable_, k0_, a_, \[Sigma]w_, numCellsX_, numCellsY_, J_, J1_, n_] :=
    Module[
      {
        ret,
        uA, uB
      },
      {uA, uB} = Which[n == 0, stateVectorABGround[k0, a, J, J1], n == 1, stateVectorABExcited[k0, a, J, J1] ];
      ret = Map[
        wannier[#, k0, 1 / \[Sigma]w, 1] * interferenceProfileSin[k0, #,
          numCellsX, numCellsY, a] * (uA + uB * Exp[I(#[[1]] - k0[[1]])a]) & ,
        fullSpaceKxKyTable, {2}
      ];
      Return[ret] (*not normalized yet*)
    ]

waveFunctionABkSpaceFermionsAbsSquared[fullSpaceKxKyTable_, k0Table_, a_, \[Sigma]w_, numCellsX_, numCellsY_, J_, J1_, n_] :=
    Module[
      {
        ret,
        uATable, uBTable
      },
      uATable = Which[n == 0, Map[stateVectorABGround[#, a, J, J1][[1]]&, k0Table, {1}], n == 1, Map[stateVectorABExcited[#, a, J, J1][[1]]&, k0Table, {1}] ];
      uBTable = Which[n == 0, Map[stateVectorABGround[#, a, J, J1][[2]]&, k0Table, {1}], n == 1, Map[stateVectorABExcited[#, a, J, J1][[2]]&, k0Table, {1}] ];
      ret = Map[
        Total@Function[k, Map[wannier[k, #, 1.4142135623730950488 / \[Sigma]w, 1] * interferenceProfileSin[#, k,
          numCellsX, numCellsY, a]^2&, k0Table, {1}] * Abs[(uATable + uBTable * Map[Exp[I(k[[1]] - #[[1]])a]&, k0Table, {1}])]^2][#] & ,
        fullSpaceKxKyTable, {2}
      ];
      Return[ret] (*not normalized yet*)
    ]

waveFunctionQ[elementaryCellXYTable_, k0x_, k0y_, \[Sigma]w_,
  stateVectorQ_, nodesExactPositions_] :=
    Map[stateVectorQ.Function[r,
      Map[Exp[-I (k0x #[[1]] + k0y #[[2]])]*
          wannier[r, #, \[Sigma]w, 1] &,
        nodesExactPositions, {1}]][#] &, elementaryCellXYTable, {2}]

fastFullSpaceWfQRSpace[k0x_, k0y_, stateVectorQ_, \[Sigma]w_, numCellsX_,
  numCellsY_, nodesExactPositions_, elementaryCellXYTable_,
  fullSpaceXYTable_, \[Delta]x_, \[Delta]y_, dimx_, dimy_, support_, ax_, ay_] :=
    Module[{
      wfQRSpaceElCell,
      wfQRSpaceFullSpace,
      norm
    },
      wfQRSpaceElCell =
          waveFunctionQ[elementaryCellXYTable, k0x, k0y, \[Sigma]w,
            stateVectorQ, nodesExactPositions];
      norm = Sqrt[Total@Total[Abs[wfQRSpaceElCell]^2] \[Delta]x \[Delta]y];
      Return[
          support*CenterArray[
            1/(norm*Sqrt[(2 numCellsX + 1) (2 numCellsY + 1)]) (*blochWave[
              fullSpaceXYTable, {k0x, k0y}]*)connectBlochPhase2DElementaryCells[k0x, k0y, numCellsX,
              numCellsY, elementaryCellXYTable, ax, ay]*
                connect2DElementaryCells[numCellsX, numCellsY,
                  wfQRSpaceElCell], {dimx, dimy}]
      ]
    ]

fastFullSpaceWfQRSpaceAbs[stateVectorQ_, \[Sigma]w_, numCellsX_,
  numCellsY_, nodesExactPositions_, elementaryCellXYTable_,
  fullSpaceXYTable_, \[Delta]x_, \[Delta]y_, dimx_, dimy_, support_, ax_, ay_] :=
    Module[{
      wfQRSpaceElCell,
      norm
    },
      wfQRSpaceElCell =
          Abs@waveFunctionQ[elementaryCellXYTable, 0, 0, \[Sigma]w,
            stateVectorQ, nodesExactPositions];
      norm = Sqrt[Total@Total[Abs[wfQRSpaceElCell]^2] \[Delta]x \[Delta]y];
      Return[
        support*CenterArray[
          1/(norm*Sqrt[(2 numCellsX + 1) (2 numCellsY + 1)])*connect2DElementaryCells[numCellsX, numCellsY,
                wfQRSpaceElCell], {dimx, dimy}]
      ]
    ]


End[] (* `Private` *)

EndPackage[]