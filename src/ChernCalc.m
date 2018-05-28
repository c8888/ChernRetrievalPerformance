(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)

(* :Title: ChernCalc *)
(* :Context: ChernCalc` *)
(* :Author: c8888 *)
(* :Date: 2017-11-23 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 11.1 *)
(* :Copyright: (c) 2017 c8888 *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["ChernCalc`"]
(* Exported symbols added here with SymbolName::usage *)

Needs["Space2D`"];
Needs["WaveFunctions`"];
Needs["HIOER`"];


FxyT::usage =
    "returns a table of Berry curvature in BZ (except egdes)"

findCkRetrSupportQ::usage =
    "findCkRetrSupportQ[wfQRSpaceFullSpace_, ckModel_, nodesNeighbourhoods_,
  wannierRectangleTableValues_, \[Delta]x_, \[Delta]y_, support_, nIterations_, nRepeats_, nHIO_, gamma_, nSets_] returns a table of model wave function in wannier basis
  for given Harper model parameters"
findCkRetrSupportAbsImposeQ::usage =
    "findCkRetrSupportAbsImposeQ[wfQRSpaceFullSpace_, ckModel_, nodesNeighbourhoods_,
  wannierRectangleTableValues_, \[Delta]x_, \[Delta]y_, support_, nIterations_, nRepeats_, nHIO_, gamma_, nSets_,
  nEREnd_, nAbsImpose_,
  nAbsImposeStart_, nAbsImposeEnd_, q_, \[Sigma]w_,
  numCellsX_, numCellsY_, nodesExactPositions_, elementaryCellXYTable_,
  fullSpaceXYTable_, ckSupportMemberTable_,  SNR_, FBZnColnRow_, \[Sigma]resolution_] returns a table of model wave function in wannier
  basis
  for given Harper model parameters"

Begin["`Private`"]

ComplexDotProduct[x_, y_] := Chop[Dot[x, Conjugate[y]]]

link[ckBZ_, i1_, j1_, i2_, j2_] :=
    Exp[I Arg@complexDotProduct[ckBZ[[i2, j2]], ckBZ[[i1, j1]]]]

FxyT[ckBZ_] :=
    Table[Log[
      link[ckBZ, i, j, i + 1, j]*link[ckBZ, i + 1, j, i + 1, j + 1]*
          link[ckBZ, i + 1, j + 1, i, j + 1]*
          link[ckBZ, i, j + 1, i, j]], {i, 1,
      Dimensions[ckBZ][[1]] - 1}, {j, 1, Dimensions[ckBZ][[2]] - 1}]

findCkRetrSupportQ[wfQRSpaceFullSpace_, ckModel_, nodesNeighbourhoods_,
  wannierRectangleTableValues_, \[Delta]x_, \[Delta]y_, support_, nIterations_, nRepeats_, nHIO_, gamma_, nSets_] :=
      Module[{
        ckRetr,
        ckRetrMirror,
        retr,
        overlapRetr,
        overlapRetrMirror,
        distKSpace
      },
        Return@Table[
            retr = phaseRetrieveSupport[Abs@Fourier[wfQRSpaceFullSpace], support, nIterations,
              nRepeats, nHIO, gamma];
            distKSpace = 1/Total[Total[Abs[Fourier[wfQRSpaceFullSpace]]^2]]*Total@Total[Abs[Abs[Fourier
            [wfQRSpaceFullSpace]]-Abs[Fourier[retr]]]];

            ckRetr = wannierProject[retr, nodesNeighbourhoods, wannierRectangleTableValues, \[Delta]x, \[Delta]y];
            overlapRetr = Abs[ckRetr.Conjugate[ckModel]]^2/Abs[ckModel.Conjugate[ckModel]]^2;

            ckRetrMirror = wannierProject[mirrorXY[retr], nodesNeighbourhoods, wannierRectangleTableValues, \[Delta]x,
              \[Delta]y];
            overlapRetrMirror = Abs[ckRetrMirror.Conjugate[ckModel]]^2;

            If[overlapRetr > overlapRetrMirror,
              {ckRetr, overlapRetr, distKSpace},
              {ckRetrMirror, overlapRetrMirror, distKSpace}
            ]
            ,
        nSets]
  ]
findCkRetrSupportAbsImposeQ[wfQRSpaceFullSpace_, ckModel_, nodesNeighbourhoods_,
  wannierRectangleTableValues_, \[Delta]x_, \[Delta]y_, support_, nIterations_, nRepeats_, nHIO_, gamma_, nSets_,
  nEREnd_, nAbsImpose_,
  nAbsImposeStart_, nAbsImposeEnd_, q_, \[Sigma]w_,
  numCellsX_, numCellsY_, nodesExactPositions_, elementaryCellXYTable_,
  fullSpaceXYTable_, ckSupportMemberTable_, SNR_, FBZnColnRow_, \[Sigma]resolution_] :=
    Module[{
      ckRetr,
      ckRetrMirror,
      retr,
      overlapRetr,
      overlapRetrMirror,
      distKSpace,
      measuredAbsSq = GaussianFilter[addNoise[Abs[Fourier[wfQRSpaceFullSpace]]^2, SNR, FBZnColnRow],
        {If[3 \[Sigma]resolution > 1, 3 \[Sigma]resolution, 1], \[Sigma]resolution}]
    },
      Return@Table[
        retr = phaseRetrieveSupportAbsImpose[Sqrt[measuredAbsSq], support, nIterations, nRepeats, nHIO, gamma, nEREnd,
          nAbsImpose,
          nAbsImposeStart, nAbsImposeEnd, nodesNeighbourhoods,
          wannierRectangleTableValues, \[Delta]x, \[Delta]y, q, \[Sigma]w,
          numCellsX, numCellsY, nodesExactPositions, elementaryCellXYTable,
          fullSpaceXYTable, ckSupportMemberTable];
        distKSpace = 1/Total[Total[measuredAbsSq]]*(Total@Total[Abs[Sqrt[measuredAbsSq]-Abs[Fourier[retr]]]]);

        ckRetr = wannierProject[retr, nodesNeighbourhoods, wannierRectangleTableValues, \[Delta]x, \[Delta]y];
        overlapRetr = Abs[ckRetr.Conjugate[ckModel]]^2/Abs[ckModel.Conjugate[ckModel]]^2;

        (*ckRetrMirror = wannierProject[mirrorXY[retr], nodesNeighbourhoods, wannierRectangleTableValues, \[Delta]x,
          \[Delta]y];
        overlapRetrMirror = Abs[ckRetrMirror.Conjugate[ckModel]]^2;

        If[overlapRetr > overlapRetrMirror,
          {ckRetr, overlapRetr, distKSpace},
          {ckRetrMirror, overlapRetrMirror, distKSpace}
        ]*)
        {ckRetr, overlapRetr, distKSpace}
        ,
        nSets]
    ]




End[] (* `Private` *)

EndPackage[]