(* Mathematica Source File  *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)
(* :Author: c8888 *)
(* :Date: 2018-03-24 *)

Needs["Space2D`", "../src/Space2D.m"];
Needs["WaveFunctions`", "../src/WaveFunctions.m"];
Needs["HIOER`", "../src/HIOER.m"];
Needs["ChernCalc`", "../src/ChernCalc.m"];
Needs["Protocoling`", "../src/Protocoling.m"];

SetSystemOptions["ParallelOptions" -> "ParallelThreadNumber" -> 1];
SetSystemOptions["ParallelOptions" -> "MKLThreadNumber" -> 1];

(**************************************************************)
(*DEFAULT PARAMETERS VALUES*)

numCellsX = 1;
numCellsY = 1;
q = 3;
gx = 2; gy = 2;
\[Sigma]w = 0.155;
rangeRectangleSizeX = 1.;
rangeRectangleSizeY = 1.;
a = 1;
J = 1;
J1 = 2;
nIterations = 5;
nRepeats = 1;
nHIO = 20;
gamma = 0.9;
npts = 3;(*points in the 1st Brillouin zone*)
n = 1; (* band number 1...q *)
nSets = 2; (* Number of separate phase retrievals. Each phase retrieval gives a single chern number. *)
nEREnd = 30;
nAbsImpose = 3;
nAbsImposeStart = 1;

trapezeRatio = 0.8;
folder = "default";
SNR = 10000000000000000;
(**************************************************************)

ToExpression[$CommandLine[[4;;-1]]]; (*Execute parameter overriding from command line. HAS TO BE RUN WITH
WOLFRAMSCRIPT: wolframscript -file chernHistogram.m npts=8*)
extraArg = StringJoin@Map[ToString[#]&, $CommandLine[[4;;-1]]];
If[folder=="default", folder=extraArg];

dir = "../out/" <> ToString[folder];

CreateDirectory[dir];
protocolSetDir[dir];
protocolSetExtraArg[extraArg];
t1 = DateList[];
protocolAdd[ToString[t1] <> " Program started."];
(**************************************************************)
protocolBar[];
protocolAdd["Parameters: "];
protocolAdd["numCellsX = "<> ToString[numCellsX] ];
protocolAdd["numCellsY = "<> ToString[numCellsY] ];
protocolAdd["q = "<> ToString[q] ];
protocolAdd["dimx = "<> ToString[dimx] ];
protocolAdd["dimy = "<> ToString[dimy] ];
protocolAdd["\[Delta]x = "<> ToString[\[Delta]x] ];
protocolAdd["\[Delta]y = "<> ToString[\[Delta]y] ];
protocolAdd["gx = "<> ToString[gx] ];
protocolAdd["gy = "<> ToString[gy] ];
protocolAdd["\[Sigma]w = "<> ToString[\[Sigma]w ] ];
protocolAdd["rangeRectangleSizeX = "<> ToString[rangeRectangleSizeX] ];
protocolAdd["rangeRectangleSizeY = "<> ToString[rangeRectangleSizeY] ];
protocolAdd["a = "<> ToString[a] ];
protocolAdd["J = "<> ToString[J] ];
protocolAdd["J1 = "<> ToString[J1] ];
protocolAdd["nIterations = "<> ToString[nIterations] ];
protocolAdd["nRepeats = "<> ToString[nRepeats] ];
protocolAdd["nHIO = "<> ToString[nHIO] ];
protocolAdd["gamma = "<> ToString[gamma] ];
protocolAdd["npts = "<> ToString[npts] ];
protocolAdd["n = "<> ToString[n] ];
protocolAdd["nSets = "<> ToString[nSets] ];
protocolAdd["nEREnd = "<> ToString[nEREnd] ];
protocolAdd["nAbsImpose = "<> ToString[nAbsImpose] ];
protocolAdd["nAbsImposeStart = "<> ToString[nAbsImposeStart] ];
protocolAdd["nAbsImposeEnd = "<> ToString[nAbsImposeEnd] ];
protocolAdd["trapezeRatio = "<> ToString[trapezeRatio] ];


protocolBar[];

(**************************************************************)
(* INITIALISATION *)
nAbsImposeEnd = nIterations;
(*---SPACE---*)
dimx = (2 numCellsX + 1)*gx*q*10; dimy = (2*numCellsY + 1)* gy*10; (* dimx, dimy must be even numbers for fast FFT, at best \
"2^N" *)
{ax, ay} = {q a, a};
{dimBZx, dimBZy} = {(2 numCellsX + 1) gx, (2 numCellsY + 1) gy}; (*in PIXELS *)
{\[Delta]kx, \[Delta]ky} = {2 Pi/ax/gx/(2 numCellsX + 1), 2 Pi/ay/gy/(2 numCellsY + 1)}; (* in 1/[a] dimensions *)
{\[Delta]x, \[Delta]y} = {2 Pi/(dimx*\[Delta]kx), 2 Pi/(dimy * \[Delta]ky)};

BZ = latticeProbingPointsBZ[npts, a, q];

nodesExactPositions = Table[{-q/2 + 1/2 + i*a, 0}, {i, 0, q - 1, 1}];

elementaryCellXYTable =
    elementaryCell2D[-ax/2, ax/2, -ay/2, ay/2, \[Delta]x, \[Delta]y];

elCellNodes =
    elementaryCell2DNodes[elementaryCellXYTable, nodesExactPositions];
protocolAdd["elCellNodes = " <> ToString[elCellNodes]];

fullSpaceXYTable =
    connectAndTranslate2DElementaryCells[numCellsX, numCellsY,
      elementaryCellXYTable, \[Delta]x, \[Delta]y];

(*---Support construction---*)
{cellDimX, cellDimY} = {ax/\[Delta]x, ay/\[Delta]y};
{cellsSpaceX, cellsSpaceY} = {cellDimX*(2 numCellsX + 1), cellDimY*(2 numCellsY + 1)};
cellsRangeX = {dimx/2 + 1 - cellsSpaceX/2, dimx/2 + 1 + cellsSpaceX/2};
cellsRangeY = {dimy/2 + 1 - cellsSpaceY/2, dimy/2 + 1 + cellsSpaceY/2};

(*---Dot products---*)
fullSpaceNodesCellsSpace =
    fullSpace2DNodes[elCellNodes, elementaryCellXYTable, numCellsX,
      numCellsY];
nodesXYTable = nodesXY[fullSpaceXYTable, fullSpaceNodesCellsSpace];
fullSpaceNodes =
    nodesTranslateAfterDimensionalizing[
      fullSpaceNodesCellsSpace, {cellsSpaceX, cellsSpaceY}, {dimx, dimy}];
nodesNeighbourhoods =
    nodesNeighbourhood[fullSpaceNodes, \[Delta]x, \[Delta]y,
      rangeRectangleSizeX, rangeRectangleSizeY, {dimx, dimy}];
wNFactor =
    wannierNormalisationFactor[\[Sigma]w, \[Delta]x, \[Delta]y,
      rangeRectangleSizeX, rangeRectangleSizeY];
wannierRectangleTableValues =
    wannierRectangleTable[fullSpaceNodes, nodesNeighbourhoods,
      wNFactor, \[Sigma]w, \[Delta]x, \[Delta]y];

{support, ckSupportMemberTable} =
    trapezeSupport[trapezeRatio, cellsRangeX, cellsRangeY, cellsSpaceX,
      cellsSpaceY, fullSpaceNodes, dimx, dimy];

(**************************************************************)
(* NOISE *)
FBZnColnRow = {Round[{dimx/2 - dimBZx/2, dimx/2 + dimBZx/2 + 1}],
  Round[{dimy/2 - dimBZy/2, dimy/2 + dimBZy/2 + 1}]};



(**************************************************************)
(* PHASE RETRIEVAL *)
protocolAdd["$ProcessorCount = "<> ToString[$ProcessorCount]];
protocolBar[];

ckRetrSupportTable =
    Map[
      Module[{wf},
        (*findCkRetrSupportQ[(wf =
          fastFullSpaceWfQRSpace[#[[1]], #[[2]],
            myES[hamiltonianHarperQ[#[[1]], #[[2]], J, J1, q]][[2,
                n]], \[Sigma]w, numCellsX, numCellsY, nodesExactPositions,
            elementaryCellXYTable,
            fullSpaceXYTable, \[Delta]x, \[Delta]y, dimx, dimy]),
        wannierProject[wf, nodesNeighbourhoods,
          wannierRectangleTableValues, \[Delta]x, \[Delta]y],
        nodesNeighbourhoods,
        wannierRectangleTableValues, \[Delta]x, \[Delta]y, support,
        nIterations, nRepeats, nHIO, gamma, nSets]*)
        protocolAdd["{kx,ky}=" <> ToString[#]];
        findCkRetrSupportAbsImposeQ[(wf =
            fastFullSpaceWfQRSpace[#[[1]], #[[2]],
              myES[hamiltonianHarperQ[#[[1]], #[[2]], J, J1, q]][[2,
                  n]], \[Sigma]w, numCellsX, numCellsY, nodesExactPositions,
              elementaryCellXYTable,
              fullSpaceXYTable, \[Delta]x, \[Delta]y, dimx, dimy, support]),
          wannierProject[wf, nodesNeighbourhoods,
            wannierRectangleTableValues, \[Delta]x, \[Delta]y],
          nodesNeighbourhoods,
          wannierRectangleTableValues, \[Delta]x, \[Delta]y, support,
          nIterations, nRepeats, nHIO, gamma, nSets, nEREnd, nAbsImpose,
          nAbsImposeStart, nAbsImposeEnd, q, \[Sigma]w,
          numCellsX, numCellsY, nodesExactPositions, elementaryCellXYTable,
          fullSpaceXYTable, ckSupportMemberTable, SNR, FBZnColnRow]] &, BZ, {2}(*,
      DistributedContexts -> All*)];

blochNodePhaseConjugateTableBZ =
    Map[blochNodePhaseConjugateTable[#[[1]], #[[2]], nodesXYTable] &,
      BZ, {2}];

kBasisRetrSupportTable = Table[
  Table[{kProject[ckRetrSupportTable[[i, j, r, 1]], blochNodePhaseConjugateTableBZ[[i, j]], q, ckSupportMemberTable], ckRetrSupportTable[[i, j, r, 2]], ckRetrSupportTable[[i, j, r, 3]]}, {r, 1, Dimensions[ckRetrSupportTable][[3]]}],
  {i, Dimensions[BZ][[1]]}, {j, Dimensions[BZ][[2]]}];

(*FxyTRetrSupportTable =
    Table[FxyT[ckRetrSupportTable[[All, All, i, 1, All]]], {i, 1,
      Dimensions[ckRetrSupportTable][[3]]}];

wRetrSupportTable =
    1/(2 Pi I)*
        Table[Chop@Total@Total[FxyTRetrSupportTable[[i]]], {i, 1,
          Dimensions[ckRetrSupportTable][[3]]}];
          *)

(**************************************************************)
(* SAVE DATA *)

Export[dir <> "/" <> ToString[extraArg] <>"_" <> ToString[Last@$CommandLine] <> "_" <> ToString[$ProcessID] <> "ckRetrSupportTable
.mx",
  ckRetrSupportTable];
protocolAdd["File saved to " <> dir <> "/" <> ToString[extraArg] <>"_" <>ToString[Last@$CommandLine] <> "_" <>
    ToString[$ProcessID] <>
    "ckRetrSupportTable.mx"];

Export[dir <> "/" <> ToString[extraArg] <>"_" <>ToString[Last@$CommandLine] <> "_" <> ToString[$ProcessID] <> "kBasisRetrSupportTable.mx",
  kBasisRetrSupportTable];
protocolAdd["File saved to " <> dir <> "/" <> ToString[extraArg] <>"_" <>ToString[Last@$CommandLine] <> "_" <>
    ToString[$ProcessID] <>
    "kBasisRetrSupportTable.mx"];

(*Export["../out/" <>ToString[Last@$CommandLine] <> "_" <> ToString[$ProcessID] <> "FxyTRetrSupportTable.mx",
  FxyTRetrSupportTable];
protocolAdd["File saved to ../out/" <>ToString[Last@$CommandLine] <> "_" <> ToString[$ProcessID] <> "FxyTRetrSupportTable.mx"];
*)

(*Export["../out/" <>ToString[Last@$CommandLine] <> "_" <> ToString[$ProcessID] <> "wRetrSupportTable.mx",
  wRetrSupportTable];
protocolAdd["File saved to ../out/" <>ToString[Last@$CommandLine] <> "_" <> ToString[$ProcessID] <> "wRetrSupportTable.mx"];
*)

(* protocolAdd["wRetrSupportTable = " <> ToString[wRetrSupportTable]];*)
protocolBar[];

t2 = DateList[];
protocolMaxMemoryUsed[];
protocolAdd[ToString[t2] <> " Program evaluated successfully. Total time taken: YMDHMS " <> ToString[t2-t1]];