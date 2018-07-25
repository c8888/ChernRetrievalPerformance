(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)

(* :Title: HIOER *)
(* :Context: HIOER` *)
(* :Author: c8888 *)
(* :Date: 2017-11-19 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 11.1 *)
(* :Copyright: (c) 2017 c8888 *)
(* :Keywords: *)
(* :Discussion: Implementations of Hybrid Input-Output Algorithm. Takes a table with absolute value in Fourier 2D space
                and after a number of iterations returns retrieved object in real space.*)

BeginPackage["HIOER`"]
(* Exported symbols added here with SymbolName::usage *)
Needs["WaveFunctions`"];
Needs["Space2D`"];

phaseRetrieveSupport::usage =
    "phaseRetrieveSupport[FTXAbs, support, nIterations, nRepeats, nHIO, gamma] returns a table of retrieved complex X from its (representation in reciprocal space) absolute value FTXAbs.
    Support is the table of the same structure as FTXAbs, it has rough values {0,1} for each {kx,ky} point. Procedure performs nIterations HIOER iterations. ER is made every nHIO
iterations not to stay in local minima."
phaseRetrieveSupportAbsImpose::usage =
    "phaseRetrieveSupportAbsImpose[FTXAbs_, support_, nIterations_, nRepeats_, nHIO_, gamma_, nEREnd_, nAbsImpose_,
  nAbsImposeStart_, nAbsImposeEnd_, nodesNeighbourhoods_, wannierRectangleTableValues_, \[Delta]x_, \[Delta]y_, q_,
  \[Sigma]w_, numCellsX_,  numCellsY_, nodesExactPositions_, elementaryCellXYTable_, fullSpaceXYTable_,
  ckSupportMemberTable_, ax_, ay_] imposes
  condition that wannier heights should be equal in the whole lattice. Also inputs the information about the symmetry
  of the lattice."
phaseRetrieveGuess::usage =
    "phaseRetrieveGuess[FTXAbs, XAbsGuess, support, nIterations, nRepeats, nHIO, gamma] returns a table of retrieved complex X from its reciprocal space representation absolute value FTXAbs.
     Uses XAbsGuess to converge faster.ER is made every nHIO iterations not to stay in local minima."
phaseRetrieveSupportStartHint::usage =
    "phaseRetrieveSupportStartHint[FTXAbs_, StartHintX_, support_, nIterations_, nRepeats_, nHIO_, gamma_,
  nIterationsOverlapCalc_, ckStartHintAbs_, nodesNeighbourhoods_,
  wannierRectangleTableValues_, \[Delta]x_, \[Delta]y_] returns table
    of a
    retrieved object and overlaps of retrieval MODULUS and startHint MODULUS every nIterationsOverlapCalc steps. Takes
     as a starting point the complex value in the real space."

Begin["`Private`"]

phaseRetrieveGuess[FTXAbs_, wavefAbs_, support_, nIterations_, nRepeats_, nHIO_, gamma_]:= (* returns table of a
    retrieved object *)
    Module[ {
      nCol,
      nRow,
      xi={},
      xiprim={},
      xierror,
      retrerror,
      retr={},
      FTxi={},
      inverseSupport
    },
      {nCol, nRow} = Dimensions[FTXAbs];
      inverseSupport = -(support-1);

      Do[
        xi=Table[RandomComplex[], nCol, nRow]; (* random initialization, different complex numbers at each repetition *)
        Do[
        (*protocolAdd[{"Repeat, Iteration: ", {k,i}}];*)
          xiprim = xi;
          FTxi = Fourier[xi];
          FTxi = FTXAbs*Exp[I*Arg[FTxi]];
          xi = InverseFourier[FTxi];
          If[Mod[i, nHIO]!=0,
          (*HIO case*)
            xi = inverseSupport * (xiprim - gamma xi) + wavefAbs * Exp[I Arg[xi] ]
            (* we assume that wavefAbs is properly "supported"*)
            ,
            (*ER case*)
            xi = wavefAbs * Exp[I Arg[xi]]
          ];
          ,
          {i, nIterations}
        ];
        xierror=Total@Total@Abs[Abs[Fourier[xi]]^2-Abs[FTXAbs]^2];
        If[k == 1, retrerror=xierror]; (*first repetition*)
        If[xierror<=retrerror, retr=xi; retrerror=xierror]; (* error estimator *)
      (* backup mod *) (*Export["retrieved_insite_nRepeat="<>ToString[k]<>"_nIterations="<>ToString[nIterations]<>"_RTF="<>ToString[RTF]<>"_sigma_n="<>ToString[\[Sigma]n]<>".dat", xi];*)
        ,
        {k, nRepeats}
      ];
      Return[retr*support] (* returned value *)

    ]


phaseRetrieveSupport[FTXAbs_, support_, nIterations_, nRepeats_, nHIO_, gamma_]:= (* returns table of a retrieved object *)
    Module[ {
      nCol,
      nRow,
      xi={},
      xiprim={},
      xierror,
      retrerror,
      retr={},
      FTxi={},
      inverseSupport
    },
      {nCol, nRow} = Dimensions[FTXAbs];
      inverseSupport = -(support-1);

      Do[
        xi=Table[RandomComplex[], nCol, nRow]; (* random initialization, different complex numbers at each repetition *)
        Do[
        (*protocolAdd[{"Repeat, Iteration: ", {k,i}}];*)
          xiprim = xi;
          FTxi = Fourier[xi];
          FTxi = FTXAbs*Exp[I*Arg[FTxi]];
          xi = InverseFourier[FTxi];

          If[Mod[i, nHIO]!=0,
          (*HIO case*)
            xi = inverseSupport * (xiprim - gamma xi) + support * xi;
            ,
          (*ER case*)
            xi *= support;
          ];,
          {i, nIterations}
        ];
        xierror=Total@Total@Abs[Abs[Fourier[xi]]^2-Abs[FTXAbs]^2];
        If[k == 1, retrerror=xierror]; (*first repetition*)
        If[xierror<=retrerror, retr=xi; retrerror=xierror]; (* error estimator *)
      (* backup mod *) (*Export["retrieved_insite_nRepeat="<>ToString[k]<>"_nIterations="<>ToString[nIterations]<>"_RTF="<>ToString[RTF]<>"_sigma_n="<>ToString[\[Sigma]n]<>".dat", xi];*)
        ,
        {k, nRepeats}
      ];
      Return[retr*support]; (* returned value *)

    ]

phaseRetrieveSupportAbsImpose[FTXAbs_, support_, nIterations_, nRepeats_, nHIO_, gamma_, nEREnd_, nAbsImpose_,
  nAbsImposeStart_, nAbsImposeEnd_, nodesNeighbourhoods_, wannierRectangleTableValues_, \[Delta]x_, \[Delta]y_, q_,
  \[Sigma]w_, numCellsX_,  numCellsY_, nodesExactPositions_, elementaryCellXYTable_, fullSpaceXYTable_,
  ckSupportMemberTable_, ax_, ay_]:= (*
    returns
    table of a
    retrieved object *)
    Module[ {
      dimx,
      dimy,
      xi={},
      xiprim={},
      xierror,
      retrerror,
      retr={},
      FTxi={},
      inverseSupport
    },
      {dimx, dimy} = Dimensions[FTXAbs];
      inverseSupport = -(support-1);

      Do[
        xi=Table[RandomComplex[], dimx, dimy]; (* random initialization, different complex numbers at each repetition *)
        Do[
        (*protocolAdd[{"Repeat, Iteration: ", {k,i}}];*)
          xiprim = xi;
          FTxi = Fourier[xi];
          FTxi = FTXAbs*Exp[I*Arg[FTxi]];
          xi = InverseFourier[FTxi];

          If[Mod[i, nHIO]!=0 && (nIterations - i) > nEREnd,
          (*HIO case*)
            xi = inverseSupport * (xiprim - gamma xi) + support * xi;
            ,
          (*ER case*)
            xi *= support;
          ];
          If[Mod[i, nAbsImpose] == 0 && i >= nAbsImposeStart && i <= nAbsImposeEnd,
            xi = Exp[I Arg[xi]]*
                fastFullSpaceWfQRSpaceAbs[
                  meanCkAbs[
                    wannierProject[xi, nodesNeighbourhoods,
                      wannierRectangleTableValues, \[Delta]x, \[Delta]y],
                    2 q, ckSupportMemberTable], \[Sigma]w, numCellsX, numCellsY, nodesExactPositions,
                  elementaryCellXYTable,
                  fullSpaceXYTable, \[Delta]x, \[Delta]y, dimx, dimy, support, ax, ay]
          ];,
          {i, nIterations}
        ];
        xierror=Total@Total[Abs[Abs[Fourier[xi]]-Abs[FTXAbs]]^2];
        If[k == 1, retrerror=xierror]; (*first repetition*)
        If[xierror<=retrerror, retr=xi; retrerror=xierror]; (* error estimator *)
        ,
        {k, nRepeats}
      ];
      Return[retr*support]; (* returned value *)
    ]

phaseRetrieveSupportStartHint[FTXAbs_, StartHintX_, support_, nIterations_, nRepeats_, nHIO_, gamma_,
  nIterationsOverlapCalc_, ckStartHintAbs_, nodesNeighbourhoods_,
  wannierRectangleTableValues_, \[Delta]x_, \[Delta]y_]:= (*
    returns table
    of a
    retrieved object and overlaps of retrieval MODULUS and startHint MODULUS every nIterationsOverlapCalc steps *)
    Module[ {
      nCol,
      nRow,
      xi={},
      xiprim={},
      xierror,
      retrerror,
      FTxi={},
      inverseSupport,
      ret={{},{}},
      ckRetr
    },
      {nCol, nRow} = Dimensions[FTXAbs];
      inverseSupport = -(support-1);

      Do[
        xi=(*Table[RandomComplex[], nCol, nRow]*)StartHintX; (* initialization *)
        Do[
        (*protocolAdd[{"Repeat, Iteration: ", {k,i}}];*)
          xiprim = xi;
          FTxi = Fourier[xi];
          FTxi = FTXAbs*Exp[I*Arg[FTxi]];
          xi = InverseFourier[FTxi];

          If[Mod[i, nHIO]!=0,
          (*HIO case*)
            xi = inverseSupport * (xiprim - gamma xi) + support * xi;
            ,
          (*ER case*)
            xi *= support;
          ];
          If[Mod[i, nIterationsOverlapCalc]==0,
            AppendTo[ret[[2]],
              {
                i,
                Max[{overlapWannier[ckStartHintAbs, wannierProject[Abs[xi], nodesNeighbourhoods,
                      wannierRectangleTableValues,
                      \[Delta]x,
                      \[Delta]y]],
                    overlapWannier[ckStartHintAbs,  wannierProject[mirrorXY[Abs[xi]], nodesNeighbourhoods,
                      wannierRectangleTableValues,
                      \[Delta]x,
                      \[Delta]y]]
                    }]
              }
            ];
          ];
          ,
          {i, nIterations}
        ];
        xierror=Total@Total@Abs[Abs[Fourier[xi]]^2-Abs[FTXAbs]^2];
        If[k == 1, retrerror=xierror]; (*first repetition*)
        If[xierror<=retrerror, ret[[1]]=xi; retrerror=xierror]; (* error estimator *)
      (* backup mod *) (*Export["retrieved_insite_nRepeat="<>ToString[k]<>"_nIterations="<>ToString[nIterations]<>"_RTF="<>ToString[RTF]<>"_sigma_n="<>ToString[\[Sigma]n]<>".dat", xi];*)
        ,
        {k, nRepeats}
      ];
      Return[{ret[[1]]*support, ret[[2]]}]; (* returned value *)

    ]


End[] (* `Private` *)

EndPackage[]