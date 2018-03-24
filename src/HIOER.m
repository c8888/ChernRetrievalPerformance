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
                and after some number of iterations returns retrieved object in real space.*)

BeginPackage["HIOER`"]
(* Exported symbols added here with SymbolName::usage *)

phaseRetrieveSupport::usage =
    "phaseRetrieveSupport[FTXAbs, support, nIterations, nRepeats, nHIO, gamma] returns a table of retrieved complex X from its (representation in reciprocal space) absolute value FTXAbs.
    Support is the table of the same structure as FTXAbs, it has rough values {0,1} for each {kx,ky} point. Procedure performs nIterations HIOER iterations. ER is made every nHIO
iterations not to stay in local minima."
phaseRetrieveGuess::usage =
    "phaseRetrieveGuess[FTXAbs, XAbsGuess, support, nIterations, nRepeats, nHIO, gamma] returns a table of retrieved complex X from its reciprocal space representation absolute value FTXAbs.
     Uses XAbsGuess to converge faster.ER is made every nHIO iterations not to stay in local minima."

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


End[] (* `Private` *)

EndPackage[]