(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)

(* :Title: Space2D *)
(* :Context: Space2D` *)
(* :Author: c8888 *)
(* :Date: 2017-11-19 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: *)
(* :Copyright: (c) 2017 c8888 *)
(* :Keywords: *)
(* :Discussion: Here all spatial structures are created. *)

BeginPackage["Space2D`"]
(* Exported symbols added here with SymbolName::usage *)

elementaryCell2D::usage =
    "elementaryCell2D[xmin_, xmax_, ymin_, ymax_, \[Delta]x_, \[Delta]y_] Constructs 2D table with point coordinates for elementary cell.
    Should me matched perfectly with a chosen wave function."

connect2DElementaryCells::usage =
    "connect2DElementaryCells[numCellsX_, numCellsY_, elementaryCell2DValues_] connects elementary cells and constructs
    table with (2numCellsX + 1) times (2numCellsY + 1) cells."

addSupportRectAndDimensionalize::usage =
    "addSupportRectAndDimensionalize[connect2DElementaryCellsTable_, marginSizeMinPercentageX_, marginSizeMinPercentageY_]:
    connect2DElementaryCellsTable - 2D array that will be padded with support
    marginSizeMinPercentageX - minimal size of total padding in percent of current size in x
    marginSizeMinPercentageY - minimal size of total padding in percent of current size in y

    This function returns new array space2DProbingPoints of size 2^N1 times 2^N2 where N1, N2 are the smallest integers with which new array
    can include old array and the minimal margins in x and y, respectively. "




Begin["`Private`"]

elementaryCell2D[xmin_, xmax_, ymin_, ymax_, \[Delta]x_, \[Delta]y_]:=
    Table[{x,y}, {x, N@xmin, N@xmax, N@\[Delta]x}, {y, N@ymin, N@ymax, N@\[Delta]y}]

(* rectSupport[numCellsX_, numCellsY_]:= *)

connect2DElementaryCells[numCellsX_, numCellsY_, elementaryCell2DValues_]:=
     Return@Flatten[Table[Transpose[Flatten[Table[elementaryCell2DValues, 2 numCellsX + 1], 1]], 2 numCellsY + 1], 1];

addSupportRectAndDimensionalize[connect2DElementaryCellsTable_, marginSizeMinPercentageX_, marginSizeMinPercentageY_]:=
    Module[
      {
        dims = Dimensions@connect2DElementaryCellsTable,
        newDims,
        N1=0,
        N2=0
      },
      If[Dimensions@dims != 2, Throw["Space2D_addSupportRectAndDimensionalize_Bad_Dimensions"]];
      newDims = Floor[{1 + marginSizeMinPercentageX / 100, 1 + marginSizeMinPercentageY / 100 } * dims];

      While[2^N1 <= newDims[[1]], N1++];
      newDims[[1]] = 2^(N1);

      While[2^N2 <= newDims[[2]], N2++];
      newDims[[2]] = 2^(N2);

      Return@CenterArray[connect2DElementaryCellsTable, newDims];
    ]

End[] (* `Private` *)

EndPackage[]