(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)

(* :Title: Protocoling *)
(* :Context: Protocoling` *)
(* :Author: c8888 *)
(* :Date: 2017-11-20 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: *)
(* :Copyright: (c) 2017 c8888 *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["Protocoling`"]
(* Exported symbols added here with SymbolName::usage *)

protocolAdd::usage =
    "protocolAdd[stringMessage] adds stringMessage at the end of protocol file"
protocolMaxMemoryUsed::usage =
    "protocolMaxMemoryUsed[] adds information about max memory used (in GB) in current kernel session"
protocolBar::usage =
    "adds a horizontal bar to the protocol"

Begin["`Private`"]

protocolAdd[stringMessage_]:=PutAppend[stringMessage, "out/" <> ToString[Last@$CommandLine] <> "_"  <> ToString[$ProcessID] <> "protocol.txt"]
protocolMaxMemoryUsed[]:=protocolAdd["Max memory used (GB): " <> ToString[MaxMemoryUsed[]/1024.^3]]
protocolBar[]:=Module[{}, protocolAdd[" "]; protocolAdd["*************************************************"]; protocolAdd[" "];]

End[] (* `Private` *)

EndPackage[]