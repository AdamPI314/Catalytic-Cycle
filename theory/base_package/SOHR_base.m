(* ::Package:: *)

Hij[i_, j_] := If[i < j, 0, If[i == 1 && j == 1, 1, 
      Product[Symbol[StringJoin["k", ToString[m]]], {m, 1, i - 1}]/
       Product[If[n == j, 1, Symbol[StringJoin["k", ToString[n]]] - 
          Symbol[StringJoin["k", ToString[j]]]], {n, 1, i}]]]
 
GetExponentialBase[q_] := Table[E^((-Symbol[StringJoin["k", ToString[i]]])*
       t), {i, 1, q}]
 
GetHMatrix[m_, n_] := Table[Hij[i, j], {i, 1, m}, {j, 1, n}]
 
GetPathwayProbabilityVector[l_] := Module[{HMatrix, ExponentialBasis, 
      output}, HMatrix = GetHMatrix[l, l]; ExponentialBasis = 
       GetExponentialBase[l]; output = HMatrix . ExponentialBasis; 
      Return[output]; ]
 
CalculateLimit[P_, i_, j_] := Limit[P, 
     Symbol[StringJoin["k", ToString[i]]] -> 
      Symbol[StringJoin["k", ToString[j]]]]
 
ToSubScriptExpression[P_, i_] := Module[{output}, 
     output = P /. {Symbol[StringJoin["k", ToString[i]]] -> 
          Subscript["k", ToString[i]]}; Return[output]; ]
