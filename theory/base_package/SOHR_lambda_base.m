(* ::Package:: *)

LambdaHij[i_, j_] := If[i < j, 0, If[i == 1 && j == 1, 1, 
      1/Product[If[n == j, 1, Symbol[StringJoin["\[Lambda]", ToString[n]]] - 
          Symbol[StringJoin["\[Lambda]", ToString[j]]]], {n, 1, i}]]]
 
LambdaGetExponentialBase[q_] := 
    Table[E^(-Symbol[StringJoin["\[Lambda]", ToString[i]]]), {i, 1, q}]
 
LambdaGetHMatrix[m_, n_] := Table[LambdaHij[i, j], {i, 1, m}, {j, 1, n}]
 
LambdaGetPathwayProbabilityVector[l_] := 
    Module[{HMatrix, ExponentialBasis, output}, 
     HMatrix = LambdaGetHMatrix[l, l]; ExponentialBasis = 
       LambdaGetExponentialBase[l]; output = HMatrix . ExponentialBasis; 
      Return[output]; ]
 
LambdaCalculateLimit[P_, i_, j_] := 
    Limit[P, Symbol[StringJoin["\[Lambda]", ToString[i]]] -> 
      Symbol[StringJoin["\[Lambda]", ToString[j]]]]
 
LambdaToSubScriptExpression[P_, i_] := Module[{output}, 
     output = P /. {Symbol[StringJoin["\[Lambda]", ToString[i]]] -> 
          Subscript["\[Lambda]", ToString[i]]}; Return[output]; ]
