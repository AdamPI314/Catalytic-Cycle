PILambdaGetExponentialBase[l_] := 
    Table[E^((-Symbol[StringJoin["\[Lambda]", ToString[i]]])*t), {i, 1, l}]
 
i = 1
 
PILambdaGetPathIntegral[N_] := Module[{output, tlist, klist}, 
     tlist = Table[Symbol[StringJoin["t", IntegerString[i]]], {i, 1, N}]; 
      klist = Table[Symbol[StringJoin["\[Lambda]", IntegerString[i]]], 
        {i, 1, N + 1}]; output = Product[With[{}, 
         Evaluate[klist[[i]]/E^(klist[[i]]*tlist[[i]])]], {i, 1, N}]; 
      Simplify[output *= E^((-klist[[-1]])*(t - Sum[tlist[[i]], 
            {i, 1, N}]))]; For[i = N, i > 1, i--, 
       output = Integrate[output, {tlist[[i]], 0, 
          t - Sum[tlist[[j]], {j, 1, i - 1}]}]]; 
      If[N == 0, output = E^((-klist[[1]])*t), 
       output = Integrate[output, {t1, 0, t}]]; output = Expand[output]; 
      output = output /. {Exp[x_] :> Exp[Together[FullSimplify[x]]]}; 
      Return[output]; ]
 
PILambdaGetHMatrix[pathProb_] := Module[{l, ExponentialVector, CoeffMatrix}, 
     l = Length[pathProb]; ExponentialVector = PIGetExponentialBase[l]; 
      CoeffMatrix = Table[Table[Simplify[Coefficient[pathProb[[i]], 
           ExponentialVector[[j]], 1]], {j, 1, l}], {i, 1, l}]; 
      Return[CoeffMatrix]; ]
 
PIGetExponentialBase[l_] := Table[E^((-Symbol[StringJoin["k", ToString[i]]])*
       t), {i, 1, l}]
 
PILambdaGetPathwayProbabilityVector[l_] := Module[{output}, 
     output = Table[PIGetPathIntegral[i], {i, 0, l}]; Return[output]; ]
 
PIGetPathIntegral[N_] := Module[{output, tlist, klist}, 
     tlist = Table[Symbol[StringJoin["t", IntegerString[i]]], {i, 1, N}]; 
      klist = Table[Symbol[StringJoin["k", IntegerString[i]]], 
        {i, 1, N + 1}]; output = Product[With[{}, 
         Evaluate[klist[[i]]/E^(klist[[i]]*tlist[[i]])]], {i, 1, N}]; 
      output *= E^((-klist[[-1]])*(t - Sum[tlist[[i]], {i, 1, N}])); 
      For[i = N, i > 1, i--, output = Integrate[output, 
         {tlist[[i]], 0, t - Sum[tlist[[j]], {j, 1, i - 1}]}]]; 
      If[N == 0, output = E^((-klist[[1]])*t), 
       output = Integrate[output, {t1, 0, t}]]; output = Expand[output]; 
      output = output /. {Exp[x_] :> Exp[Together[FullSimplify[x]]]}; 
      Return[output]; ]
 
PILambdaToSubScriptExpression[P_, i_] := Module[{output}, 
     output = P /. {Symbol[StringJoin["\[Lambda]", ToString[i]]] -> 
          Subscript["\[Lambda]", ToString[i]]}; Return[output]; ]
 
Attributes[Subscript] = {NHoldRest}
