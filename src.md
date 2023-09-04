```Mathematica
mynormf[{a_, b_}] := a^2 + b^2; 
ClearAll[a, b, h, H];
(*major axes*)
a[t_] = 2 + 0.1 Sin[2 \[Pi] t]; 
b[t_] = 1.5 + 0.2 Cos[2 \[Pi] t];
(*parameter curve*)
para[{s_, t_}] = {a[t] Cos[2 \[Pi] s], b[t] Sin[2 \[Pi] s]};
(*the frequency*)
p1 = 4; p2 = 2; q = 19;
(*the generation function*)
h[{{s0_, t0_}, {s1_, t1_}}] = 
  mynormf[para[{s1, t1}] - para[{s0, t0}]]*(1/(t1 - t0));
(*variable: {{s_0,t_0},...,{s_{q-1},t_{q-1}}}*)
var = Table[{ToExpression[ToString[s] <> ToString[i]], 
    ToExpression[ToString[t] <> ToString[i]]}, {i, 0, q - 1}];
(*the total action*)
H = Total[h /@ Partition[var, 2, 1]] + 
   h[{Last[var], {s0 + p1, t0 + p2}}];
(*the grad of the total action*)
dd = Grad[H, Flatten@var];
(*{0,...,0} length = 2q*)
zeros = Table[0, {i, 1, Length@Flatten@var}];
(*the initial condition*)
initial = 
 Join[Transpose@{var[[All, 1]], Table[i*(p1/q) + 0.2, {i, 0, q - 1}]},
   Transpose@{var[[All, 2]], Table[i*(p2/q) + 0.7, {i, 0, q - 1}]}];
(*generate the equations*)
eqgen[l_, l1_] := Table[l[[i]] == l1[[i]], {i, 1, Length@l}];
(*the equations*)
eq = eqgen[dd, zeros];
(*solve dd=0*)
sol = Quiet@
   FindRoot[eq, initial, MaxIterations -> 10000, AccuracyGoal -> 16];
points = (para /@ ((Join[var, {{s0 + p1, t0 + p2}}]) /. sol));

(*We want to verify that our sol is a real solution.
So we pick the first point of sol as the initial condition and then check whether it is a real solution*)

(*our solution*)
data = var /. sol;
(*the initial conditions x0 and v0*)
x0 = ((para@data[[2]]) + (para@data[[1]]))/
  2; v0 = ((para@data[[2]]) - (para@data[[1]]))*(1/(t1 - t0)) /. sol;
(*region of the billiard boundary*)
region[x_, y_, s_] := x^2/(a[s])^2 + y^2/(b[s])^2 - 1;
(*the unit normal vector*)
normal[x_, y_, s_] = Normalize[Grad[region[x, y, s], {x, y}]];
(*the unit tangent vector*)
tangent[x_, y_, s_] = RotationMatrix[\[Pi]/2] . normal[x, y, s];
(*the normal velocity of the billiard boundary*)
normalv[x_, y_, s_] = 
  normal[x, y, s] (normal[x, y, s] . {a'[s] (x/a[s]), b'[s] (y/b[s])});
(*the initial time*)
tt0 = (t0 + (t1 - t0)/2) /. sol;
(*integration time*)
time = p2*2;
(*the numerical billiard solution*)
aa = NDSolve[{x''[t] == 0, y''[t] == 0, s'[t] == 1, s[0] == tt0, 
    x[0] == x0[[1]], y[0] == x0[[2]], x'[0] == v0[[1]], 
    y'[0] == v0[[2]],
    WhenEvent[region[x[t], y[t], s[t]] == 0,
     {x'[t], y'[t]} -> (
       ({x'[t], y'[t]} . tangent[x[t], y[t], s[t]]) tangent[x[t], 
          y[t], s[t]] +
        (-({x'[t], y'[t]} . normal[x[t], y[t], s[t]]) normal[x[t], 
            y[t], s[t]] + 2 normalv[x[t], y[t], s[t]])
       )]
    
    }, {x, y, s}, {t, 0, time}, MaxStepSize -> 0.001, 
   MaxSteps -> Infinity, AccuracyGoal -> 15];
(*Plot both the variational solution and the billiard solution*)
{Show[ParametricPlot[{x[t], y[t]} /. aa, {t, 0, time}, 
    PlotStyle -> {Thickness[0.004], Red}, PlotPoints -> 1000, 
    Axes -> False], 
   Graphics[{Arrowheads[{{.02, 0.55}}], 
     Arrow /@ Partition[points, 2, 1]}], PlotRange -> All, 
   ImageSize -> Medium],  
  boundaries = 
   Table[ParametricPlot[para[{s, i}], {s, 0, 1}, 
     PlotStyle -> {Red, Thickness[0.001]}, Axes -> False, 
     Frame -> False], {i, var[[All, 2]] /. sol}];
  Show[Graphics[{Arrowheads[{{.02, 0.55}}], 
     Arrow /@ Partition[points, 2, 1]}], boundaries, PlotRange -> All,
    ImageSize -> Medium]} // Row
```
