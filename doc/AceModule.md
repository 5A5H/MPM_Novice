# Docu for adding a AceMate

## Generating a Material Subroutine with AceGen:

* first we generate a routine
```mathematica
<< AceGen`;
SMSInitialize["SmallStrainHookePlaneStress2D", "Language" -> "C++",
  "Mode" -> "Optimal"];
SMSModule["SmallStrainHookePlaneStress2D",
  Real[MaterialData$$[2], Eps$$[3], Exp$$[16]]];
{Emod, \[Nu]} \[RightTee]
  SMSReal[{MaterialData$$[1], MaterialData$$[2]}];
{\[Lambda], \[Mu]} \[DoubleRightTee] SMSHookeToLame[Emod, \[Nu]];
\[CurlyEpsilon]vec \[RightTee] SMSReal[{Eps$$[1], Eps$$[2], Eps$$[3]}];
\[CurlyEpsilon]33 = -\[Lambda] (\[CurlyEpsilon]vec[[
       1]] + \[CurlyEpsilon]vec[[2]])/(\[Lambda] + 2 \[Mu]);
SMSFreeze[\[CurlyEpsilon], {{\[CurlyEpsilon]vec[[
     1]], \[CurlyEpsilon]vec[[3]],
    0}, {\[CurlyEpsilon]vec[[3]], \[CurlyEpsilon]vec[[2]], 0}, {0,
    0, \[CurlyEpsilon]33}}, "Symmetric" -> True];
\[CapitalPi] \[DoubleRightTee] \[Lambda]/
    2 Tr[\[CurlyEpsilon]]^2 + \[Mu] \
Tr[\[CurlyEpsilon].\[CurlyEpsilon]];
\[Sigma] \[DoubleRightTee]
  SMSD[\[CapitalPi], \[CurlyEpsilon], "Symmetric" -> True];
\[Sigma]vec = {\[Sigma][[1, 1]], \[Sigma][[2, 2]], \[Sigma][[3,
    3]], \[Sigma][[1, 2]]};
SMSExport[\[Sigma]vec, Table[Exp$$[i], {i, 4}]];
\[DoubleStruckD]\[Sigma]\[DoubleStruckD]\[CurlyEpsilon] =
  Table[SMSD[\[Sigma]vec[[i]], \[CurlyEpsilon]vec[[j]]], {i, 4}, {j,
    3}];
SMSExport[
  Flatten[\[DoubleStruckD]\[Sigma]\[DoubleStruckD]\[CurlyEpsilon]],
  Table[Exp$$[i], {i, 5, 16}]];
SMSWrite["LocalAuxiliaryVariables" -> True];
```
* then we need to add the headder manually.
```
void SmallStrainHooke(double MaterialData[2],double Eps[3],double Exp[16]);
```
* headder need to be included in the main file (or wehere it is called)
```c++
#include <SmallStrainHookePlaneStress2D.hpp>
```
to use this include procedure the cmake file was extended by the directory
```
target_include_directories(
  mdm2dp
  PUBLIC
    "include"
    "ace_mate")
```
