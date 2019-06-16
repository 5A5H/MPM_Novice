(*************************************************************
* AceGen    6.923 MacOSX (19 Apr 19)                         *
*           Co. J. Korelc  2013           15 Jun 19 21:56:58 *
**************************************************************
User     : Full professional version
Notebook : mate_planestrain_linearelasticity_ctest
Evaluation time                 : 0 s     Mode  : Optimal
Number of formulae              : 6       Method: Automatic
Module                          : mate01 size: 205
Total size of Mathematica  code : 205 subexpressions       *)
(*********************** M O D U L E **************************)
SetAttributes[mate01,HoldAll];
mate01[MaterialData$$_,Eps$$_,Exp$$_]:=Module[{},
$VV[28]=MaterialData$$[[1]]/(1+MaterialData$$[[2]]);
$VV[3]=(MaterialData$$[[2]]*$VV[28])/(1-2*MaterialData$$[[2]]);
$VV[25]=$VV[3]+$VV[28];
$VV[20]=(Eps$$[[1]]+Eps$$[[3]])*$VV[3];
Exp$$[[1]]=$VV[20]+Eps$$[[1]]*$VV[28];
Exp$$[[2]]=$VV[20]+Eps$$[[3]]*$VV[28];
Exp$$[[3]]=$VV[20];
Exp$$[[4]]=Eps$$[[2]]*$VV[28];
Exp$$[[5]]=$VV[25];
Exp$$[[6]]=0;
Exp$$[[7]]=$VV[3];
Exp$$[[8]]=$VV[3];
Exp$$[[9]]=0;
Exp$$[[10]]=$VV[25];
Exp$$[[11]]=$VV[3];
Exp$$[[12]]=0;
Exp$$[[13]]=$VV[3];
Exp$$[[14]]=0;
Exp$$[[15]]=$VV[28];
Exp$$[[16]]=0;
];
