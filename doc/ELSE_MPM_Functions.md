# Else Functions:
## Function of an MPM Material:
An MPM Material is pre-defined by the Material Base class. This lies in the namespace `ELSE::MPM::Material::`. The base interface provides the following functions:

 `addMaterialParameter()`

 This function is used to add a material parameter. A material parameter is defined by a name of type `string` and a value which may be `int`, `double` or `bool`.

 The following represents an example of how to use and interact with a Material object.
 ```
 ELSE::MPM::Material Mate1("Steel");
 Mate1.addMaterialParameter("Emod",21000.0e0);
 Mate1.addMaterialParameter("nue",0.3e0);
 Mate1.addMaterialParameter("rho",1e3);
 Mate1.addMaterialParameter("integervalue",1);
 Mate1.addMaterialParameter("boolvalue",true);
 Mate1.dumpMaterialParameter(ELSE::LogFile);

 int testint = 0;
 Mate1.getMaterialParameter("integervalue",testint);
 std::cout << "Int   : " << testint << std::endl;
 double testdbl = 1.0;
 Mate1.getMaterialParameter("Emod",testdbl);
 std::cout << "Double: " << testdbl << std::endl;
 bool testbool = false;
 Mate1.getMaterialParameter("boolvalue",testbool);
 std::cout << "Bool  : " << testbool << std::endl;

 std::array<double, 6> Sig = {0,0,0,   0,0,     0};
 std::array<double, 9> F   = {1,0,0, 0,1,0, 0,0,1};
 std::map<std::string, double> MaterialHistory;
 std::map<std::string, int>    IntegerMaterialIO;
 std::map<std::string, double> DoubleMaterialIO;
 Mate1.getCauchyStress(Sig,F,MaterialHistory,IntegerMaterialIO,DoubleMaterialIO);
 ```
The first lines make use of the
