// Class to provide and handle Ansatz/Shape functions for a Q4
// Q4 -> four noded polygon in 2d

#include <iostream>
#include <iomanip>

class MPMSHPQ4 {
  public:
      MPMSHPQ4(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]);
      MPMSHPQ4();
      ~MPMSHPQ4();
      double N1;                                                                // Shape Function at Node 1
      double N2;                                                                // Shape Function at Node 2
      double N3;                                                                // Shape Function at Node 3
      double N4;                                                                // Shape Function at Node 4

      double dN1dX;                                                             // Shape Function at Node 1 differentiated w.r.t X
      double dN2dX;                                                             // Shape Function at Node 2 differentiated w.r.t X
      double dN3dX;                                                             // Shape Function at Node 3 differentiated w.r.t X
      double dN4dX;                                                             // Shape Function at Node 4 differentiated w.r.t X

      double dN1dY;                                                             // Shape Function at Node 1 differentiated w.r.t Y
      double dN2dY;                                                             // Shape Function at Node 2 differentiated w.r.t Y
      double dN3dY;                                                             // Shape Function at Node 3 differentiated w.r.t Y
      double dN4dY;                                                             // Shape Function at Node 4 differentiated w.r.t Y

      double *N[4] = {&N1, &N2, &N3, &N4};                                      // Introduce vector structure to loop over ansatz functions ! use *() to get value
      double *dNdX[4] = {&dN1dX, &dN2dX, &dN3dX, &dN4dX};                       // Introduce vector structure to loop over ansatz functions ! use *() to get value
      double *dNdY[4] = {&dN1dY, &dN2dY, &dN3dY, &dN4dY};                       // Introduce vector structure to loop over ansatz functions ! use *() to get value

      void evaluate(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]);
      //double getMass(void);                       // A Member Function that returns the Current Mass of the Particle  MemberFunction: hass access to all data of the object
      //void Report(void);                          // A Member Function to print out a report of this object
      //std::cout << MyShape.N1 << std::endl;
      //std::cout << *MyShape.SHP[0] << std::endl;

  private:
    double SHPN1(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]);
    double SHPN2(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]);
    double SHPN3(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]);
    double SHPN4(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]);

    double dSHPN1dX(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]);
    double dSHPN2dX(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]);
    double dSHPN3dX(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]);
    double dSHPN4dX(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]);

    double dSHPN1dY(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]);
    double dSHPN2dY(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]);
    double dSHPN3dY(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]);
    double dSHPN4dY(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]);

    void SHPTest(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]);
};

MPMSHPQ4::~MPMSHPQ4(){}
MPMSHPQ4::MPMSHPQ4(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]){
  bool DetailedOutput = true;
  N1 = SHPN1(X1, X2, X3, X4, XP);
  N2 = SHPN2(X1, X2, X3, X4, XP);
  N3 = SHPN3(X1, X2, X3, X4, XP);
  N4 = SHPN4(X1, X2, X3, X4, XP);

  dN1dX = dSHPN1dX(X1, X2, X3, X4, XP);
  dN2dX = dSHPN2dX(X1, X2, X3, X4, XP);
  dN3dX = dSHPN3dX(X1, X2, X3, X4, XP);
  dN4dX = dSHPN4dX(X1, X2, X3, X4, XP);

  dN1dY = dSHPN1dY(X1, X2, X3, X4, XP);
  dN2dY = dSHPN2dY(X1, X2, X3, X4, XP);
  dN3dY = dSHPN3dY(X1, X2, X3, X4, XP);
  dN4dY = dSHPN4dY(X1, X2, X3, X4, XP);

  if (DetailedOutput) SHPTest(X1, X2, X3, X4, XP);
}
MPMSHPQ4::MPMSHPQ4(){
}

// Method to evaluate the shape functions
void MPMSHPQ4::evaluate(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]){
  N1 = SHPN1(X1, X2, X3, X4, XP);
  N2 = SHPN2(X1, X2, X3, X4, XP);
  N3 = SHPN3(X1, X2, X3, X4, XP);
  N4 = SHPN4(X1, X2, X3, X4, XP);

  dN1dX = dSHPN1dX(X1, X2, X3, X4, XP);
  dN2dX = dSHPN2dX(X1, X2, X3, X4, XP);
  dN3dX = dSHPN3dX(X1, X2, X3, X4, XP);
  dN4dX = dSHPN4dX(X1, X2, X3, X4, XP);

  dN1dY = dSHPN1dY(X1, X2, X3, X4, XP);
  dN2dY = dSHPN2dY(X1, X2, X3, X4, XP);
  dN3dY = dSHPN3dY(X1, X2, X3, X4, XP);
  dN4dY = dSHPN4dY(X1, X2, X3, X4, XP);
}

// Shape Function For Node 1 Evaluated at XP
double MPMSHPQ4::SHPN1(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]){
  return (X2[0]*(-(XP[0]*(X3[1] - X4[1])*(X2[1] - XP[1])) + X4[0]*(X2[1] - X4[1])*(X3[1] - XP[1]) - X3[0]*(X2[1] - X3[1])*(X4[1] - XP[1])) -
     X3[0]*X4[0]*(X3[1] - X4[1])*(X2[1] - XP[1]) + X3[0]*XP[0]*(X2[1] - X4[1])*(X3[1] - XP[1]) - X4[0]*XP[0]*(X2[1] - X3[1])*(X4[1] - XP[1]))/
   (X1[0]*(X4[0]*(X2[1] - X3[1])*(X1[1] - X4[1]) - X3[0]*(X1[1] - X3[1])*(X2[1] - X4[1]) + X2[0]*(X1[1] - X2[1])*(X3[1] - X4[1])) +
     X2[0]*X3[0]*(X2[1] - X3[1])*(X1[1] - X4[1]) - X2[0]*X4[0]*(X1[1] - X3[1])*(X2[1] - X4[1]) + X3[0]*X4[0]*(X1[1] - X2[1])*(X3[1] - X4[1]));
}
// Shape Function For Node 1 Evaluated at XP
double MPMSHPQ4::SHPN2(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]){
  return (X1[0]*(XP[0]*(X3[1] - X4[1])*(X1[1] - XP[1]) - X4[0]*(X1[1] - X4[1])*(X3[1] - XP[1]) + X3[0]*(X1[1] - X3[1])*(X4[1] - XP[1])) +
     X3[0]*X4[0]*(X3[1] - X4[1])*(X1[1] - XP[1]) - X3[0]*XP[0]*(X1[1] - X4[1])*(X3[1] - XP[1]) + X4[0]*XP[0]*(X1[1] - X3[1])*(X4[1] - XP[1]))/
   (X1[0]*(X4[0]*(X2[1] - X3[1])*(X1[1] - X4[1]) - X3[0]*(X1[1] - X3[1])*(X2[1] - X4[1]) + X2[0]*(X1[1] - X2[1])*(X3[1] - X4[1])) +
     X2[0]*X3[0]*(X2[1] - X3[1])*(X1[1] - X4[1]) - X2[0]*X4[0]*(X1[1] - X3[1])*(X2[1] - X4[1]) + X3[0]*X4[0]*(X1[1] - X2[1])*(X3[1] - X4[1]));
}
// Shape Function For Node 1 Evaluated at XP
double MPMSHPQ4::SHPN3(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]){
  return (X1[0]*(-(XP[0]*(X2[1] - X4[1])*(X1[1] - XP[1])) + X4[0]*(X1[1] - X4[1])*(X2[1] - XP[1]) - X2[0]*(X1[1] - X2[1])*(X4[1] - XP[1])) -
     X2[0]*X4[0]*(X2[1] - X4[1])*(X1[1] - XP[1]) + X2[0]*XP[0]*(X1[1] - X4[1])*(X2[1] - XP[1]) - X4[0]*XP[0]*(X1[1] - X2[1])*(X4[1] - XP[1]))/
   (X1[0]*(X4[0]*(X2[1] - X3[1])*(X1[1] - X4[1]) - X3[0]*(X1[1] - X3[1])*(X2[1] - X4[1]) + X2[0]*(X1[1] - X2[1])*(X3[1] - X4[1])) +
     X2[0]*X3[0]*(X2[1] - X3[1])*(X1[1] - X4[1]) - X2[0]*X4[0]*(X1[1] - X3[1])*(X2[1] - X4[1]) + X3[0]*X4[0]*(X1[1] - X2[1])*(X3[1] - X4[1]));
}
// Shape Function For Node 1 Evaluated at XP
double MPMSHPQ4::SHPN4(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]){
  return (X1[0]*(XP[0]*(X2[1] - X3[1])*(X1[1] - XP[1]) - X3[0]*(X1[1] - X3[1])*(X2[1] - XP[1]) + X2[0]*(X1[1] - X2[1])*(X3[1] - XP[1])) +
     X2[0]*X3[0]*(X2[1] - X3[1])*(X1[1] - XP[1]) - X2[0]*XP[0]*(X1[1] - X3[1])*(X2[1] - XP[1]) + X3[0]*XP[0]*(X1[1] - X2[1])*(X3[1] - XP[1]))/
   (X1[0]*(X4[0]*(X2[1] - X3[1])*(X1[1] - X4[1]) - X3[0]*(X1[1] - X3[1])*(X2[1] - X4[1]) + X2[0]*(X1[1] - X2[1])*(X3[1] - X4[1])) +
     X2[0]*X3[0]*(X2[1] - X3[1])*(X1[1] - X4[1]) - X2[0]*X4[0]*(X1[1] - X3[1])*(X2[1] - X4[1]) + X3[0]*X4[0]*(X1[1] - X2[1])*(X3[1] - X4[1]));
}

// Derived Shape Functions
double MPMSHPQ4::dSHPN1dX(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]){
  return (-(X2[0]*(X3[1] - X4[1])*(X2[1] - XP[1])) + X3[0]*(X2[1] - X4[1])*(X3[1] - XP[1]) - X4[0]*(X2[1] - X3[1])*(X4[1] - XP[1]))/(X1[0]*(X4[0]*(X2[1] - X3[1])*(X1[1] - X4[1]) - X3[0]*(X1[1] - X3[1])*(X2[1] - X4[1]) + X2[0]*(X1[1] - X2[1])*(X3[1] - X4[1])) + X2[0]*X3[0]*(X2[1] - X3[1])*(X1[1] - X4[1]) - X2[0]*X4[0]*(X1[1] - X3[1])*(X2[1] - X4[1]) + X3[0]*X4[0]*(X1[1] - X2[1])*(X3[1] - X4[1]));
}
double MPMSHPQ4::dSHPN2dX(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]){
  return (X1[0]*(X3[1] - X4[1])*(X1[1] - XP[1]) - X3[0]*(X1[1] - X4[1])*(X3[1] - XP[1]) + X4[0]*(X1[1] - X3[1])*(X4[1] - XP[1]))/(X1[0]*(X4[0]*(X2[1] - X3[1])*(X1[1] - X4[1]) - X3[0]*(X1[1] - X3[1])*(X2[1] - X4[1]) + X2[0]*(X1[1] - X2[1])*(X3[1] - X4[1])) + X2[0]*X3[0]*(X2[1] - X3[1])*(X1[1] - X4[1]) - X2[0]*X4[0]*(X1[1] - X3[1])*(X2[1] - X4[1]) + X3[0]*X4[0]*(X1[1] - X2[1])*(X3[1] - X4[1]));
}
double MPMSHPQ4::dSHPN3dX(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]){
  return (-(X1[0]*(X2[1] - X4[1])*(X1[1] - XP[1])) + X2[0]*(X1[1] - X4[1])*(X2[1] - XP[1]) - X4[0]*(X1[1] - X2[1])*(X4[1] - XP[1]))/(X1[0]*(X4[0]*(X2[1] - X3[1])*(X1[1] - X4[1]) - X3[0]*(X1[1] - X3[1])*(X2[1] - X4[1]) + X2[0]*(X1[1] - X2[1])*(X3[1] - X4[1])) + X2[0]*X3[0]*(X2[1] - X3[1])*(X1[1] - X4[1]) - X2[0]*X4[0]*(X1[1] - X3[1])*(X2[1] - X4[1]) + X3[0]*X4[0]*(X1[1] - X2[1])*(X3[1] - X4[1]));
}
double MPMSHPQ4::dSHPN4dX(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]){
  return (X1[0]*(X2[1] - X3[1])*(X1[1] - XP[1]) - X2[0]*(X1[1] - X3[1])*(X2[1] - XP[1]) + X3[0]*(X1[1] - X2[1])*(X3[1] - XP[1]))/(X1[0]*(X4[0]*(X2[1] - X3[1])*(X1[1] - X4[1]) - X3[0]*(X1[1] - X3[1])*(X2[1] - X4[1]) + X2[0]*(X1[1] - X2[1])*(X3[1] - X4[1])) + X2[0]*X3[0]*(X2[1] - X3[1])*(X1[1] - X4[1]) - X2[0]*X4[0]*(X1[1] - X3[1])*(X2[1] - X4[1]) + X3[0]*X4[0]*(X1[1] - X2[1])*(X3[1] - X4[1]));
}
double MPMSHPQ4::dSHPN1dY(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]){
  return (X4[0]*XP[0]*(X2[1] - X3[1]) + X2[0]*(X3[0]*(X2[1] - X3[1]) - X4[0]*(X2[1] - X4[1]) + XP[0]*(X3[1] - X4[1])) - X3[0]*XP[0]*(X2[1] - X4[1]) + X3[0]*X4[0]*(X3[1] - X4[1]))/(X1[0]*(X4[0]*(X2[1] - X3[1])*(X1[1] - X4[1]) - X3[0]*(X1[1] - X3[1])*(X2[1] - X4[1]) + X2[0]*(X1[1] - X2[1])*(X3[1] - X4[1])) + X2[0]*X3[0]*(X2[1] - X3[1])*(X1[1] - X4[1]) - X2[0]*X4[0]*(X1[1] - X3[1])*(X2[1] - X4[1]) + X3[0]*X4[0]*(X1[1] - X2[1])*(X3[1] - X4[1]));
}
double MPMSHPQ4::dSHPN2dY(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]){
  return (-(X4[0]*XP[0]*(X1[1] - X3[1])) + X1[0]*(-(X3[0]*(X1[1] - X3[1])) + X4[0]*(X1[1] - X4[1]) - XP[0]*(X3[1] - X4[1])) + X3[0]*XP[0]*(X1[1] - X4[1]) - X3[0]*X4[0]*(X3[1] - X4[1]))/(X1[0]*(X4[0]*(X2[1] - X3[1])*(X1[1] - X4[1]) - X3[0]*(X1[1] - X3[1])*(X2[1] - X4[1]) + X2[0]*(X1[1] - X2[1])*(X3[1] - X4[1])) + X2[0]*X3[0]*(X2[1] - X3[1])*(X1[1] - X4[1]) - X2[0]*X4[0]*(X1[1] - X3[1])*(X2[1] - X4[1]) + X3[0]*X4[0]*(X1[1] - X2[1])*(X3[1] - X4[1]));
}
double MPMSHPQ4::dSHPN3dY(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]){
  return (X4[0]*XP[0]*(X1[1] - X2[1]) + X1[0]*(X2[0]*(X1[1] - X2[1]) - X4[0]*(X1[1] - X4[1]) + XP[0]*(X2[1] - X4[1])) - X2[0]*XP[0]*(X1[1] - X4[1]) + X2[0]*X4[0]*(X2[1] - X4[1]))/(X1[0]*(X4[0]*(X2[1] - X3[1])*(X1[1] - X4[1]) - X3[0]*(X1[1] - X3[1])*(X2[1] - X4[1]) + X2[0]*(X1[1] - X2[1])*(X3[1] - X4[1])) + X2[0]*X3[0]*(X2[1] - X3[1])*(X1[1] - X4[1]) - X2[0]*X4[0]*(X1[1] - X3[1])*(X2[1] - X4[1]) + X3[0]*X4[0]*(X1[1] - X2[1])*(X3[1] - X4[1]));
}
double MPMSHPQ4::dSHPN4dY(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]){
  return (-(X3[0]*XP[0]*(X1[1] - X2[1])) + X1[0]*(-(X2[0]*(X1[1] - X2[1])) + X3[0]*(X1[1] - X3[1]) - XP[0]*(X2[1] - X3[1])) + X2[0]*XP[0]*(X1[1] - X3[1]) - X2[0]*X3[0]*(X2[1] - X3[1]))/(X1[0]*(X4[0]*(X2[1] - X3[1])*(X1[1] - X4[1]) - X3[0]*(X1[1] - X3[1])*(X2[1] - X4[1]) + X2[0]*(X1[1] - X2[1])*(X3[1] - X4[1])) + X2[0]*X3[0]*(X2[1] - X3[1])*(X1[1] - X4[1]) - X2[0]*X4[0]*(X1[1] - X3[1])*(X2[1] - X4[1]) + X3[0]*X4[0]*(X1[1] - X2[1])*(X3[1] - X4[1]));
}

// Shape Function Test
void MPMSHPQ4::SHPTest(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]){
  std::cout << "Shape Function Created! Test: :" << std::endl;
  std::cout << "X1 :" << X1[0] << ", " << X1[1] << ", " << X1[2] << std::endl;
  std::cout << "X2 :" << X2[0] << ", " << X2[1] << ", " << X2[2] << std::endl;
  std::cout << "X3 :" << X3[0] << ", " << X3[1] << ", " << X3[2] << std::endl;
  std::cout << "X4 :" << X4[0] << ", " << X4[1] << ", " << X4[2] << std::endl;
  std::cout << "XP :" << XP[0] << ", " << XP[1] << ", " << XP[2] << std::endl;

  //Test for Partition OF Unity
  std::cout << "N(X1) :" << SHPN1(X1, X2, X3, X4, X1) << ", " << SHPN2(X1, X2, X3, X4, X1) << ", " << SHPN3(X1, X2, X3, X4, X1) << ", " << SHPN4(X1, X2, X3, X4, X1) << std::endl;
  std::cout << "N(X2) :" << SHPN1(X1, X2, X3, X4, X2) << ", " << SHPN2(X1, X2, X3, X4, X2) << ", " << SHPN3(X1, X2, X3, X4, X2) << ", " << SHPN4(X1, X2, X3, X4, X2) << std::endl;
  std::cout << "N(X3) :" << SHPN1(X1, X2, X3, X4, X3) << ", " << SHPN2(X1, X2, X3, X4, X3) << ", " << SHPN3(X1, X2, X3, X4, X3) << ", " << SHPN4(X1, X2, X3, X4, X3) << std::endl;
  std::cout << "N(X4) :" << SHPN1(X1, X2, X3, X4, X4) << ", " << SHPN2(X1, X2, X3, X4, X4) << ", " << SHPN3(X1, X2, X3, X4, X4) << ", " << SHPN4(X1, X2, X3, X4, X4) << std::endl;

  //Test Call of dx
  std::cout << "dNdx(X1) :" << dSHPN1dX(X1, X2, X3, X4, X1) << ", " << dSHPN2dX(X1, X2, X3, X4, X1) << ", " << dSHPN3dX(X1, X2, X3, X4, X1) << ", " << dSHPN4dX(X1, X2, X3, X4, X1) << std::endl;
  std::cout << "dNdx(X2) :" << dSHPN1dX(X1, X2, X3, X4, X2) << ", " << dSHPN2dX(X1, X2, X3, X4, X2) << ", " << dSHPN3dX(X1, X2, X3, X4, X2) << ", " << dSHPN4dX(X1, X2, X3, X4, X2) << std::endl;
  std::cout << "dNdx(X3) :" << dSHPN1dX(X1, X2, X3, X4, X3) << ", " << dSHPN2dX(X1, X2, X3, X4, X3) << ", " << dSHPN3dX(X1, X2, X3, X4, X3) << ", " << dSHPN4dX(X1, X2, X3, X4, X3) << std::endl;
  std::cout << "dNdx(X4) :" << dSHPN1dX(X1, X2, X3, X4, X4) << ", " << dSHPN2dX(X1, X2, X3, X4, X4) << ", " << dSHPN3dX(X1, X2, X3, X4, X4) << ", " << dSHPN4dX(X1, X2, X3, X4, X4) << std::endl;

  //Test Call of dx
  std::cout << "dNdy(X1) :" << dSHPN1dY(X1, X2, X3, X4, X1) << ", " << dSHPN2dY(X1, X2, X3, X4, X1) << ", " << dSHPN3dY(X1, X2, X3, X4, X1) << ", " << dSHPN4dY(X1, X2, X3, X4, X1) << std::endl;
  std::cout << "dNdy(X2) :" << dSHPN1dY(X1, X2, X3, X4, X2) << ", " << dSHPN2dY(X1, X2, X3, X4, X2) << ", " << dSHPN3dY(X1, X2, X3, X4, X2) << ", " << dSHPN4dY(X1, X2, X3, X4, X2) << std::endl;
  std::cout << "dNdy(X3) :" << dSHPN1dY(X1, X2, X3, X4, X3) << ", " << dSHPN2dY(X1, X2, X3, X4, X3) << ", " << dSHPN3dY(X1, X2, X3, X4, X3) << ", " << dSHPN4dY(X1, X2, X3, X4, X3) << std::endl;
  std::cout << "dNdy(X4) :" << dSHPN1dY(X1, X2, X3, X4, X4) << ", " << dSHPN2dY(X1, X2, X3, X4, X4) << ", " << dSHPN3dY(X1, X2, X3, X4, X4) << ", " << dSHPN4dY(X1, X2, X3, X4, X4) << std::endl;

  //Test For Interpolation of a constant field
  // interpolate nodal values TI over element defined with TestXI and evalueate the interpolation as well as its gradient at points TestXP
  double TI[4] = {2.0,2.0,2.0,2.0};
  double TestXI[4][2] = {{0.0 ,0.0},{1.0,0.0},{1.0,1.0},{0.0,1.0}};
  double TestXP[4][2] = {{0.2 ,0.2},{0.8,0.2},{0.8,0.8},{0.2,0.8}};
  double TXP[4];
  double Grad_TXP[4][2];
  for(int TestPoint=0;TestPoint<4;TestPoint++){
    evaluate(TestXI[0], TestXI[1], TestXI[2], TestXI[3], TestXP[TestPoint]);
    // interpolate TI
    TXP[TestPoint] = 0; Grad_TXP[TestPoint][0] = 0.0; Grad_TXP[TestPoint][1] = 0.0;
    for(int i=0;i<4;i++) TXP[TestPoint] += (*N[i]) * TI[i];
    for(int i=0;i<4;i++) Grad_TXP[TestPoint][0] += (*dNdX[i]) * TI[i];
    for(int i=0;i<4;i++) Grad_TXP[TestPoint][1] += (*dNdY[i]) * TI[i];
  }
  std::cout << "Constant Field TI = 2 at Test Points: "  << std::endl;
  std::cout << "           P1     P2     P3     P4 " << std::endl;
  std::cout << " TI     " << std::setw(5) << TXP[0] << " ," << std::setw(5) << TXP[1] << " ," << std::setw(5) << TXP[2] << " ," << std::setw(5) << TXP[3] << std::endl;
  std::cout << " dTidX  " << std::setw(5) << Grad_TXP[0][0] << " ," << std::setw(5) << Grad_TXP[1][0] << " ," << std::setw(5) << Grad_TXP[2][0] << " ," << std::setw(5) << Grad_TXP[3][0] << std::endl;
  std::cout << " dTidX  " << std::setw(5) << Grad_TXP[0][1] << " ," << std::setw(5) << Grad_TXP[1][1] << " ," << std::setw(5) << Grad_TXP[2][1] << " ," << std::setw(5) << Grad_TXP[3][1] << std::endl;

  //Test For Interpolation of a bilinear field
  // interpolate nodal values TI over element defined with TestXI and evalueate the interpolation as well as its gradient at points TestXP
  double TIBiLin[4] = {-2.0,3.0,5.0,7.239};
  double TestXIBiLin[4][2] = {{0.0 ,0.0},{1.0,0.0},{1.0,1.0},{0.0,1.0}};
  double TestXPBiLin[4][2] = {{0.2 ,0.2},{0.8,0.2},{0.8,0.8},{0.2,0.8}};
  double TXPBiLin[4];
  double Grad_TXPBiLin[4][2];
  for(int TestPoint=0;TestPoint<4;TestPoint++){
    evaluate(TestXIBiLin[0], TestXIBiLin[1], TestXIBiLin[2], TestXIBiLin[3], TestXPBiLin[TestPoint]);
    // interpolate TI
    TXPBiLin[TestPoint] = 0; Grad_TXPBiLin[TestPoint][0] = 0.0; Grad_TXPBiLin[TestPoint][1] = 0.0;
    for(int i=0;i<4;i++) TXPBiLin[TestPoint] += (*N[i]) * TIBiLin[i];
    for(int i=0;i<4;i++) Grad_TXPBiLin[TestPoint][0] += (*dNdX[i]) * TIBiLin[i];
    for(int i=0;i<4;i++) Grad_TXPBiLin[TestPoint][1] += (*dNdY[i]) * TIBiLin[i];
  }
  std::cout << "Constant Field T(x,y) = -2+x(5-7.239y)+ 9.239y at Test Points: "  << std::endl;
  std::cout << "           P1     P2     P3     P4 " << std::endl;
  std::cout << " TI     " << std::setw(5) << TXPBiLin[0] << " ," << std::setw(5) << TXPBiLin[1] << " ," << std::setw(5) << TXPBiLin[2] << " ," << std::setw(5) << TXPBiLin[3] << std::endl;
  std::cout << " dTidX  " << std::setw(5) << Grad_TXPBiLin[0][0] << " ," << std::setw(5) << Grad_TXPBiLin[1][0] << " ," << std::setw(5) << Grad_TXPBiLin[2][0] << " ," << std::setw(5) << Grad_TXPBiLin[3][0] << std::endl;
  std::cout << " dTidX  " << std::setw(5) << Grad_TXPBiLin[0][1] << " ," << std::setw(5) << Grad_TXPBiLin[1][1] << " ," << std::setw(5) << Grad_TXPBiLin[2][1] << " ," << std::setw(5) << Grad_TXPBiLin[3][1] << std::endl;

}
