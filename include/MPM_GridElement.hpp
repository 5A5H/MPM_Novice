class MPMGridElement {
  public:
      MPMGridElement();
      MPMGridElement(int id, int n1, int n2, int n3, int n4);
      ~MPMGridElement();
      int ID;
      int N1;
      int N2;
      int N3;
      int N4;
      int *N[4] = {&N1, &N2, &N3, &N4};
      void Report();

  private:
};
