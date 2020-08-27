#include <cstdlib>
#include <iostream>
#include <TH1D.h>
#include <TRandom.h>
using namespace std;
void jonathon(){
  TH1* h1 = new TH1D("h1", "h1 title", 100.0, 0.0, 100000000004.0);
  for (int i=0; i<100; i++){
    h1->Fill(    rand(   )  );
    cout << rand(   ) << endl;
  }
  h1-> Draw();

}
