#include "PackSet.h"
#include "Random.hpp"
#include <iostream>
#include <cmath>

using namespace std;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    exit(EXIT_FAILURE);
  }
  string temp(argv[1]);
  size_t num=strtol(temp.c_str(), NULL, 10);
  Random::Global.Reseed();
  cout << "Seed set to: " << Random::Global.SeedString() << endl;
  PackSet s;
  while (s.Number() < num)
    s.Add(Quaternion(Random::Global.Normal<double>(),
		     Random::Global.Normal<double>(),
		     Random::Global.Normal<double>(),
		     Random::Global.Normal<double>()));
  for (size_t iter = 0; iter < 200; ++iter) {
    for (size_t i = 0; i < 2 * num; ++num)
      s.MinMaxRadius(Random::Global(s.Number()), 0.01);
    double d = s.MaxRadiusA(0.01);
    cout << iter << " "
	 << d << endl;
  }
  return 0;
}
