#ifndef PACKSET_H
#define PACKSET_H

#include <vector>
#include "Quaternion.h"

using namespace std;

class PackSet {
public:
  PackSet(Quaternion pre = Quaternion(1,0,0,0),
	  Quaternion post = Quaternion(1,0,0,0))
    : m_pre(pre)
    , m_post(post)
    , m_maxcat(0) {
    m_pre.Normalize();
    m_post.Normalize();
  }
  void Add(const Quaternion& q, unsigned char cat = 0);
  size_t Number() const {
    return m_set.size();
  }
  void Analyze(size_t num) const;
  void Analyze0(size_t num) const;
  void Analyze1() const;
  double MaxRadius(const Quaternion& g, double d, size_t num) const;
  double MaxRadiusA(size_t k, double eps) const;
  double MaxRadiusA(double eps) const;
  double MinDistance() const;
  double MinDistance(size_t k) const;
  double MinMaxRadius(size_t k, double eps);
  size_t MonteCarlo(size_t num, double delta, double beta);
  size_t MonteCarloA(size_t num, double delta, double beta);
private:
  vector<Quaternion> m_set;
  vector<unsigned char> m_cat;
  Quaternion m_pre;
  Quaternion m_post;
  unsigned char m_maxcat;
  size_t FindClosest(size_t test) const;
  size_t FindClosest(const Quaternion& q) const;
  static double Closeness(const Quaternion& q1, const Quaternion& q2);
  static double Dist(double closeness);
};

#endif
