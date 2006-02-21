#include "PackSet.h"
#include "Random.hpp"
#include <iostream>
#include <iomanip>
//#include <cstdlib>

void PackSet::Add(const Quaternion& q, unsigned char cat, bool skipcheck) {
  Quaternion qq = m_pre;
  qq *= q;
  qq *= m_post;
  qq.Normalize();
  qq.Canonicalize();
  // Eliminate duplicates and one of q/-q.
  size_t k = FindClosest(qq);
  if (k == m_set.size() || skipcheck
      || Closeness(qq, m_set[k]) < 1 - 1.0e-9) {
    m_set.push_back(qq);
    m_cat.push_back(cat);
    m_maxcat = std::max(++cat, m_maxcat);
  }
}

void PackSet::Analyze(size_t num) const {
  vector<size_t> count(m_set.size(), 0);
  vector<double> close(m_set.size(), 1);
  for (size_t i = 0; i < num; ++i) {
    Quaternion q(Random::Global.Normal<double>(),
		 Random::Global.Normal<double>(),
		 Random::Global.Normal<double>(),
		 Random::Global.Normal<double>());
    q.Normalize();
    size_t j = FindClosest(q);
    count[j]++;
    double c = Closeness(q, m_set[j]);
    if (c < close[j])
      close[j] = c;
  }
  cout << "Analysis with " << num << " probes" << endl;
  cout << "Point Neighbor Dist Count Radius" << endl;
  double mindist = 6, maxdist = 0, minradius = 6, maxradius = 0;
  size_t mincount = num, maxcount = 0;
  for (size_t i = 0; i < m_set.size(); ++i) {
    size_t j = FindClosest(i);
    double c = Closeness(m_set[i], m_set[j]);
    double d = Dist(c);
    double r = Dist(close[i]);
    cout << i << " "
	 << j << " "
	 << d << " "
	 << count[i] << " "
	 << r << endl;
    if (d < mindist)
      mindist = d;
    if (d > maxdist)
      maxdist = d;
    if (r < minradius)
      minradius = r;
    if (r > maxradius)
      maxradius = r;
    if (count[i] < mincount)
      mincount = count[i];
    if (count[i] > maxcount)
      maxcount = count[i];
  }
  cout << "Mindist = " << mindist
       << " Maxdist = " << maxdist
       << " Minradius = " << minradius
       << " Maxradius = " << maxradius
       << " Mincount = " << mincount
       << " Maxcount = " << maxcount << endl;
}

void PackSet::Analyze1() const {
  vector<unsigned long long> sumcount(m_maxcat, 0);
  unsigned long long num = 0;
  vector<double> cnt(m_maxcat, 0);
  for (size_t i = 0; i < m_set.size(); ++i)
    cnt[m_cat[i]]++;
  while (true) {
    num++;
    Quaternion q(Random::Global.Normal<double>(),
		 Random::Global.Normal<double>(),
		 Random::Global.Normal<double>(),
		 Random::Global.Normal<double>());
    q.Normalize();
    size_t j = FindClosest(q);
    sumcount[m_cat[j]]++;
    if (num % 1000000 == 0) {
      cout << setw(11) << num;
      for (size_t i = 0; i < m_maxcat; ++i)
	cout << " " << setw(11) << sumcount[i];
      for (size_t i = 0; i < m_maxcat; ++i)
	cout << " " << setw(8) << 1000.0*double(sumcount[i])/double(num)/cnt[i];
      cout << endl;
    }
  }
}

void PackSet::Analyze0(size_t num) const {
  vector<double> count(m_set.size(), 0);
  vector<double> close(m_set.size(), 1);
  double d = 1.0/num;
  int n = num;
  double mult = 32;
  for (int i = 0; i < n; ++i) {
    double x = (i + 0.5) * d;
      for (int j = 0; j < n; ++j) {
	double y = (j + 0.5) * d;
	for (int k = 0; k < n; ++k) {
	  double z = (k + 0.5) * d;
	  Quaternion q(1.0, x, y, z);
	  double s = q.Magnitude();
	  s = 1/s;
	  s *= s;
	  s *= s;
	  q.Normalize();
	  for (int m = 0; m < 1; ++m) {
	    size_t j = FindClosest(q);
	    count[j] += s;
	    double c = Closeness(q, m_set[j]);
	    if (c < close[j])
	      close[j] = c;
	    q.CircularRotate(1);
	  }
	}
      }
  }
  cout << "Analysis with " << num << "^3 probes" << endl;
  if (0)
    cout << "Point Cat Neighbor Dist Count Radius" << endl;
  vector<double> mindist(m_maxcat, 6);
  vector<double> maxdist(m_maxcat, 0);
  vector<double> maxradius(m_maxcat, 0);
  vector<double> sumcount(m_maxcat, 0);
  vector<double> cnt(m_maxcat, 0);
  for (size_t i = 0; i < m_set.size(); ++i) {
    size_t j = FindClosest(i);
    double c = Closeness(m_set[i], m_set[j]);
    double d = Dist(c);
    double r = Dist(close[i]);
    size_t cat = m_cat[i];
    cnt[cat]++;
    if (0)
      cout << i << " "
	   << cat << " "
	   << j << " "
	   << d << " "
	   << count[i] << " "
	   << r << endl;
    if (d < mindist[cat])
      mindist[cat] = d;
    if (d > maxdist[cat])
      maxdist[cat] = d;
    if (r > maxradius[cat])
      maxradius[cat] = r;
    sumcount[cat] += count[i];
  }
  cout << "Cat count mindist maxdist maxradius weight totw" << endl;
  for (size_t i = 0; i < m_maxcat; ++i) {
    double weight = mult*sumcount[i]/cnt[i]*d*d*d;
    cout << i << " "
	 << cnt[i] << " "
	 << mindist[i] << " "
	 << maxdist[i] << " "
	 << maxradius[i] << " "
	 << weight << " "
	 << weight * cnt[i] << endl;
  }
}

double PackSet::Closeness(const Quaternion& q1, const Quaternion& q2) {
  return abs(q1.DotProduct(q2));
}

double PackSet::Dist(double closeness) {
  return 2 * acos(closeness);
}

size_t PackSet::FindClosest(size_t test) const {
  double x = -1;
  size_t res = m_set.size();
  for (size_t i = 0; i < m_set.size(); ++i) {
    if (i == test)
      continue;
    double y = Closeness(m_set[i], m_set[test]);
    if (y > x) {
      x = y;
      res = i;
    }
  }
  return res;
}

size_t PackSet::FindClosest(const Quaternion& q) const {
  double x = -1;
  size_t res = m_set.size();
  for (size_t i = 0; i < m_set.size(); ++i) {
    double y = Closeness(m_set[i], q);
    if (y > x) {
      x = y;
      res = i;
    }
  }
  return res;
}
  
double PackSet::MaxRadiusA(size_t k, double eps) const {
  double s0 = MinDistance(k) * 0.5;
  double res = 0;
  for (size_t n = 0; n < 12; ++n) {
    double s = s0;
    Quaternion  qq(1,
		   s * Random::Global.Normal<double>(),
		   s * Random::Global.Normal<double>(),
		   s * Random::Global.Normal<double>());
    qq *= m_set[k];
    qq.Normalize();
    double d = Dist(Closeness(qq, m_set[FindClosest(qq)]));
    double saved = d;
    s *= 0.5;
    size_t iter = 0;
    while (true) {
      iter++;
      Quaternion nq(1,
		    s * Random::Global.Normal<double>(),
		    s * Random::Global.Normal<double>(),
		    s * Random::Global.Normal<double>());
      nq *= qq;
      nq.Normalize();
      double nd = Dist(Closeness(nq, m_set[FindClosest(nq)]));
      if (nd > d) {
	qq = nq;
	d = nd;
      }
      if (iter % 20 == 0) {
	if (d <= saved + eps)
	  break;
	saved = d;
	s *= 0.9;
      }
    }
    if (d > res)
      res = d;
  }
  return res;
}

double PackSet::MaxRadiusA(double eps) const {
  double res = 0;
  size_t num = Number();
  for (size_t i = 0; i < num; ++i)
    res = max(res, MaxRadiusA(i, eps));
  return res;
}

double PackSet::MaxRadius(const Quaternion& q, double d, size_t num) const {
  Quaternion qq(q);
  qq.Normalize();
  size_t j = FindClosest(qq);
  double c = Closeness(qq, m_set[j]);
  while (num > 0) {
    num--;
    Quaternion nq(1,
		  d * Random::Global.Normal<double>(),
		  d * Random::Global.Normal<double>(),
		  d * Random::Global.Normal<double>());
    nq *= qq;
    nq.Normalize();
    size_t nj = FindClosest(nq);
    double nc = Closeness(nq, m_set[nj]);
    if (nc < c) {
      qq = nq;
      c = nc;
      j = nj;
    }
    if (num % 100 == 0) {
      d *= 0.95;
      if (false && num % 100000 == 0)
	cout << num << " " << setprecision(15) << Dist(c) << endl;
    }
  }
  return Dist(c);
}

double PackSet::MinDistance() const {
  size_t n = Number();
  double d = 10;
  for (size_t i = 0; i < n; ++i) {
    double d1 = MinDistance(i);
    if (d1 < d)
      d = d1;
  }
  return d;
}

double PackSet::MinDistance(size_t k) const {
  return Dist(Closeness(m_set[FindClosest(k)], m_set[k]));
}

size_t PackSet::MonteCarlo(size_t num, double delta, double beta) {
  size_t count = 0;
  size_t n = Number();
  for (size_t i = 0; i < num; ++i) {
    size_t k = Random::Global(n);
    double d = MinDistance(k);
    Quaternion old(m_set[k]);
    m_set[k] *= Quaternion(1,
			   delta * Random::Global.Normal<double>(),
			   delta * Random::Global.Normal<double>(),
			   delta * Random::Global.Normal<double>());
    m_set[k].Normalize();
    double d1 = MinDistance(k);
    // Distance acts like -E
    if (d1 > d || Random::Global.Prob(exp(beta*(d1-d))))
      ++count;			// Accept
    else
      m_set[k] = old;		// Reject
  }
  return count;
}

size_t PackSet::MonteCarloA(size_t num, double delta, double beta) {
  size_t count = 0;
  size_t n = Number();
  for (size_t i = 0; i < num; ++i) {
    size_t k = Random::Global(n);
    double d = MaxRadiusA(k,0.01);
    Quaternion old(m_set[k]);
    m_set[k] *= Quaternion(1,
			   delta * Random::Global.Normal<double>(),
			   delta * Random::Global.Normal<double>(),
			   delta * Random::Global.Normal<double>());
    m_set[k].Normalize();
    double d1 = MaxRadiusA(k,0.01);
    // Distance acts like +E
    if (d1 < d || Random::Global.Prob(exp(-beta*(d1-d))))
      ++count;			// Accept
    else
      m_set[k] = old;		// Reject
  }
  return count;
}

double PackSet::MinMaxRadius(size_t k, double eps) {
  double s0 = MinDistance(k) * 0.5;
  double res = 0;
  vector<Quaternion> boundary;
  for (size_t n = 0; n < 12; ++n) {
    double s = s0;
    Quaternion  qq(1,
		   s * Random::Global.Normal<double>(),
		   s * Random::Global.Normal<double>(),
		   s * Random::Global.Normal<double>());
    qq *= m_set[k];
    qq.Normalize();
    double d = Dist(Closeness(qq, m_set[FindClosest(qq)]));
    double saved = d;
    s *= 0.5;
    size_t iter = 0;
    while (true) {
      iter++;
      Quaternion nq(1,
		    s * Random::Global.Normal<double>(),
		    s * Random::Global.Normal<double>(),
		    s * Random::Global.Normal<double>());
      nq *= qq;
      nq.Normalize();
      double nd = Dist(Closeness(nq, m_set[FindClosest(nq)]));
      if (nd > d) {
	qq = nq;
	d = nd;
      }
      if (iter % 20 == 0) {
	if (d <= saved + eps)
	  break;
	saved = d;
	s *= 0.9;
      }
    }
    boundary.push_back(qq);
    if (d > res)
      res = d;
  }
  double sep = 0;
  for (size_t i = 0; i < boundary.size(); ++i)
    sep = max(sep, Dist(Closeness(m_set[k], boundary[i])));
  double s = s0 * 0.1;
  for (size_t n = 0; n < 100; ++n) {
    Quaternion qq(1,
		   s * Random::Global.Normal<double>(),
		   s * Random::Global.Normal<double>(),
		   s * Random::Global.Normal<double>());
    qq *= m_set[k];
    qq.Normalize();
    double nsep = 0;
    for (size_t i = 0; i < boundary.size(); ++i)
      nsep = max(nsep, Dist(Closeness(qq, boundary[i])));
    if (nsep < sep) {
      m_set[k] = qq;
      sep = nsep;
    }
  }
  return sep;
}

double PackSet::Volume(size_t ind, double radius, size_t n,
		       double& maxrad, const Quaternion& tx) {
  PackSet temp;
  temp.Add(Quaternion(1,0,0,0));
  Quaternion q = Member(ind);
  Quaternion qc = q.Conjugate();
  for (size_t i = 0; i < Number(); ++i) {
    if (i == ind)
      continue;
    if (Dist(Closeness(q, Member(i))) <= 2 * radius)
      temp.Add(Member(i) * qc, 0, true);
  }
  // Integrate over sphere of radius radt by constructing cubical grid
  // and count cubes.  V = int WithinCell(x) dV.  Better to integrate
  // over the surface of the unit sphere and compute V = int l(Omega)^3
  // dOmega/3 where l(Omega) is the extent of the Voronoi cell in
  // direction Omega which can be found by the bisection method.
#if 1
  double radt = Quaternion::AngleToTurn(radius);
  double vol = 0;
  double maxturn = 0;
  size_t ndiv = 4 + int(ceil(2*log(double(2*n))/log(2.0)));
  for (size_t i = 0; i < n; ++i) {
    double x = (2.0 * i + 1.0 - n)/ n;
    for (size_t j = 0; j < n; ++j) {
      double y = (2.0 * j + 1.0 - n)/ n;
      for (size_t k = 0; k < 2; ++k) {
	// z = +/- 1
	double z = 2.0 * k - 1.0;
	// Convert to unit vector
	double s = sqrt(x*x + y*y + 1);
	double x0 = x/s, y0 = y/s, z0 = z/s;
	s = 1/(s*s*s);		// Metric factor for this element.
	for (size_t m = 0; m < 3; ++m) {
	  // Rotate through 3 axes
	  {
	    double t = x0;
	    x0 = y0;
	    y0 = z0;
	    z0 = t;
	  }
	  Vector3D<double> v = tx.transformPoint(Vector3D<double>(x0, y0, z0));
	  double g = 0.5, del = 0.25;
	  for (size_t it = 0; it < ndiv; ++it) {
	    Quaternion t(cbrt(g) * radt * v, true);
	    // Quaternion t(g * radt * v, true);
	    // If closest is nonzero move closer
	    g += temp.FindClosest(t) ? -del : del;
	    del /= 2.0;
	  }
	  maxturn = max(g, maxturn);
	  vol += s * g;
	  // vol += s * g * g * g;
	}
      }
    }
  }
  // mult by area element and cone factor
  vol *= 4.0*radt*radt*radt/(3.0*n*n);
  // Correct for systematic discretization error.
  vol /= 1 + 0.25/(n*n);
  // Further possible correction:
  // vol /= 1 - 1.0/(n*n*n);
  maxrad = Quaternion::TurnToAngle(cbrt(maxturn)*radt);
  return vol;
#else
  double radt = Quaternion::AngleToTurn(radius);
  size_t count = 0;
  double mincloseness = 2;
  for (size_t i = 0; i < n; ++i) {
    double x = (2.0 * i + 1.0 - n) * radt / n;
    for (size_t j = 0; j < n; ++j) {
      double y = (2.0 * j + 1.0 - n) * radt / n;
      if (x * x + y * y > radt * radt)
	continue;
      for (size_t k = 0; k < n; ++k) {
	double z = (2.0 * k + 1.0 - n) * radt / n;
	if (x * x + y * y + z * z > radt * radt)
	  continue;
	Quaternion t(tx.transformPoint(Vector3D<double>(x, y, z)), true);
	size_t q = temp.FindClosest(t);
	if (q == 0) {
	  count++;
	  mincloseness = min(mincloseness, abs(t.w()));
	}
      }
    }
  }
  radt *= 2.0/n;
  maxrad = Dist(mincloseness);
  //  cout << temp.Number() - 1 << " " << Dist(mincloseness)*180/M_PI << endl;
  return count * radt * radt * radt;
#endif
}
    
  
  
