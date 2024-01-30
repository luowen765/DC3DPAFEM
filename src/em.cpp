
#include "em.h"

// Print time in seconds(sec), minutes(min), and hours(hr)
void print_time(double time_in_seconds) {
  int seconds = (int)time_in_seconds;
  int hr = seconds / 3600;
  seconds %= 3600;
  int min = seconds / 60;
  double sec = time_in_seconds - hr * 3600 - min * 60;
  std::cout << hr << " hr " << min << " min " << sec << " sec\n";
}

double length_two_point(double *coord1, double *coord2) {
  double length = (coord2[0] - coord1[0]) * (coord2[0] - coord1[0]) +
                  (coord2[1] - coord1[1]) * (coord2[1] - coord1[1]) +
                  (coord2[2] - coord1[2]) * (coord2[2] - coord1[2]);
  length = std::sqrt(length);
  return length;
}

double length_two_point(Vector &coord1, Vector &coord2) {

  myassert(coord1.Size() == 3);
  myassert(coord2.Size() == 3);
  double c1[3];
  double c2[3];
  for (int i = 0; i < 3; i++) {
    c1[i] = coord1[i];
    c2[i] = coord2[i];
  }
  double length = length_two_point(c1, c2);
  return length;
}

bool cmp(std::vector<double> a, std::vector<double> b) {
  if (a[0] != b[0])
    return a[0] < b[0];
  else if (a[1] != b[1])
    return a[1] < b[1];
  else
    return a[2] < b[2];
}

bool cmp2(std::vector<double> a, std::vector<double> b) {
  if (a[0] != b[0])
    return a[0] < b[0];
  else if (a[1] != b[1])
    return a[1] < b[1];
  else if (a[2] != b[2])
    return a[2] < b[2];
  else if (a[3] != b[3])
    return a[3] < b[3];
  else if (a[4] != b[4])
    return a[4] < b[4];
  else
    return a[5] < b[5];
}

void remove_duplicate_vertex(std::vector<mfem::Vertex> &all_point) {
  struct Point {
    double x;
    double y;
    double z;
    bool operator==(const Point &other) const {
      return (x == other.x) && (y == other.y) && (z == other.z);
    }
    bool operator<(const Point &other) const {
      if (x != other.x) {
        return x < other.x;
      }
      if (y != other.y) {
        return y < other.y;
      }
      return z < other.z;
    }
  };
  std::vector<Point> points;
  points.resize(all_point.size());
  for (int i = 0; i < all_point.size(); i++) {
    points[i].x = all_point[i](0);
    points[i].y = all_point[i](1);
    points[i].z = all_point[i](2);
  }

  std::sort(points.begin(), points.end());
  auto unique_end = std::unique(points.begin(), points.end());
  points.erase(unique_end, points.end());

  all_point.resize(points.size());
  for (int i = 0; i < points.size(); i++) {
    all_point[i](0) = points[i].x;
    all_point[i](1) = points[i].y;
    all_point[i](2) = points[i].z;
  }
}

double B_value(Vector p, Vector pi, DenseMatrix sigma) { // r*rho*r
  Vector r(3);
  r = 0.0;
  subtract(p, pi, r); 
  Vector temp(3);
  temp = 0.0;
  DenseMatrix rho = sigma;
  rho.Invert(); // attention！ Replaces the current matrix with its inverse

  rho.AddMult(r, temp);     
  double Bvalue = r * temp; // ri*rho*ri
  return Bvalue;
}

double Beta_value(Vector p, Vector pi, DenseMatrix sigma, Vector normal) {
  double Bvalue = B_value(p, pi, sigma);

  Vector rc(3);
  rc = 0.0;
  subtract(p, pi, rc); 

  return (normal * rc) / Bvalue;
}

double U_i_s(mfem::Vector p, mfem::Vector pi, mfem::DenseMatrix sigma0) {

  if (p(0) == pi(0) && p(1) == pi(1) && p(2) == pi(2)) {
    std::cout << "U_i_s(): p!=pi failed!";
    std::abort();
  }
  DenseMatrix rho0 = sigma0;
  rho0.Invert();
  double B0value = B_value(p, pi, sigma0);
  double Us;
  Us = (1.0 / (2 * EM::PI)) * std::pow(rho0.Det(), 0.5) *
       std::pow(B0value, -0.5);
  return Us;
}

void gradient_U_i_s(mfem::Vector p, mfem::Vector pi, mfem::DenseMatrix sigma0,
                    mfem::Vector &deltaUs) {

  DenseMatrix rho0 = sigma0;
  rho0.Invert();
  double B0value = B_value(p, pi, sigma0);

  double coef = (-1.0 / (2 * EM::PI)) * std::pow(rho0.Det(), 0.5) *
                std::pow(B0value, -1.5);
  Vector temp(3);
  temp = 0.0; 
  Vector r(3);
  r = 0.0;
  subtract(p, pi, r); // r = p - pi
  rho0.AddMult(r, temp);
  temp *= coef;
  deltaUs = temp;
}

void compute_normal(mfem::ElementTransformation &T,
                    const mfem::IntegrationPoint &ip, mfem::Vector &normal) {
  T.SetIntPoint(&ip);
  CalcOrtho(T.Jacobian(), normal);
  normal /= normal.Norml2(); // normal.unit();
}

void get_tetId_bdr(mfem::ParMesh *pmesh, int ElementNo, int &tet_id) {
  int flag = -1; 
  mfem::FaceElementTransformations *face_trans =
      pmesh->GetBdrFaceTransformations(ElementNo);
  mfem::ElementTransformation *face = face_trans->Face;
  mfem::ElementTransformation *e1 = face_trans->Elem1;
  if (e1 != NULL) {
    tet_id = e1->Attribute - 1;
    flag *= -1;
  }
  mfem::ElementTransformation *e2 = face_trans->Elem2;
  if (e2 != NULL) {
    tet_id = e2->Attribute - 1;
    flag *= -1;
  }
  myassert(flag == 1 && tet_id != -1);
}

// Function to get the current process's memory usage (in bytes)
double GetCurrentMemoryUsage() {
  rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  return double(usage.ru_maxrss / 1024.0 /
                1024.0); // ru_maxrss is in KB, converting to MB
}

bool equalDM(DenseMatrix a, DenseMatrix b) {
  myassert(a.Size() == b.Size());
  for (int i = 0; i < a.Size(); i++) {
    for (int j = 0; j < a.Size(); j++) {
      if (a(i, j) != b(i, j))
        return false;
    }
  }
  return true;
}

bool equalVector(Vector a, Vector b) {
  myassert(a.Size() == b.Size());
  for (int i = 0; i < a.Size(); i++) {
    if (a(i) != b(i))
      return false;
  }
  return true;
}

void Mult_DenseMatrix3(DenseMatrix &a, DenseMatrix &b, DenseMatrix &c) {
  myassert(a.Size() == b.Size());
  myassert(a.Size() == c.Size());
  myassert(a.Size() == 3);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        c(i, j) += a(i, k) * b(k, j);
      }
    }
  }
}