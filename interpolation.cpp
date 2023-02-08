#include <math.h>
#include <vector>
#include <iostream>
#include <set>
#include <fstream>
#include <string>
typedef struct {
  double x, y, z;
} DATA;


template<typename T>
struct Point {
  T x;
  T y;
  Point(T x, T y) : x(x), y(y) {}
  Point() = default;

  template<typename U>
  Point(const Point<U>& other) {
    x = static_cast<T>(other.x);
    y = static_cast<T>(other.y);
  }

  Point operator+=(const Point& other) {
    x += other.x;
    y += other.y;
    return *this;
  }

  Point operator-=(const Point& other) {
    x -= other.x;
    y -= other.y;
    return *this;
  }

  template<typename U>
  Point<U> operator-(const Point<U>& other) const {
    Point<U> copy(*this);
    copy -= other;
    return copy;
  }

  bool operator<(const Point& other) const {
    if (x < other.x) {
      return true;
    }
    return static_cast<bool>(x == other.x && y < other.y);
  }

  bool operator==(const Point& other) const {
    return x == other.x && y == other.y;
  }
};

template<typename T>
T ScalarMult(const Point<T>& v1, const Point<T>& v2) {
  return v1.x * v2.x + v1.y * v2.y;
}

template<typename T>
class Segment;

template<typename T>
class Line;

template<typename T>
class Line {
 public:
  Line(const Segment<T>& segment) {
    Point<T> vector = segment.point_b_ - segment.point_a_;
    c_ = -vector.y * segment.point_a_.x + vector.x * segment.point_a_.y;
    a_ = vector.y;
    b_ = -vector.x;
  }

  bool Parallel(const Line& other) { return a_ * other.b_ == other.a_ * b_; }

  bool Vertical() { return a_ == 0; }

  bool Same(const Line& other) {
    if (Parallel(other)) {
      return b_ * other.c_ == c_ * other.b_;
    }
    return false;
  }

  bool Contains(const Point<T>& point) const {
    return a_ * point.x + b_ * point.y + c_ == 0;
  }

  template<typename U>
  Point<U> PointOfIntersection(const Line<T>& second) {
    U denominator = static_cast<U>(a_ * second.b_ - b_ * second.a_);
    U x = (static_cast<U>(-c_ * second.b_ + second.c_ * b_)) / denominator;
    U y = (static_cast<U>(c_ * second.a_ - second.c_ * a_)) / denominator;
    return Point<U>(x, y);
  }

 private:
  T a_;
  T b_;
  T c_;
  friend Segment<T>;
};

template<typename T>
class Segment {
 public:
  Segment(const Point<T>& point_a, const Point<T>& point_b)
      : point_a_(point_a), point_b_(point_b) {
    ChangePlace();
  };

  bool SameBound(const Segment<T>& other) const {
    return static_cast<bool>((point_a_ == other.point_a_) ||
        (point_b_ == other.point_b_));
  }

  void ChangePlace() {
    if (point_b_ < point_a_) {
      std::swap(point_b_, point_a_);
    }
  }

  bool operator==(const Segment<T>& other) const {
    return static_cast<bool>(point_a_ == other.point_a_ &&
        point_b_ == other.point_b_);
  }

  bool IsPoint() const { return point_a_ == point_b_; }

  bool Left(const Segment& other) const { return point_b_ < other.point_a_; }
  bool Lower(const Segment& other) const { return point_b_ < other.point_a_; }

  bool IntersectionInside(const Point<long double>& intersection,
                          const Segment<T>& second) {
    Point<long double> vector1 = point_a_ - intersection;
    Point<long double> vector2 = point_b_ - intersection;
    Point<long double> vector3 = second.point_a_ - intersection;
    Point<long double> vector4 = second.point_b_ - intersection;
    return ScalarMult(vector1, vector2) <= 0 &&
        ScalarMult(vector3, vector4) <= 0;
  }

  bool Intersect(const Segment<T>& second,
                 bool boundaries) { //boundaries = true если при нахождении скважины на разломе НЕ нужно ее учитывать
    if (SameBound(second) && boundaries) {
      return true;
    }
    if (SameBound(second)) {
      return false;
    }
    Line<T> first_line(*this);
    Line<T> second_line(second);
    if (!boundaries) {
      if (first_line.Contains(second.point_a_)
          || first_line.Contains(second.point_b_)
          || second_line.Contains(point_a_) || second_line.Contains(point_b_)) {
        return false;
      }
    }
    if (first_line.Same(second_line)) {
      if (Left(second) || second.Left(*this)) {
        return false;
      }
      if (first_line.Vertical() && (Lower(second) || second.Lower(*this))) {
        return false;
      }
      if (IsPoint()) {
        return second_line.Contains(point_a_);
      }
      if (second.IsPoint()) {
        return first_line.Contains(second.point_a_);
      }
      return true;
    }
    if (first_line.Parallel(second_line)) {
      return false;
    }
    Point<long double> intersection =
        first_line.template PointOfIntersection<long double>(second_line);
    return IntersectionInside(intersection, second);
  }

 private:
  friend Line<T>;
  Point<T> point_a_;
  Point<T> point_b_;
};

struct Polygon {
  std::vector<Point<double>> points;
  Polygon(const std::vector<Point<double>>& points) : points(points) {};
};

double INF = 1000000;
class Graph {
 public:
  std::vector<double> dist;
  int n;
  std::vector<std::vector<std::pair<int, double>>> g;
  Graph(int n, std::vector<std::vector<std::pair<int, double>>>& g)
      : n(n), g(g) {
    dist.assign(n, INF);
  }

  void Djikstra(int start) {
    dist[start] = 0;
    std::set<std::pair<double, int>> near;
    near.insert(std::make_pair(dist[start], start));
    while (!near.empty()) {
      int vert = near.begin()->second;
      near.erase(near.begin());
      for (int i = 0; i < g[vert].size(); ++i) {
        int next = g[vert][i].first;
        double length = g[vert][i].second;
        if (dist[vert] + length < dist[next]) {
          near.erase(std::make_pair(dist[next], next));
          dist[next] = dist[vert] + length;
          near.insert(std::make_pair(dist[next], next));
        }
      }
    }
  }

  void ClearDist() {
    dist.assign(n, INF);
  }
};

struct Interpolation {
  DATA* aData;
  int anData;
  int nx;
  int ny;
  int dx;
  int dy;
  int xleft;
  int ytop;
  double apower;
  double asmoosnes;
  double* aout;
  std::vector<Polygon> polygons;
  int count_of_polygons;
  std::vector<std::vector<double>> distances;
  std::vector<Point<double>> vertice;

  Interpolation(DATA* aData,
                int anData,
                int nx,
                int ny,
                int dx,
                int dy,
                int xleft,
                int ytop,
                double apower,
                double asmoosmes,
                double* aout,
                const std::vector<Polygon>& polygons)
      : aData(aData),
        anData(anData),
        nx(nx),
        ny(ny),
        dx(dx),
        dy(dy),
        xleft(xleft),
        ytop(ytop),
        apower(apower),
        asmoosnes(asmoosmes),
        aout(aout),
        polygons(polygons) {
    distances =
        std::vector<std::vector<double>>(nx * ny, std::vector<double>(anData));
    count_of_polygons = polygons.size();
    vertice = std::vector<Point<double>>(nx * ny + anData);
    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        vertice[i * ny + j].x = xleft + dx * i;
        vertice[i * ny + j].y = ytop + dy * j;
      }
    }
    for (int i = 0; i < anData; ++i) {
      vertice[nx * ny + i].x = aData[i].x;
      vertice[nx * ny + i].y = aData[i].y;
    }
  }

  void InverseDistance() { //сама функция интерполяции
    int i, j, k;
    double h, ph, xn, yn, sum, sumz;
    xn = xleft;
    for (i = 0; i < nx; i++) {
      yn = ytop;
      for (j = 0; j < ny; j++) {
        sum = 0;
        sumz = 0;
        for (k = 0; k < anData; k++) {
          h = sqrt(distances[i * ny + j][k] + asmoosnes * asmoosnes);
          if (h == fabs(asmoosnes)) {
            sumz = aData[k].z;
            sum = 1;
            break;
          }
          ph = 1 / pow(h, apower);
          sum += ph;
          sumz += aData[k].z * ph;
        }
        if (fabs(sum) > 0)
          aout[i * ny + j] = sumz / sum;
        yn += dy;
      }
      xn += dx;
    }
  }

  void FindDistances(bool boundaries) {
    //первый этап включает в себя построение графа по сетке и пробуренным вершинам, убирая те ребра, которые проходят через разломы
    // из каждого узла сетки идет 8 ребер в ближайшие узлы (вверх, вверх-влево, влево, влево-вниз и т.д.
    //  из каждой скважины идет 4 ребра в узлы клетки, в которой она находится
    std::vector<std::vector<std::pair<int, double>>>
        graph(nx * ny + anData, std::vector<std::pair<int, double>>());
    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        std::vector<bool> possible_connections(8, true);
        for (int k = 0; k < count_of_polygons; ++k) {
          for (int l = 1; l < polygons[k].points.size(); ++l) {
            Point cut_begin
                (polygons[k].points[l - 1].x, polygons[k].points[l - 1].y);
            Point cut_end(polygons[k].points[l].x, polygons[k].points[l].y);
            Segment curr_cut(cut_begin, cut_end);
            if (j > 0) {
              if (curr_cut.Intersect(Segment<double>(vertice[i * ny + j],
                                                     vertice[i * ny + j - 1]),
                                     boundaries)) {
                possible_connections[0] = false;
              }
            }
            if (i > 0) {
              if (curr_cut.Intersect(Segment<double>(vertice[i * ny + j],
                                                     vertice[(i - 1) * ny + j]),
                                     boundaries)) {
                possible_connections[1] = false;
              }
            }
            if (j < ny - 1) {
              if (curr_cut.Intersect(Segment<double>(vertice[i * ny + j],
                                                     vertice[i * ny + j + 1]),
                                     boundaries)) {
                possible_connections[2] = false;
              }
            }
            if (i < nx - 1) {
              if (curr_cut.Intersect(Segment<double>(vertice[i * ny + j],
                                                     vertice[(i + 1) * ny + j]),
                                     boundaries)) {
                possible_connections[3] = false;
              }
            }
            if (j > 0 && i > 0) {
              if (curr_cut.Intersect(Segment<double>(vertice[i * ny + j],
                                                     vertice[(i - 1) * ny + j
                                                         - 1]), boundaries)) {
                possible_connections[4] = false;
              }
            }
            if (j > 0 && i < nx - 1) {
              if (curr_cut.Intersect(Segment<double>(vertice[i * ny + j],
                                                     vertice[(i + 1) * ny + j
                                                         - 1]), boundaries)) {
                possible_connections[5] = false;
              }
            }
            if (j < ny - 1 && i > 0) {
              if (curr_cut.Intersect(Segment<double>(vertice[i * ny + j],
                                                     vertice[(i - 1) * ny + j
                                                         + 1]), boundaries)) {
                possible_connections[6] = false;
              }
            }
            if (j < ny - 1 && i < nx - 1) {
              if (curr_cut.Intersect(Segment<double>(vertice[i * ny + j],
                                                     vertice[(i + 1) * ny + j
                                                         + 1]), boundaries)) {
                possible_connections[7] = false;
              }
            }
          }
        }
        if (j > 0 && possible_connections[0]) {
          graph[i * ny + j].push_back(std::make_pair(i * ny + j - 1, dy));
        }
        if (i > 0 && possible_connections[1]) {
          graph[i * ny + j].push_back(std::make_pair((i - 1) * ny + j, dx));
        }
        if (j < ny - 1 && possible_connections[2]) {
          graph[i * ny + j].push_back(std::make_pair(i * ny + j + 1, dy));
        }
        if (i < nx - 1 && possible_connections[3]) {
          graph[i * ny + j].push_back(std::make_pair((i + 1) * ny + j, dx));
        }
        if (j > 0 && i > 0 && possible_connections[4]) {
          graph[i * ny + j].push_back(std::make_pair((i - 1) * ny + j - 1,
                                                     pow(dy * dy + dx * dx,
                                                         0.5)));
        }
        if (j > 0 && i < nx - 1 && possible_connections[5]) {
          graph[i * ny + j].push_back(std::make_pair((i + 1) * ny + j - 1,
                                                     pow(dy * dy + dx * dx,
                                                         0.5)));
        }
        if (j < ny - 1 && i > 0 && possible_connections[6]) {
          graph[i * ny + j].push_back(std::make_pair((i - 1) * ny + j + 1,
                                                     pow(dy * dy + dx * dx,
                                                         0.5)));
        }
        if (j < ny - 1 && i < nx - 1 && possible_connections[7]) {
          graph[i * ny + j].push_back(std::make_pair((i + 1) * ny + j + 1,
                                                     pow(dy * dy + dx * dx,
                                                         0.5)));
        }
      }
    }
    for (int i = 0; i < anData; ++i) {
      int left_x_ind = ((int) aData[i].x - xleft) / dx;
      int lowest_y_ind = ((int) aData[i].y - ytop) / dy;
      std::vector<bool> possible_connections = {true, true, true, true};
      double dist1 = pow(pow(aData[i].x - xleft - dx * left_x_ind, 2)
                             + pow(aData[i].y - ytop - dy * lowest_y_ind, 2),
                         0.5);
      double dist2 = pow(pow(dx - (aData[i].x - xleft - dx * left_x_ind), 2)
                             + pow((aData[i].y - ytop - dy * lowest_y_ind), 2),
                         0.5);
      double dist3 = pow(pow((aData[i].x - xleft - dx * left_x_ind), 2)
                             + pow((dy
                                 - (aData[i].y - ytop - dy * lowest_y_ind)), 2),
                         0.5);
      double dist4 = pow(pow((dx - (aData[i].x - xleft - dx * left_x_ind)), 2)
                             + pow(dy - (aData[i].y - ytop - dy * lowest_y_ind),
                                   2), 0.5);
      for (int k = 0; k < count_of_polygons; ++k) {
        for (int l = 1; l < polygons[k].points.size(); ++l) {
          Point
              cut_begin
              (polygons[k].points[l - 1].x, polygons[k].points[l - 1].y);
          Point cut_end(polygons[k].points[l].x, polygons[k].points[l].y);
          Segment curr_cut(cut_begin, cut_end);

          if (curr_cut.Intersect(Segment<double>(Point(
                                                     (double) left_x_ind * dx + xleft,
                                                     (double) lowest_y_ind * dy + ytop),
                                                 Point(aData[i].x,
                                                       aData[i].y)),
                                 boundaries)) {
            possible_connections[0] = false;
          }
          if (left_x_ind + 1 < nx && curr_cut.Intersect(Segment<double>(Point(
                                                                            (double) (left_x_ind + 1) * dx + xleft,
                                                                            (double) lowest_y_ind * dy + ytop),
                                                                        Point(
                                                                            aData[i].x,
                                                                            aData[i].y)),
                                                        boundaries)) {
            possible_connections[1] = false;
          }
          if (lowest_y_ind + 1 < ny && curr_cut.Intersect(Segment<double>(Point(
                                                                              (double) left_x_ind * dx + xleft,
                                                                              (double) (lowest_y_ind + 1) * dy + ytop),
                                                                          Point(
                                                                              aData[i].x,
                                                                              aData[i].y)),
                                                          boundaries)) {
            possible_connections[2] = false;
          }
          if (left_x_ind + 1 < nx && lowest_y_ind + 1 < ny
              && curr_cut.Intersect(Segment<double>(Point(
                                                        (double) (left_x_ind + 1) * dx + xleft,
                                                        (double) (lowest_y_ind + 1) * dy + ytop),
                                                    Point(aData[i].x,
                                                          aData[i].y)),
                                    boundaries)) {
            possible_connections[3] = false;
          }
        }
      }
      if (possible_connections[0]) {
        graph[nx * ny + i].push_back(std::make_pair(
            left_x_ind * ny + lowest_y_ind, dist1));
        graph[left_x_ind * ny + lowest_y_ind].push_back(std::make_pair(
            nx * ny + i, dist1));
      }
      if (possible_connections[1] && left_x_ind + 1 < nx) {
        graph[nx * ny + i].push_back(std::make_pair(
            (left_x_ind + 1) * ny + lowest_y_ind, dist2));
        graph[(left_x_ind + 1) * ny + lowest_y_ind].push_back(std::make_pair(
            nx * ny + i,
            dist2));
      }
      if (possible_connections[2] && lowest_y_ind + 1 < ny) {
        graph[nx * ny + i].push_back(std::make_pair(
            left_x_ind * ny + lowest_y_ind + 1, dist3));
        graph[left_x_ind * ny + lowest_y_ind + 1].push_back(std::make_pair(
            nx * ny + i, dist3));
      }
      if (possible_connections[3] && left_x_ind + 1 < nx
          && lowest_y_ind + 1 < ny) {
        graph[nx * ny + i].push_back(std::make_pair(
            (left_x_ind + 1) * ny + lowest_y_ind + 1, dist4));
        graph[(left_x_ind + 1) * ny + lowest_y_ind
            + 1].push_back(std::make_pair(nx * ny + i, dist4));
      }
    }
    //построение закончилось
    //и теперь запускаем дейкстру(поиск кратчайших путей) по существующим ребрам из всех вершин aData
    Graph new_graph(nx * ny + anData, graph);
    for (int i = nx * ny; i < nx * ny + anData; ++i) {
      new_graph.Djikstra(i);
      for (int j = 0; j < nx; ++j) {
        for (int k = 0; k < ny; ++k) {
          distances[j * ny + k][i - nx * ny] =
              pow(new_graph.dist[j * ny + k], 2);
          bool straight = true;
          for (int p = 0; p < count_of_polygons; ++p) {
            for (int ind = 1; ind < polygons[p].points.size(); ++ind) {
              Point cut_begin
                  (polygons[p].points[ind - 1].x,
                   polygons[p].points[ind - 1].y);
              Point
                  cut_end(polygons[p].points[ind].x, polygons[p].points[ind].y);
              Segment curr_cut(cut_begin, cut_end);
              if (curr_cut.Intersect(Segment<double>(vertice[j * ny + k],
                                                     vertice[i]), boundaries)) {
                straight = false;
              }
            }
          }
          if (straight) {
            distances[j * ny + k][i - nx * ny] =
                pow((vertice[j * ny + k] - vertice[i]).x, 2)
                    + pow((vertice[j * ny + k] - vertice[i]).y, 2);
          }
        }
      }
      new_graph.ClearDist();
    }
  }

};

int main(int argc, char* argv[]) {
  std::string params = argv[1];
  std::string data;
  std::string faults;
  std::string line;
  std::string heights;
  std::ifstream in(params);
  int anData;
  int nx;
  int ny;
  int dx;
  int dy;
  int xleft;
  int ytop;
  double apower;
  double asmoosnes;
  int boundaries;
  if (in.is_open()) {
    in >> anData >> nx >> ny >> dx >> dy >> xleft >> ytop >> apower >> asmoosnes
       >> boundaries >> faults >> data >> heights;
  }
  in.close();
  std::ifstream in1(data);
  std::string well, name;
  DATA* aData = new DATA[anData];
  double x, y, z;
  if (in1.is_open()) {
    in1 >> well >> name >> name >> name;
    for (int i = 0; i < anData; ++i) {
      getline(in1, line);
      in1 >> well >> x >> y >> z;
      aData[i].x = x;
      aData[i].y = y;
      aData[i].z = z;
    }
  }
  in1.close();
  std::fstream in2(faults);
  std::vector<Polygon> polygons;
  std::vector<Point<double>> curr_polygon;
  if (in2.is_open()) {
    while (in2 >> x >> y >> z) {
      if (x != 999.0) {
        curr_polygon.push_back(Point<double>(x, y));
      } else {
        polygons.push_back(curr_polygon);
        curr_polygon.clear();
      }
    }
  }
  in2.close();
  double* aout = new double[nx * ny];
  Interpolation interpolation(aData,
                              anData,
                              nx,
                              ny,
                              dx,
                              dy,
                              xleft,
                              ytop,
                              apower,
                              asmoosnes,
                              aout,
                              polygons);
  interpolation.FindDistances(boundaries);
  interpolation.InverseDistance();
  std::ofstream out;
  out.open(heights);
  if (out.is_open()) {
    out << -996 << " " << ny << " " << dx << " " << dy << "\n";
    out << xleft << " " << xleft + dx * (nx - 1) << " " << ytop << " "
        << ytop + dy * (ny - 1) << "\n";
    out << nx << " " << 0 << " " << xleft << " " << ytop << "\n";
    out << "0   0   0   0   0   0   0\n";
    for (int i = 0; i < ny; ++i) {
      for (int j = 0; j < nx; ++j) {
//        out<<xleft+j*dx<<" "<<ytop+i*dy<<" "<<aout[j*ny+i]<<"\n";
        out << aout[j * ny + i] << " ";
      }
      out << "\n";
    }
  }
}
