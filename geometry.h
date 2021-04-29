//
// Created by Savichev Dmitrii on 30.11.2020.
//

#ifndef GEOMETRY_GEOMETRY_H
#define GEOMETRY_GEOMETRY_H
#include <iostream>
#include <cmath>
#include <vector>
const double precision = 1e-6;

class Line;

struct Point {
    double x = 0;
    double y = 0;
    Point() = default;
    Point(const double& x, const double& y): x(x), y(y) {};

    Point rotate(const Point& center, double angle) const;
    void scale(const Point& center, double coefficient);
    void reflex(const Line axis);
};

class Line {
public:
    double a;
    double b;
    double c;

    Line(const double k, const double b): a(k), b(-1), c(b) {};
    Line(const double a, const double b, const double c): a(a), b(b), c(c) {};
    Line(const Point& first, const Point& second): a(second.y - first.y), b(first.x - second.x), c(second.x * first.y - second.y * first.x) {};
    Line(const Point& point, const double k): a(k), b(-1), c(point.y - k * point.x) {};
};

class Shape {
public:
    virtual double perimeter() const = 0;
    virtual double area() const = 0;
    virtual bool operator ==(const Shape& another) const = 0;
    bool operator != (const Shape& another) const {
        return !(*this == another);
    }
    virtual bool isCongruentTo(const Shape& another) const = 0;
    virtual bool isSimilarTo(const Shape& another) const = 0;
    virtual bool containsPoint(Point point) const = 0;
    virtual void rotate(Point center, double angle) = 0;
    virtual void reflex(Point center) = 0;
    virtual void reflex(Line axis) = 0;
    virtual void scale(Point center, double coefficient) = 0;
    virtual ~Shape() = default;
};

class Polygon: public Shape {
public:
    Polygon(const std::vector <Point>& v): vertices(v) {};
    Polygon (std::initializer_list <Point> v): vertices(v) {};
    int verticesCount() const;

    std::vector <Point> getVertices() const;

    bool isConvex() const;
    double perimeter() const override;
    double area() const override;
    bool operator ==(const Shape& another) const override;
    bool isCongruentTo(const Shape& another) const override;
    bool isSimilarTo(const Shape& another) const override;
    bool containsPoint(Point point) const override;
    void rotate(Point center, double angle) override;
    void reflex(Point center) override;
    void reflex(Line axis) override;
    void scale(Point center, double coefficient) override;
protected:
    std::vector <Point> vertices;
private:
    void move(double dx, double dy);
};

class Ellipse: public Shape {
public:
    Ellipse(const Point& first, const Point& second, double x);

    std::pair<Point, Point> focuses() const;
    std::pair<Line, Line> directrices() const;
    double eccentricity() const;
    Point center() const;

    double perimeter() const override;
    double area() const override;
    bool operator ==(const Shape& another) const override;
    bool isSimilarTo(const Shape& another) const override;
    bool isCongruentTo(const Shape& another) const override;
    bool containsPoint(Point point) const override;
    void rotate(Point center, double angle) override;
    void reflex(Point center) override;
    void reflex(Line axis) override;
    void scale(Point center, double coefficient) override;
protected:
    std::pair<Point, Point> foci;
    double a;
    double b;
};

class Circle: public Ellipse {
public:
    Circle(const Point& center, const double radius): Ellipse(center, center, 2 * radius) {};
    double radius() const {
        return a;
    }
};

class Rectangle: public Polygon {
public:
    Rectangle(const Point& first, const Point& second, double ratio);

    Point center() const;
    std::pair <Line, Line> diagonals() const;
};

class Square: public Rectangle {
public:
    Square(const Point& first, const Point& second);

    Circle circumscribedCircle() const;
    Circle inscribedCircle() const;
};

class Triangle: public Polygon {
public:
    Triangle(const Point& first, const Point& second, const Point& third): Polygon{first, second, third} {}

    Circle circumscribedCircle() const;
    Circle inscribedCircle() const;

    Point centroid() const;
    Point orthocenter() const;
    Line EulerLine() const;
    Circle ninePointsCircle() const;
};

//=========Arithmetics===========
bool equal(double first, double second) {
    return std::abs(first - second) <= precision;
}

bool between(double x, double first, double second) {
    return equal(std::abs(first - x) + std::abs(second - x), std::abs(first - second));
}

double sqr(double a) {
    return a * a;
}
//========Points&Lines===========
double dist(const Point& first, const Point& second) {
    return sqrt((first.x - second.x) * (first.x - second.x) + (first.y - second.y) * (first.y - second.y));
}

Point Point::rotate(const Point& center, double angle) const {
    double ans_x = center.x + cos(angle) * (x - center.x) - sin(angle) * (y - center.y);
    double ans_y = center.y + sin(angle) * (x - center.x) + cos(angle) * (y - center.y);
    return Point(ans_x, ans_y);
}

void Point::scale(const Point& center, double coefficient) {
    double ans_x = (x - center.x) * coefficient + center.x;
    double ans_y = (y - center.y) * coefficient + center.y;
    *this = Point(ans_x, ans_y);
}

bool operator == (const Point& first, const Point& second) {
    return equal(first.x, second.x) && equal(first.y, second.y);
}

bool operator != (const Point& first, const Point& second) {
    return !(first == second);
}

bool operator == (const Line& first, const Line& second) {
    double coef;
    if (second.a != 0) {
        coef = first.a / second.a;
    } else {
        coef = first.b / second.b;
    }
    return equal(first.a, second.a * coef) && equal(first.b, second.b * coef) && equal(first.c, second.c * coef);
}

bool operator != (const Line& first, const Line& second) {
    return !(first == second);
}

Point intersection(const Line& first, const Line& second) {
    Point tmp;
    tmp.x = -(first.c * second.b - second.c * first.b) / (first.a * second.b - second.a * first.b);
    tmp.y = -(first.a * second.c - second.a * first.c) / (first.a * second.b - second.a * first.b);
    return tmp;
}

Line perpendicular(const Line& line) {
    return Line(line.b, -line.a, 0);
}

Line middlePerpendecular(const Point& first, const Point& second) {
    Point middle((first.x + second.x) / 2, (first.y + second.y) / 2);
    Line axis(first, second);
    Line result = perpendicular(axis);
    result.c = -result.a * middle.x - result.b * middle.y;
    return result;
}

Point middlePoint(const Point& first, const Point& second) {
    return Point((first.x + second.x) / 2, (first.y + second.y) / 2);
}

double angle(const Point first, const Point second) {
    double cos = (first.x * second.x + first.y * second.y) / (sqrt(sqr(first.x) + sqr(first.y)) * sqrt(sqr(second.x) + sqr(second.y)));
    double sin = (first.x * second.y - first.y * second.x) / (sqrt(sqr(first.x) + sqr(first.y)) * sqrt(sqr(second.x) + sqr(second.y)));
    if (sin < 0) {
        return -acos(cos) * 180 / M_PI;
    } else {
        return acos(cos) * 180 / M_PI;
    }
}

void Point::reflex(const Line axis) {
    Line tmp = perpendicular(axis);
    tmp.c = -tmp.a * x - tmp.b * y;
    scale(intersection(tmp, axis), -1);
}

//===========Polygon=============

int Polygon::verticesCount() const {
    return vertices.size();
}

std::vector <Point> Polygon::getVertices() const {
    return vertices;
}

bool Polygon::isConvex() const {
    double det = vertices[0].x * vertices[1].y + vertices[0].y * vertices[2].x + vertices[1].x * vertices[2].y - vertices[0].x * vertices[2].y - vertices[0].y * vertices[1].x - vertices[1].y * vertices[2].x;
    int sign;
    if (det > 0) {
        sign = 1;
    } else {
        sign = -1;
    }
    int n = verticesCount();
    for (int i = 1; i < n; ++i) {
        det = vertices[i].x * vertices[(i + 1) % n].y + vertices[i].y * vertices[(i + 2) % n].x + vertices[(i + 1) % n].x * vertices[(i + 2) % n].y - vertices[i].x * vertices[(i + 2) % n].y - vertices[i].y * vertices[(i + 1) % n].x - vertices[(i + 1) % n].y * vertices[(i + 2) % n].x;
        if ((det > 0 && sign == -1) || (det < 0 && sign == 1)) {
            return false;
        }
    }
    return true;
}

double Polygon::perimeter() const {
    double p = 0;
    for (int i = 0; i < verticesCount(); ++i) {
        p += dist(vertices[i], vertices[(i + 1) % verticesCount()]);
    }
    return p;
}

double Polygon::area() const {
    double s = 0;
    for (int i = 0; i < verticesCount(); ++i) {
        s += vertices[i].x * vertices[(i + 1) % verticesCount()].y - vertices[i].y * vertices[(i + 1) % verticesCount()].x;
    }
    return std::abs(s / 2);
}

bool Polygon::operator ==(const Shape& another) const {
    const Polygon* polygon_ptr = dynamic_cast<const Polygon*>(&another);
    if (polygon_ptr == nullptr || verticesCount() != polygon_ptr->verticesCount()) {
        return false;
    }
    Point start_point = polygon_ptr->vertices[0];
    int point_index = -1;
    for (int i = 0; i < verticesCount(); ++i) {
        if (vertices[i] == start_point) {
            point_index = i;
            break;
        }
    }
    if (point_index == -1) {
        return false;
    }
    bool ok = true;
    for(int i = 0; i < verticesCount(); ++i) {
        if (polygon_ptr->vertices[i] != vertices[(point_index + i) % verticesCount()]) {
            ok = false;
            break;
        }
    }
    if (ok) {
        return true;
    }
    ok = true;
    for(int i = 0; i < verticesCount(); ++i) {
        if (polygon_ptr->vertices[i] != vertices[(verticesCount() + point_index - i) % verticesCount()]) {
            ok = false;
            break;
        }
    }
    return ok;
}

bool Polygon::isCongruentTo(const Shape& another) const {
    const Polygon* polygon_ptr = dynamic_cast<const Polygon*>(&another);
    if (polygon_ptr == nullptr || polygon_ptr->verticesCount() != verticesCount()) {
        return false;
    }

    for (int i = 0; i < verticesCount(); ++i) {
        if (!equal(dist(vertices[i], vertices[(i + 1) % verticesCount()]), dist(polygon_ptr->vertices[0], polygon_ptr->vertices[1]))) {
            continue;
        }
        Polygon tmp(polygon_ptr->vertices);
        tmp.move(vertices[i].x - tmp.vertices[0].x, vertices[i].y - tmp.vertices[0].y);
        Line axis(tmp.vertices[0], vertices[(i + 1) % verticesCount()]);
        Point tmp_vec(tmp.vertices[1].x - tmp.vertices[0].x, tmp.vertices[1].y - tmp.vertices[0].y);
        Point vec(vertices[(i + 1) % verticesCount()].x - tmp.vertices[0].x, vertices[(i + 1) % verticesCount()].y - tmp.vertices[0].y);
        tmp.rotate(tmp.vertices[0], angle(tmp_vec, vec));
        if (*this == tmp) {
            return true;
        }
        tmp.reflex(middlePerpendecular(tmp.vertices[0], tmp.vertices[1]));
        if (*this == tmp) {
            return true;
        }
        tmp.reflex(axis);
        if (*this == tmp) {
            return true;
        }
        tmp.reflex(middlePerpendecular(tmp.vertices[0], tmp.vertices[1]));
        if (*this == tmp) {
            return true;
        }
    }
    return false;
}

bool Polygon::isSimilarTo(const Shape& another) const {
    const Polygon* polygon_ptr = dynamic_cast<const Polygon*>(&another);
    if (polygon_ptr == nullptr || polygon_ptr->verticesCount() != verticesCount()) {
        return false;
    }
    double coefficient = perimeter() / polygon_ptr->perimeter();
    Polygon tmp(polygon_ptr->vertices);
    tmp.scale(tmp.vertices[0], coefficient);
    if (isCongruentTo(tmp)) {
        return true;
    }

    return false;
}

bool Polygon::containsPoint(Point point) const {
    Line axis(point, 0);
    int count = 0;
    for (int i = 0; i < verticesCount(); ++i) {
        if (point == vertices[i]) {
            return true;
        }
        Point min_point = vertices[i].y < vertices[(i + 1) % verticesCount()].y ? vertices[i] : vertices[(i + 1) % verticesCount()];
        Line tmp(vertices[i], vertices[(i + 1) % verticesCount()]);
        if (tmp != axis) {
            Point tmp_point = intersection(axis, tmp);
            if (tmp_point != min_point && tmp_point.x > point.x && between(tmp_point.x, vertices[i].x, vertices[(i + 1) % verticesCount()].x) && between(tmp_point.y, vertices[i].y, vertices[(i + 1) % verticesCount()].y)) {
                ++count;
            }
        }
    }
    return count % 2;
}

void Polygon::rotate(Point center, double angle) {
    angle = angle / 180 * M_PI;
    for (int i = 0; i < verticesCount(); ++i) {
        vertices[i] = vertices[i].rotate(center, angle);
    }
};
void Polygon::reflex(Point center) {
    scale(center, -1);
}

void Polygon::reflex(Line axis) {
    for (int i = 0; i < verticesCount(); ++i) {
    vertices[i].reflex(axis);
    }
}

void Polygon::scale(Point center, double coefficient) {
    for (int i = 0; i < verticesCount(); ++i) {
        vertices[i].scale(center, coefficient);
    }
}

void Polygon::move(double dx, double dy) {
    for (int i = 0; i < verticesCount(); ++i) {
        vertices[i].x += dx;
        vertices[i].y += dy;
    }
}

//==========Ellipse==============

Ellipse::Ellipse(const Point& first, const Point& second, double x): foci({first, second}), a(x / 2), b(sqrt(sqr(a) - sqr(dist(first, second) / 2))) {};

std::pair<Point, Point> Ellipse::focuses() const {
    return foci;
}

std::pair<Line, Line> Ellipse::directrices() const {
    Line axis(foci.first, foci.second);
    std::pair <Line, Line> result {perpendicular(axis), perpendicular(axis)};
    result.first.c = a * a / sqrt(a * a - b * b) * sqrt(result.first.a * result.first.a + result.first.b * result.first.b) - result.first.a * center().x - result.first.b * center().y;
    result.second.c = - a * a / sqrt(a * a - b * b) * sqrt(result.second.a * result.second.a + result.second.b * result.second.b) - result.second.a * center().x - result.second.b * center().y;
    return result;
}

double Ellipse::eccentricity() const {
    return sqrt(1 - sqr(b) / sqr(a));
}

Point Ellipse::center() const {
    return Point((foci.first.x + foci.second.x) / 2, (foci.first.y + foci.second.y) / 2);
}

double Ellipse::perimeter() const {
    return std::comp_ellint_2(eccentricity()) * 4 * a;
}

double Ellipse::area() const {
    return M_PI * a * b;
}

bool Ellipse::operator ==(const Shape& another) const {
    const Ellipse* ellipse_ptr = dynamic_cast<const Ellipse*>(&another);
    if (ellipse_ptr != nullptr) {
        return ((foci.first == ellipse_ptr->foci.first && foci.second == ellipse_ptr->foci.second) ||
                (foci.first == ellipse_ptr->foci.second && foci.second == ellipse_ptr->foci.first)) && equal(a, ellipse_ptr->a);
    }
    return false;
}

bool Ellipse::isSimilarTo(const Shape& another) const {
    const Ellipse* ellipse_ptr = dynamic_cast<const Ellipse*>(&another);
    if (ellipse_ptr != nullptr) {
        return equal(eccentricity(), ellipse_ptr->eccentricity());
    }
    return false;
}

bool Ellipse::isCongruentTo(const Shape& another) const {
    const Ellipse* ellipse_ptr = dynamic_cast<const Ellipse*>(&another);
    if (ellipse_ptr != nullptr) {
        return equal(a, ellipse_ptr->a) && equal(b, ellipse_ptr->b);
    }
    return false;
}

bool Ellipse::containsPoint(Point point) const {
    double foci_dist = dist(point, foci.first) + dist(point, foci.second);
    return foci_dist < 2 * a || equal(foci_dist, 2 * a);
}

void Ellipse::rotate(Point center, double angle) {
    angle = angle / 180 * M_PI;
    foci.first = foci.first.rotate(center, angle);
    foci.second = foci.second.rotate(center, angle);
}

void Ellipse::reflex(Point center) {
    scale(center, -1);
}

void Ellipse::reflex(Line axis) {
    foci.first.reflex(axis);
    foci.second.reflex(axis);
}

void Ellipse::scale(Point center, double coefficient) {
    foci.first.scale(center, coefficient);
    foci.second.scale(center, coefficient);
    a *= std::abs(coefficient);
    b *= std::abs(coefficient);
}

//===========Rectangle============
Rectangle::Rectangle(const Point& first, const Point& second, double ratio): Polygon{first, Point(0, 0), second, Point(0, 0)} {
    double angle = 180 - 2 * atan(ratio) * 180 / M_PI;
    Point middle = middlePoint(first, second);
    vertices[1] = first.rotate(middle, angle);
    vertices[3] = second.rotate(middle, angle);
}

Point Rectangle::center() const {
    return middlePoint(vertices[0], vertices[2]);
}

std::pair <Line, Line> Rectangle::diagonals() const {
    return {Line(vertices[0], vertices[2]), Line(vertices[1], vertices[3])};
}

//=======Square===========

Square::Square(const Point& first, const Point& second): Rectangle(first, second, 1){}

Circle Square::circumscribedCircle() const {
    return Circle(middlePoint(vertices[0], vertices[2]), dist(vertices[0], vertices[2]) / 2);
}

Circle Square::inscribedCircle() const {
    return Circle(middlePoint(vertices[0], vertices[2]), dist(vertices[0], vertices[1]) / 2);
}

//=======Triangle==========

Circle Triangle::circumscribedCircle() const {
    Point center = intersection(middlePerpendecular(vertices[0],vertices[1]), middlePerpendecular(vertices[0], vertices[2]));
    return Circle(center, dist(center, vertices[0]));
}

Circle Triangle::inscribedCircle() const {
    double a = dist(vertices[1], vertices[2]);
    double b = dist(vertices[2], vertices[0]);
    double c = dist(vertices[0], vertices[1]);
    Point tmp;
    tmp.x = (vertices[0].x * a + vertices[1].x * b + vertices[2].x * c) / (a + b + c);
    tmp.y = (vertices[0].y * a + vertices[1].y * b + vertices[2].y * c) / (a + b + c);
    double radius = 2 * area() / perimeter();
    return Circle(tmp, radius);
}

Point Triangle::centroid() const {
    return Point((vertices[0].x + vertices[1].x + vertices[2].x) / 3, (vertices[0].y + vertices[1].y + vertices[2].y) / 3);
}

Point Triangle::orthocenter() const {
    Line a(vertices[0], vertices[1]);
    Line ha = perpendicular(a);
    ha.c = -ha.a * vertices[2].x - ha.b * vertices[2].y;
    Line b(vertices[1], vertices[2]);
    Line hb = perpendicular(b);
    hb.c = -hb.a * vertices[0].x - hb.b * vertices[0].y;
    return intersection(ha, hb);
}

Line Triangle::EulerLine() const {
    Point center = circumscribedCircle().center();
    return Line(center, orthocenter());
}

Circle Triangle::ninePointsCircle() const {
    return Triangle{middlePoint(vertices[0], vertices[1]), middlePoint(vertices[1], vertices[2]), middlePoint(vertices[2], vertices[0])}.circumscribedCircle();
}

#endif //GEOMETRY_GEOMETRY_H
