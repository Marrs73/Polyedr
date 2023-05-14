from math import pi, acos
from functools import reduce
from operator import add
from r3 import R3
from tk_drawer import TkDrawer


# 40.0	45.0	-30.0	-60.0 # изменённая строка в run_shadow
class Segment:
    """ Одномерный отрезок """
    # Параметры конструктора: начало и конец отрезка (числа)

    def __init__(self, beg, fin):
        self.beg, self.fin = beg, fin

    def __repr__(self):
        return f"({self.beg}, {self.fin})"

    # Отрезок вырожден?
    def is_degenerate(self):
        return self.beg >= self.fin - 1e-14

    # Пересечение с отрезком
    def intersect(self, other):
        if other.beg > self.beg:
            self.beg = other.beg
        if other.fin < self.fin:
            self.fin = other.fin
        return self

    # Разность отрезков
    # Разность двух отрезков всегда является списком из двух отрезков!
    def subtraction(self, other):
        return [Segment(
            self.beg, self.fin if self.fin < other.beg else other.beg),
            Segment(self.beg if self.beg > other.fin else other.fin, self.fin)]


class Edge:
    """ Ребро полиэдра """
    # Начало и конец стандартного одномерного отрезка
    SBEG, SFIN = 0.0, 1.0

    # Параметры конструктора: начало и конец ребра (точки в R3)
    def __init__(self, beg, fin):
        self.beg, self.fin = beg, fin
        # Список «просветов»
        self.gaps = [Segment(Edge.SBEG, Edge.SFIN)]

    # Учёт тени от одной грани
    def shadow(self, facet):
        # «Вертикальная» грань не затеняет ничего
        if facet.is_vertical():
            return
        # Нахождение одномерной тени на ребре
        shade = Segment(Edge.SBEG, Edge.SFIN)
        for u, v in zip(facet.vertexes, facet.v_normals()):
            shade.intersect(self.intersect_edge_with_normal(u, v))
            if shade.is_degenerate():
                return

        shade.intersect(
            self.intersect_edge_with_normal(
                facet.vertexes[0], facet.h_normal()))
        if shade.is_degenerate():
            return
        # Преобразование списка «просветов», если тень невырождена
        gaps = [s.subtraction(shade) for s in self.gaps]
        self.gaps = [
            s for s in reduce(add, gaps, []) if not s.is_degenerate()]

    # Преобразование одномерных координат в трёхмерные
    def r3(self, t):
        return self.beg * (Edge.SFIN - t) + self.fin * t

    # Пересечение ребра с полупространством, задаваемым точкой (a)
    # на плоскости и вектором внешней нормали (n) к ней
    def intersect_edge_with_normal(self, a, n):
        f0, f1 = n.dot(self.beg - a), n.dot(self.fin - a)
        if f0 >= 0.0 and f1 >= 0.0:
            return Segment(Edge.SFIN, Edge.SBEG)
        if f0 < 0.0 and f1 < 0.0:
            return Segment(Edge.SBEG, Edge.SFIN)
        x = - f0 / (f1 - f0)
        return Segment(Edge.SBEG, x) if f0 < 0.0 else Segment(x, Edge.SFIN)

    # Проверяет видимо ли ребро хотя бы частично
    def is_seen(self):
        if len(self.gaps) == 0:
            return 0
        if self.gaps[0].beg == 0.0 and self.gaps[0].fin == 1.0:
            return 1
        return 1j


class Facet:
    """ Грань полиэдра """
    # Параметры конструктора: список вершин

    def __init__(self, vertexes):
        self.vertexes = vertexes

    # «Вертикальна» ли грань?
    def is_vertical(self):
        return self.h_normal().dot(Polyedr.V) == 0.0

    # Нормаль к «горизонтальному» полупространству
    def h_normal(self):
        n = (
            self.vertexes[1] - self.vertexes[0]).cross(
            self.vertexes[2] - self.vertexes[0])
        return n * (-1.0) if n.dot(Polyedr.V) < 0.0 else n

    # Нормали к «вертикальным» полупространствам, причём k-я из них
    # является нормалью к грани, которая содержит ребро, соединяющее
    # вершины с индексами k-1 и k
    def v_normals(self):
        return [self._vert(x) for x in range(len(self.vertexes))]

    # Вспомогательный метод
    def _vert(self, k):
        n = (self.vertexes[k] - self.vertexes[k - 1]).cross(Polyedr.V)
        return n * \
            (-1.0) if n.dot(self.vertexes[k - 1] - self.center()) < 0.0 else n

    # Центр грани
    def center(self):
        return sum(self.vertexes, R3(0.0, 0.0, 0.0)) * \
            (1.0 / len(self.vertexes))

    def perimeter(self):
        perimeter = 0

        for i in range(len(self.vertexes)-1):
            perimeter += (self.vertexes[i] - self.vertexes[i+1]).abs()
        perimeter += (self.vertexes[0] -
                      self.vertexes[len(self.vertexes)-1]).abs()

        return perimeter


class Polyedr:
    """ Полиэдр """
    # вектор проектирования
    V = R3(0.0, 0.0, 1.0)

    # Параметры конструктора: файл, задающий полиэдр
    def __init__(self, file):

        # списки вершин, рёбер и граней полиэдра
        self.vertexes, self.edges, self.facets = [], [], []
        self.starter_vertexes, self.starter_edges,\
            self.starter_facets = [], [], []

        # список строк файла
        with open(file) as f:
            for i, line in enumerate(f):
                if i == 0:
                    # обрабатываем первую строку; buf - вспомогательный массив
                    buf = line.split()
                    # коэффициент гомотетии
                    c = float(buf.pop(0))
                    # углы Эйлера, определяющие вращение
                    alpha, beta, gamma = (float(x) * pi / 180.0 for x in buf)
                elif i == 1:
                    # во второй строке число вершин, граней и рёбер полиэдра
                    nv, nf, ne = (int(x) for x in line.split())
                elif i < nv + 2:
                    # задание всех вершин полиэдра
                    x, y, z = (float(x) for x in line.split())
                    self.vertexes.append(R3(x, y, z).rz(
                        alpha).ry(beta).rz(gamma) * c)
                    self.starter_vertexes.append(R3(x, y, z))
                else:
                    # вспомогательный массив
                    buf = line.split()
                    # количество вершин очередной грани
                    size = int(buf.pop(0))
                    # массив вершин этой грани
                    vertexes = list(self.vertexes[int(n) - 1] for n in buf)
                    starter_vertexes = list(
                        self.starter_vertexes[int(n) - 1] for n in buf)
                    # задание рёбер грани
                    for n in range(size):
                        self.edges.append(Edge(vertexes[n - 1], vertexes[n]))
                        self.starter_edges.append(
                            Edge(starter_vertexes[n - 1], starter_vertexes[n]))
                    # задание самой грани
                    self.facets.append(Facet(vertexes))
                    self.starter_facets.append(Facet(starter_vertexes))

    def draw(self, tk):
        perimeter = 0
        tk.clean()

        for e in self.edges:
            for f in self.facets:
                e.shadow(f)

        for num, f in enumerate(self.facets):
            print(f"\nStarting examining the facet {num}")
            print("perimeter = ", self.starter_facets[num].perimeter())
            facet_edges = [(f.vertexes[i-1], f.vertexes[i])
                           for i in range(len(f.vertexes))]

            first = None
            ind = True
            for i in self.edges:
                if (i.beg, i.fin) in facet_edges:
                    if first is None:
                        first = i.is_seen()
                    print("is_seen() = ", i.is_seen())
                    if i.is_seen() != first or i.is_seen() == 1j:
                        ind = False
                        print("Facet is partly visible, accepted.")
                        break

            if ind is True:
                continue

            angle = acos(self.starter_facets[num].h_normal().dot(self.V) /
                         self.starter_facets[num].h_normal().abs())
            center_x = self.starter_facets[num].center().x
            center_y = self.starter_facets[num].center().y

            print(f"angle = {angle}, cent_x = {center_x}, cent_y = {center_y}")
            # perimeter += self.starter_facets[num].perimeter()  # test line
            if (angle <= pi/7 and center_x > 1 and center_y > 1):
                print("All conditions were completed")
                perimeter += self.starter_facets[num].perimeter()
            else:
                if (acos(self.starter_facets[num].h_normal().dot(self.V) /
                         self.starter_facets[num].h_normal().abs()) <= pi/7):
                    print("Angle condition was completed")
                if (self.starter_facets[num].center().x > 1 and
                        self.starter_facets[num].center().y > 1):
                    print("Center condition was completed")

        print(f"\nСумма периметров подходящих граней равна {perimeter}")
        for e in self.edges:
            print(e.gaps)
            for s in e.gaps:
                tk.draw_line(e.r3(s.beg), e.r3(s.fin))

    def perimeter(self):
        perimeter = 0

        for e in self.edges:
            for f in self.facets:
                e.shadow(f)

        for num, f in enumerate(self.facets):
            print("perimeter = ", f.perimeter())
            facet_edges = [(f.vertexes[i-1], f.vertexes[i])
                           for i in range(len(f.vertexes))]

            first = None
            ind = True
            for i in self.edges:
                if (i.beg, i.fin) in facet_edges:
                    if first is None:
                        first = i.is_seen()
                    if i.is_seen() != first or i.is_seen() == 1j:
                        ind = False
                        break

            if ind is True:
                continue

            angle = acos(self.starter_facets[num].h_normal().dot(self.V) /
                         self.starter_facets[num].h_normal().abs())
            center_x = self.starter_facets[num].center().x
            center_y = self.starter_facets[num].center().y

            print(f"angle = {angle}, cent_x = {center_x}, cent_y = {center_y}")
            if (angle <= pi/7 and center_x > 1 and center_y > 1):
                perimeter += self.starter_facets[num].perimeter()

        return perimeter


# Testing
if __name__ == "__main__":
    p1, p2 = R3(0, 1, 1), R3(1, 0, 0)
    k1, k2 = R3(5, 4, 4), R3(0, 1, 1)
    print(k2 in [k1, p1])
