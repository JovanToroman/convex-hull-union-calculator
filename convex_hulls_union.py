from __future__ import division
from math import *
from shapely.geometry import LineString
from shapely.geometry import Point
from shapely.geometry import MultiLineString
from shapely.geometry import MultiPoint
from shapely.geometry import Polygon

# distance between two points
s = lambda point1, point2 : atan2(point2[1] - point1[1], point2[0] - point1[0])

# compute convex hull from given points
def h(l):
    r, t, p = [], pi / 2, min(l)
    while 1:
        q = min(set(l) - {p}, key=lambda q: (s(p, q) - t) % (2 * pi));
        m = s(p, q);
        r += [p] * (m != t);
        p = q;
        t = m
        if p in r: return r


def greater(tuple1, tuple2):
    # return true if both components of first tuple are strictly greater, otherwise false
    return tuple1[0]>=tuple2[0] and tuple1[1]>=tuple2[1]


def lower(tuple1, tuple2):
    # return true if both components of first tuple are strictly lower, otherwise false
    return tuple1[0]<=tuple2[0] and tuple1[1]<=tuple2[1]


def lines(points):
    linez = []
    for index, value in enumerate(points):
        if index < len(points) - 1 and points[index] != points[index + 1]:
            linez.append(LineString([points[index], points[index + 1]]))
    return linez


# for two given convex hulls find points of first c.h. that are inside of the other convex hull
def ch_union_find_inside_points(ch1, ch2):
    poly2 = Polygon(ch2)
    inside_points = []
    for i in ch1:
        # if lower(i, max(ch2)) and greater(i, min(ch2)):
        #     inside_points.append(i)
        if (Point(i)).within(poly2):
            inside_points.append(i)
    return inside_points


# find points which are adjacent to the given point in the given convex hull
def find_adjacent_points(ch, inside_point):
    if ch.index(inside_point) == 0:
        return [ch[len(ch) - 1], inside_point, ch[1]]
    elif ch.index(inside_point) == len(ch) - 1:
        return [ch[len(ch) - 2], ch[len(ch) - 1], inside_point]
    else:
        return [ch[ch.index(inside_point) - 1], inside_point, ch[ch.index(inside_point) + 1]]


# find global intersection of two c.h.'s
def global_intersection(ch1,ch2):
    ml1 = MultiLineString([LineString([ch1[len(ch1) - 1], ch1[0]])] + lines(ch1))
    ml2 = MultiLineString([LineString([ch2[len(ch2) - 1], ch2[0]])] + lines(ch2))
    return ml1.intersection(ml2)


# for given two convex hulls as lists of point tuples returns a set of points of union of these two convex hulls
def ch_union(ch1, ch2):
    inside_points_ch1 = ch_union_find_inside_points(ch1, ch2)
    inside_points_ch2 = ch_union_find_inside_points(ch2, ch1)
    # find lines which contain inside points of ch1
    lines1 = []
    lines2 = []
    for point in inside_points_ch1:
        lines1 = lines(find_adjacent_points(ch1, point))
    # find lines which contain inside points of ch2
    for point in inside_points_ch2:
        lines2 = lines(find_adjacent_points(ch2, point))
    # add all points from c. hulls to a set
    all_vertexes = set(ch1)
    all_vertexes.update(ch2)
    # remove all inside points from the set
    all_vertexes -= set(inside_points_ch1)
    all_vertexes -= set(inside_points_ch2)
    # find points of intersection of two convex hulls and add them to the set
    for l1 in lines1:
        for l2 in lines2:
            intersec = l1.intersection(l2)
            if type(intersec) == Point:
                all_vertexes.update([(intersec.x, intersec.y)])
    # find global intersection of two convex hulls
    intersection_ch = global_intersection(ch1,ch2)
    if type(intersection_ch) == MultiLineString:
        # remove points lying on all touching lines
        print("Convex hulls are touching by multiple lines")
        for l in intersection_ch:
            for i in l.xy:
                all_vertexes.discard((tuple(i)))
    elif type(intersection_ch) == LineString:
        print("Convex hulls are touching by single line")
        # remove points lying on the touching line
        for i in intersection_ch.xy:
            all_vertexes.discard((tuple(i)))
    elif type(intersection_ch) == Point or type(intersection_ch) == MultiPoint:
        print("Convex hulls are intersecting")
        if type(intersection_ch) == MultiPoint:
            for p in range(len(intersection_ch)):
                # add points of intersection to all points
                all_vertexes.update([(intersection_ch[p].xy[0][0], intersection_ch[p].xy[1][0])])
    else:
        print("Convex hulls are neither touching nor intersecting")
    return all_vertexes


# test cases
res1 = ch_union([(-2,2), (2,2), (2,-2),(-2,-2)],[(0,0), (4,0), (4,-4), (0,-4)])
print(res1)
assert res1 == {(4, -4), (2, 2), (0, -4), (-2, 2), (2.0, 0.0), (-2, -2), (0.0, -2.0), (4, 0)}
res2 = ch_union([(0,0),(3,0),(2,1.5),(0,1.5)],[(2.5,0),(4,-1),(5,0),(4,1)])
print(res2)
assert res2 == {(0, 1.5), (0, 0), (2, 1.5), (2.5, 0.0), (2.8461538461538463, 0.23076923076923078), (5, 0), (4, 1), (4, -1)}
res3 = ch_union([(-2,2), (2,2), (2,-2),(-2,-2)],[(-2,-2), (2,-2), (2,-6), (-2,-6)])
print(res3)
assert res3 == {(-2, -6), (2, 2), (-2, 2), (2, -6)}
res4 = ch_union([(-3,0),(0,3),(3,0),(0,0),(0,-3)],[(3,0),(0,0),(0,-3),(3,-6),(6,-3)])
print(res4)
assert res4 == {(3, -6), (6, -3), (-3, 0), (0, 3)}
res5 = ch_union([(-3,0),(0,3),(3,0),(4,-4),(0,-3)],[(5,0),(8,-3),(3,-8),(0,-5)])
print(res5)