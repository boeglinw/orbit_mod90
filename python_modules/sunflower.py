#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 21:20:36 2024


Sunflower algorithm from : https://stackoverflow.com/questions/28567166/uniformly-distribute-x-points-inside-a-circle

This is modified algorithm that adds boundary points

@author: boeglinw
"""

from math import sqrt, sin, cos, pi
phi = (1 + sqrt(5)) / 2  # golden ratio

def sunflower(n, alpha=0, geodesic=False):
    points = []
    angle_stride = 360 * phi if geodesic else 2 * pi / phi ** 2
    b = round(alpha * sqrt(n))  # number of boundary points
    for k in range(1, n + 1):
        r = radius(k, n, b)
        theta = k * angle_stride
        points.append((r * cos(theta), r * sin(theta)))
    return points

def radius(k, n, b):
    if k > n - b:
        return 1.0
    else:
        return sqrt(k - 0.5) / sqrt(n - (b + 1) / 2)


#%% example
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    points = sunflower(500, alpha=2.5, geodesic=False)
    xs = [point[0] for point in points]
    ys = [point[1] for point in points]
    ax.scatter(xs, ys)
    ax.set_aspect('equal') # display as square plot with equal axes
    plt.show()