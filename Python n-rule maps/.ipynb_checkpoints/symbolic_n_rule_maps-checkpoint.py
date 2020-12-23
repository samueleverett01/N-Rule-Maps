import math
import matplotlib.pyplot as plt
import random
import numpy as np

def lineAngle(slope1, slope2):
    #for angle between two lines use
    ang = math.atan(abs((slope2-slope1) / (1 + (slope2*slope1))))
    return ang

def lineZero(m, b):
    # find the root of the passed line, want when y=0
    x = (-b)/m
    y = 0
    return x, y

def findFunc(point1, point2):
    """find the linear function that fits two points"""
    m = ((point2[1] - point1[1]) / (point2[0] - point1[0]))
    b = (((point2[0] * point1[1]) - (point1[0] * point2[1])) / (point2[0] - point1 [0]))
    
    return m, b


def pFunc(m):
    """find function perpindicular to passed function slope"""
    if(m != 0):
        mperp = (-1/m)
        
    if(m == 0):
        mperp = 0
        
    return mperp

def newB(m, point):
    """find the y-intersection that the function/point pair passes through"""
    b = (point[1] - (m * point[0]))
    return b

def dist(point1, point2):
    """find the distance between the two passed points"""
    dist = math.sqrt((pow((point2[0] - point1[0]), 2)) + ((pow((point2[1] - point1[1]), 2))))
    return dist

def inter(mb1, mb2):
    """find the point of intersection between two lines"""
    x = ((mb2[1] - mb1[1]) / (mb1[0] - mb2[0]))
    y = ((mb1[0] * x) + mb1[1])
    return x, y

def dist_point_line(line, point):
    """return the distance between a point of form [x, y] and line of form [slope, intercept]"""
    # first find orthogonal line
    s = math.tan(math.radians(line[0]))
    p = pFunc(s)
    perp_line = [p, newB(p, point)]
    # now find intersection of to lines, where lp is converted from angle to usual slope
    lp = (s, line[1])
    # find line intersection point
    intersection = inter(perp_line, lp)
    # determine distance and return
    d = dist(point, intersection)
    return d

def mapping_options(target_line_angle, step_angle):
    """takes the target line slope angle, and the given rule step angle, and returns the two lines
    that intersect at the target line at the step angle"""
    
    # first get base step angle slopes
    base = [math.tan(math.radians(step_angle)), math.tan(math.radians(-step_angle))]
    
    # now construct other lines and give data, in form [line angle 1, line angle 2]
    lines = [math.tan(math.radians(step_angle + target_line_angle)), math.tan(math.radians(target_line_angle - step_angle))]
    return lines

def plot(dpi):
    fig, ax=plt.subplots(dpi=dpi)
    ax.legend()


# find indices of list duplicates
def list_duplicates_of(seq,item):
    start_at = -1
    locs = []
    while True:
        try:
            loc = seq.index(item,start_at+1)
        except ValueError:
            break
        else:
            locs.append(loc)
            start_at = loc
    return locs

def nRuleMap(X_m, nrule, num_iter, points_to_plot, start_point, orbit_color="red"):
    """The following code takes an n-rule, and plots the path it generates
    in the specified space X_m.  If the n-rule map converges to a periodic orbit
    within num_iter iterations, the program will return the period.
    Inputs:
    X_m: input a list of m tuples, where m is the number of lines plus the x axis, so 
        X_m is of form [(slope1, intercept1), (slope2, intercept2), ...] for line slope and intercept
        where line slope is an angle in degrees.  Further, we require that len(X_m) >= 2
    n-rule: takes a list of n length three tuples of form 
        [(angle1, orientation1, distance1), (angle2, orientation2, distance2), ...]
        where each angle parameter is between 0 and 90, orientation is either 1 or 0, and distance 
        is a positive integer from 1 to m inclusive.  If angle=90, orientation can be either 0 or 1 with no change
        of map behavior
    num_iter: number of times to iterate the map
    points_to_plot: an interval [start_point, end_point] that determines the interval of iterations
        of the map to be plotted.
    start_point: list of length 2 [x, y] specifying the x and y coordinates of the point from which to
        begin iteration
    orbit_color: the color of the orbit plot, default is red.
    
    *** In this case of symbolic n-rule maps, the d value of each rule will be taken from [1, m],
    *** where d values of 1 to m-1 will correspond with the line in at the index d of the set X_m,
    *** and if d=m then will map to the x-axis line.
    """
    
    curr_point = start_point
    visited_points = [start_point]
    n = len(nrule)
    m = len(X_m) + 1 # add 1 as include x-axis as a line.
    
    # ensure basic setup is correct    
    assert m >= 2
    
    for i in nrule:
        assert 0 < i[0] <= 90
        assert i[1] == 0 or i[1] == 1
        assert type(i[2]) == int and 1 <= i[2] <= m
    
    # begin iteration
    for i in range(0, num_iter):
        # get rule information
        rule_num = i % n
        curr_rule = nrule[rule_num]
        step_ang = curr_rule[0]
        o = curr_rule[1]
        d = curr_rule[2]
        
        if(d == (m)):
            # then is mapping to x-axis
            mapping_line1 = [0, 0]
            mapping_line2 = [0, 0]
        else:  
            mapping_line1 = X_m[(d-1)]
            s = math.tan(math.radians(mapping_line1[0]))
            mapping_line2 = (s, mapping_line1[1])      
        
        # upon picking line, determine two possible step lines via mapping angles,
        # and pick one given the orientation.
        possible_lines = mapping_options(mapping_line1[0], step_ang)
        iteration_line_slope = possible_lines[o]
        
        # generate the line to find the next point in the iteration of the map
        iteration_line = [iteration_line_slope, newB(iteration_line_slope, curr_point)]

        # determine next point in the orbit, add to visited points
        next_point = inter(iteration_line, mapping_line2)
        curr_point = next_point
        visited_points.append(curr_point)
    
    # now plot the lines in X_m, and the orbit generated by the n-rule map in X_m
    # first plot the m lines
    xLine = np.linspace(-1000, 1000, 1000)
    plt.hlines(0, -1000, 1000, linewidth=2)
    
    for k in X_m:
        s = k[0]
        s = math.tan(math.radians(s))
        b = k[1]
        plt.plot(xLine, (s*xLine + b), color='black', linewidth=2)
    
    # plot the orbit    
    x = []
    for i in visited_points[points_to_plot[0]:points_to_plot[1]]:
        x.append(i[0])
    y = []
    for i in visited_points[points_to_plot[0]:points_to_plot[1]]:
        y.append(i[1])
    
    plt.plot(x, y, color=orbit_color, linewidth=2)
    plt.xlim(-8, 16)
    plt.ylim(-8, 16)
    plt.gca().set_aspect('equal', adjustable='box')
    
    #plot final periodic orbit
    """
    x1 = []
    y1 = []
    
    for i in visited_points[-1000:]:
        x1.append(i[0])

    for i in visited_points[-1000:]:
        y1.append(i[1])
            
    plt.plot(x1, y1, color='blue', linewidth=2)
    
    plt.draw"""
        
    return visited_points
