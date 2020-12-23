### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# ╔═╡ b26a3f0c-448e-11eb-1cf0-25438afc670d
using Plots

# ╔═╡ 98aeef06-44a5-11eb-1095-adc7773fd362
using LinearAlgebra

# ╔═╡ 9e6d6934-448e-11eb-21f3-69467bda2410
function lineAngle(slope1, slope2)
    a = slope2 - slope1
    b = 1.0 + (slope2*slope1)
    c = abs(a / b)
    d = atan(c)
    return d
end

# ╔═╡ c3fb3d48-448e-11eb-3dff-6fccd893312c
function lineZero(m, b)
    # find root of the passed line
    x = (-b)/m
    y = 0
    return [x, y]
end

# ╔═╡ c3fbc362-448e-11eb-36b3-a5c2694955da
function findFunc(point1, point2)
    # find the linear function that fits two points, returns m, b
    a = point2[2] - point1[2]
    b = point2[1] - point1[1]
    c = (point2[1] * point1[2]) - (point1[1] * point2[2])
    m = a/b
    d = c/b
    return [m, d]
end

# ╔═╡ c3fc5fd4-448e-11eb-2347-9ddee5b2cc8a
function pFunc(m)
    # find the function perpindicular to passed function flope
    if m != 0.0
        mperp = -1.0/m
    end
    if m == 0
        mperp = 0.0
    end

    return mperp
end

# ╔═╡ c40b9c4c-448e-11eb-1cbc-45324d95420d
function newB(m, point)
    # find the y-intersection that the function point pair pass through
    a = m*point[1]
    b = point[2] - a
    return b
end

# ╔═╡ c41174be-448e-11eb-334d-f97165d225d0
"""Find the distance between points p1 = [x, y] and p2"""
function dist(p1, p2)
    # find the distance between two points
    p = p1-p2
    d = norm(p)
    return d
end

# ╔═╡ c418c974-448e-11eb-01b9-8167c73c6e24
"""Find intersection point between lines of form [slope, intercept]"""
function inter(mb1, mb2)
    # find line intersection point between two lines
    a = mb2[2] - mb1[2]
    b = mb1[1] - mb2[1]
    
    if b==0
        return 0, mb1[2]
    end
        
    x = a / b
    y = (mb1[1] * x) + mb1[2]
    return [x, y]
end

# ╔═╡ c41eb1ec-448e-11eb-288d-eb4fecc22f89
"""Find the distance between a line [slope, intercept] and a point [x, y], 
where the slope of the line is an angle in degrees."""
function distpointline(line, point)
    s = tan(deg2rad(line[1]))
    p = pFunc(s)
    perp_line = [p, newB(p, point)]
    # now find intersection of the two lines, where lp is not an angle but usual slope
    lp = (s, line[2])
    # find line intersection point
    intersection = inter(perp_line, lp)
    # determine distance and return
    d = dist(point, intersection)
    return d
end

# ╔═╡ c425aa88-448e-11eb-067f-132bc8887e8a
"""Takes the target lineslope angle, and the given rule step angle,
    and returns the two lines that intersect at the target line at the step angle
"""
function mappingOptions(target_line_angle, step_angle)
    # construct base step angle slopes
    b = tan(deg2rad(step_angle + target_line_angle))
    b2 = tan(deg2rad(target_line_angle - step_angle))
    lines = [b, b2]
    return lines
end

# ╔═╡ 2b64bd10-448f-11eb-361c-1798fccb5fde
# plot the visited points
function plotOrbit(X_m, visitedPoints, points_to_plot, orbit_color)
    x = -800:800
    plot(x, zeros(1601), color="black", xlims=(-20,20), ylims=(-20, 20))

    # now plot lines composing X_m
    @simd for k in X_m
        # define y coordinates
        m = k[1]
        m = tan(deg2rad(m))
        b = k[2]
        plot!(x, (x*m).+b, color="black")
    end

    # now plot the orbit
    tx = []
    for i in visitedPoints[points_to_plot[1]+1:points_to_plot[2]+1]
        push!(tx, i[1])
    end
    ty = []
    for i in visitedPoints[points_to_plot[1]+1:points_to_plot[2]+1]
        push!(ty, i[2])
    end
    # plot
    plot!(tx, ty, color=orbit_color, label="Orbit")
end

# ╔═╡ ef6fb954-448e-11eb-0915-299e31bffd95
"""
The following code takes an n-rule, and plots the path it generates
    in the specified space X_m.  If the n-rule map converges to a periodic orbit
    within num_iter iterations, the program will return the period.
    Inputs:
    X_m: input a list of m tuples, where m is the number of lines plus the x axis, so 
        X_m is of form [(slope1, intercept1), (slope2, intercept2), ...] for line slope and intercept
        where line slope is an angle in degrees.  Further, we require that len(X_m) >= 2
    n-rule: takes a list of n length three tuples of form 
        [(angle1, orientation1, distance1), (angle2, orientation2, distance2), ...]
        where each angle parameter is between 0 and 90, orientation is either 1 or 0, and distance 
        is a positive integer from 1 to m inclusive (dealing with m+1 lines as add the x axis).  If angle=90, orientation can be either 0 or 1 with no change
        of map behavior
    num_iter: number of times to iterate the map
    points_to_plot: an interval [start_point, end_point] that determines the interval of iterations
        of the map to be plotted.
    start_point: list of length 2 [x, y] specifying the x and y coordinates of the point from which to
        begin iteration
    orbit_color: the color of the orbit plot, default is red.
"""
function nRuleMap(X_m, nrule, num_iter, points_to_plot, start_point, orbit_color="red")
    currPoint = start_point
    visitedPoints = [start_point]
    n = length(nrule)
    m = length(X_m) # m-1 as we include x-axis as a line
    
    # begin iteration
    for i in 1:num_iter
        # get rule information
        rulenum = (i%n) + 1
        currRule = nrule[rulenum]
        step_ang = currRule[1]
        o = currRule[2]
        d = currRule[3]
        
        # order liens via distance and choose line at specified distance
        lineIndx = 1
        lineDistOrdering = []
        for j in X_m
            d_j = distpointline(j, currPoint)
            push!(lineDistOrdering, d_j)
        end
        
        # add distance from x axis
        xdist = abs(currPoint[2])
        push!(lineDistOrdering, xdist)
        
        @simd for j in 1:d-1
            # we must find the index of the line in X_m
            # with the dth closest line, so set all other lines
            # to the max distance
            minIndx = argmin(lineDistOrdering)
            lineDistOrdering[minIndx] = Inf
        end
        
        # now determine the index of the line to map to
        lineIndx = argmin(lineDistOrdering)
        # now if index is m+1, then it is the x axis, as we added that to the end
        if lineIndx == m+1
            mappingLine1 = [0, 0]
            mappingLine2 = [0, 0]
		else
			mappingLine1 = X_m[lineIndx]
			s = tan(deg2rad(mappingLine1[1]))
			mappingLine2 = (s, mappingLine1[2])
		end
		
		# upon picking a line, determine two possible step lines via mapping angles
		# and pick one given the orientation
		possibleLines = mappingOptions(mappingLine1[1], step_ang)
		iterationLineSlope = possibleLines[o+1]
		
		# generate the line to find the next point in iteration
		iterationLine = [iterationLineSlope, newB(iterationLineSlope, currPoint)]
		
		# determine the next point in the orbit and add to visited points
		nextPoint = inter(iterationLine, mappingLine2)
		currPoint = nextPoint
		push!(visitedPoints, currPoint)
	end
	
	## Now we plot the lines composing X_m and the iterations of the map
	# first plot the x axis. We plot lines on fairly wide axis
    plotOrbit(X_m, visitedPoints, points_to_plot, orbit_color)
end

# ╔═╡ 0e77beee-449b-11eb-39dd-6f2a6335dee0
"""
A symbolic implementation of n-rule maps. The function symbolicnrulemap takes an n-rule, and plots the path it generates in the specified space X_m.  If the n-rule map converges to a periodic orbit within num_iter iterations, the program will return the period.
Inputs:
X_m: input a list of m tuples, where m is the number of lines plus the x axis, so 
	X_m is of form [(slope1, intercept1), (slope2, intercept2), ...] for line slope 	and intercept where line slope is an angle in degrees.  Further, we require that 	len(X_m) >= 2.  The function is designed so that the index in the list X_m is the symbol value of the corresponding line, and m+1 index is the symbol value of the x-axis.
n-rule: takes a list of n length three tuples of form 
	[(angle1, orientation1, distance1), (angle2, orientation2, distance2), ...]
	where each angle parameter is between 0 and 90, orientation is either 1 or 0, and distance is a positive integer from 1 to m+1 inclusive.  If angle=90, orientation can be either 0 or 1 with no change of map behavior
num_iter: number of times to iterate the map
points_to_plot: an interval [start_point, end_point] that determines the interval of iterations
	of the map to be plotted.
start_point: list of length 2 [x, y] specifying the x and y coordinates of the point from which to
	begin iteration
orbit_color: the color of the orbit plot, default is red.
"""
function symbolicnRuleMap(X_m, nrule, num_iter, points_to_plot, start_point, orbit_color="red")
	
	currPoint = start_point
    visitedPoints = [start_point]
    n = length(nrule)
    m = length(X_m) # m-1 as we include x-axis as a line
    
    # begin iteration
    for i in 1:num_iter
        # get rule information
        rulenum = (i%n) + 1
        currRule = nrule[rulenum]
        step_ang = currRule[1]
        o = currRule[2]
        d = currRule[3]
        
        # determine line
        if d == m+1
            mappingLine1 = [0, 0]
            mappingLine2 = [0, 0]
		else
			mappingLine1 = X_m[d]
			s = tan(deg2rad(mappingLine1[1]))
			mappingLine2 = (s, mappingLine1[2])
		end
		
		# upon picking a line, determine two possible step lines via mapping angles
		# and pick one given the orientation
		possibleLines = mappingOptions(mappingLine1[1], step_ang)
		iterationLineSlope = possibleLines[o+1]
		
		# generate the line to find the next point in iteration
		iterationLine = [iterationLineSlope, newB(iterationLineSlope, currPoint)]
		
		# determine the next point in the orbit and add to visited points
		nextPoint = inter(iterationLine, mappingLine2)
		currPoint = nextPoint
		push!(visitedPoints, currPoint)
	end
	
	## Now we plot the lines composing X_m and the iterations of the map
	# first plot the x axis. We plot lines on fairly wide axis
    plotOrbit(X_m, visitedPoints, points_to_plot, orbit_color)
end

# ╔═╡ Cell order:
# ╠═b26a3f0c-448e-11eb-1cf0-25438afc670d
# ╠═98aeef06-44a5-11eb-1095-adc7773fd362
# ╠═9e6d6934-448e-11eb-21f3-69467bda2410
# ╠═c3fb3d48-448e-11eb-3dff-6fccd893312c
# ╠═c3fbc362-448e-11eb-36b3-a5c2694955da
# ╠═c3fc5fd4-448e-11eb-2347-9ddee5b2cc8a
# ╠═c40b9c4c-448e-11eb-1cbc-45324d95420d
# ╠═c41174be-448e-11eb-334d-f97165d225d0
# ╠═c418c974-448e-11eb-01b9-8167c73c6e24
# ╠═c41eb1ec-448e-11eb-288d-eb4fecc22f89
# ╠═c425aa88-448e-11eb-067f-132bc8887e8a
# ╠═2b64bd10-448f-11eb-361c-1798fccb5fde
# ╠═ef6fb954-448e-11eb-0915-299e31bffd95
# ╠═0e77beee-449b-11eb-39dd-6f2a6335dee0
