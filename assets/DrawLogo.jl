using Luxor, Colors

@drawsvg begin
    # Draw background:
    setcolor("antiquewhite")
    circle(Point(0, 0), 300, action=:fill)
    setcolor("black")
    circle(Point(0, 0), 300, action=:stroke)

    # Draw concentric circles
    Radii = range(25, 200, 8)
    JuliaColors = [Luxor.julia_blue, Luxor.julia_red, Luxor.julia_green, Luxor.julia_purple]
    for (i, R) in enumerate(reverse(Radii))
        setcolor(JuliaColors[mod1(i, length(JuliaColors))])
        circle(Point(0, 0), R, action=:fill)
        setcolor("black")
        circle(Point(0, 0), R, action=:stroke)
    end

    # Draw axes
    line(Point(0, 250), Point(0, -250), action=:stroke)
    line(Point(-250, 0), Point(250, 0), action=:stroke)

    # Add text
    fontsize(60)
    fontface("PTMono-Bold")
    textcurve("STRELA", 0, 210, Point(0, 0))
    textcurve("STRELA", π, 210, Point(0, 0))
end 600 600