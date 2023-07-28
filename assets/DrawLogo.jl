using Luxor, Colors

# @drawsvg begin
#     # Background:
#     setcolor("antiquewhite")
#     circle(Point(0, 0), 300, action=:fill)

#     # Outer ring:
#     setcolor("goldenrod")
#     circle(Point(0, 0), 200, action=:fill)
#     setcolor("antiquewhite")
#     circle(Point(0, 0), 150, action=:fill)
#     setcolor("black")
#     circle(Point(0, 0), 200, action=:stroke)
#     circle(Point(0, 0), 150, action=:stroke)

#     # Lines:
#     setcolor("black")
#     for i = 1:3
#         line(Point(-150, 0), Point(+150, 0), action=:stroke)
#         rotate(π / 3)
#     end

#     # Reset:
#     origin()

#     # Inner ring:
#     setcolor("goldenrod")
#     circle(Point(0, 0), 50, action=:fill)
#     setcolor("black")
#     circle(Point(0, 0), 50, action=:stroke)

#     # Text:
#     setcolor("black")

#     fontsize(50)
#     text("β", Point(0, 0), halign=:center, valign=:middle)
# end 600 600

@drawsvg begin
    # Draw background:
    setcolor("antiquewhite")
    circle(Point(0, 0), 300, action=:fill)

    # Draw concentric circles
    Radii = range(25, 200, 8)
    JuliaColors = [Luxor.julia_blue, Luxor.julia_red, Luxor.julia_green, Luxor.julia_purple]
    for (i, R) in enumerate(reverse(Radii))
        setcolor(JuliaColors[mod1(i, length(JuliaColors))])
        circle(Point(0, 0), R, action=:fill)
        setcolor("black")
        circle(Point(0, 0), R, action=:stroke)
    end

    # Draw axes:
    setcolor("black")
    line(Point(0, 250), Point(0, -250), action=:stroke)
    line(Point(-250, 0), Point(250, 0), action=:stroke)

    # Add text
    fontsize(50)
    fontface("PTMono-Regular")
    textcurve("FORTUNA", 0, 210, Point(0, 0))
    textcurve("FORTUNA", π, 210, Point(0, 0))
end 600 600