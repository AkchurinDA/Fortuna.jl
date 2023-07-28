using Luxor, Colors

@drawsvg begin
    # Box:
    gsave()
    setcolor("navajowhite")
    circle(O, 450, action=:fill)
    grestore()

    # Outer wheel:
    gsave()
    setline(5)
    setcolor("black")
    circle(O, 300, action=:stroke)
    for i = 1:6
        line(Point(0, -100), Point(0, -300), action=:stroke)
        rotate(π / 3)
    end
    circle(O, 200, action=:stroke)
    grestore()

    # Inner wheel:
    gsave()
    setline(5)
    setcolor("goldenrod")
    circle(O, 100, action=:fill)
    setcolor("black")
    circle(O, 100, action=:stroke)
    grestore()

    # Julia cicles:
    gsave()
    rotate(0)
    juliacircles(300, outercircleratio=0.15)
    rotate(π)
    juliacircles(300, outercircleratio=0.15)
    grestore()

    gsave()
    setline(5)
    setcolor("black")
    for i = 1:6
        circle(Point(0, -300), 45, action=:stroke)
        rotate(π / 3)
    end
    grestore()

    # # FORTUNA:
    # gsave()
    # setcolor("black")
    # fontsize(125)
    # fontface("Copperplate-Light")
    # textcurvecentered("FORTUNA", π / 2, 450, O, clockwise=false, letter_spacing=15, baselineshift=25)
    # grestore()
end 900 900