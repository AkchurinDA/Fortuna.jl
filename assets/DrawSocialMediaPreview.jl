using Luxor, Colors

H = 640
W = 1280

@draw begin
    # Draw background box
    setcolor("antiquewhite")
    box(O, W, H, action=:fill)

    begin
        # Translate:
        translate(0, -50)

        # Scale:
        scale(340 / 500)

        # Draw concentric circles
        Radii = range(25, 200, 8)
        JuliaColors = [Luxor.julia_blue, Luxor.julia_red, Luxor.julia_green, Luxor.julia_purple]
        for (i, R) in enumerate(reverse(Radii))
            setcolor(JuliaColors[mod1(i, length(JuliaColors))])
            circle(Point(0, 0), R, action=:fill)
            setopacity(1.00)
            setcolor("black")
            circle(Point(0, 0), R, action=:stroke)
        end

        # Draw axes
        line(Point(0, 250), Point(0, -250), action=:stroke)
        line(Point(-250, 0), Point(250, 0), action=:stroke)

        # Add text
        fontsize(70)
        fontface("PTMono-Regular")
        textcurve("STRELA", 0, 210, Point(0, 0))
        textcurve("STRELA", π, 210, Point(0, 0))
    end

    origin()
    fontsize(30)
    fontface("PTMono-Bold")
    text("Structural and System", Point(0, 140), halign=:center, valign=:top)
    text("Reliability Analysis", Point(0, 170), halign=:center, valign=:top)
    text("Made Easy", Point(0, 200), halign=:center, valign=:top)

end W H