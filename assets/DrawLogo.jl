using Luxor

@drawsvg begin
    setcolor("black")
    box(O, 900, 600, 100, action=:fill)
    box(O, 900, 600, 100, action=:clip)

    @layer begin
        gsave()
        TShape = [Point(x, -250 + 25 * cos(3 * x / 450 * π)) for x in range(-450, +450, 1000)]
        BShape = [Point(x, +250 - 25 * cos(3 * x / 450 * π)) for x in range(-450, +450, 1000)]
        for i in range(0, 1, 75)
            setcolor(1 - i, i / 4, 1 - i / 4)
            setopacity(0.75)
            Shape = polymorph(BShape, TShape, i,
                easingfunction=easeinoutcubic, closed=(false, false))
            poly(first(Shape), action=:stroke)
        end
        grestore()
    end

    @layer begin
        gsave()
        setcolor("white")
        setfont("PT Mono", 75)
        settext("FORTUNA", Point(0, 0), halign="center", valign="bottom")
        setfont("PT Mono", 25)
        settext("STRUCTURAL AND SYSTEM", Point(0, 0), halign="center", valign="top")
        settext("RELIABILITY ANALYSIS", Point(0, +25), halign="center", valign="top")
        grestore()
    end

    @layer begin
        gsave()
        translate(+200, 0)
        juliacircles(20)
        grestore()

        gsave()
        translate(-200, 0)
        juliacircles(20)
        grestore()
    end
end 900 600