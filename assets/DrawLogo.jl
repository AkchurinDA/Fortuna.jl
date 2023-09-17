using Luxor

@drawsvg begin
    box(O, 1200, 600, 100, action=:clip)

    @layer begin
        gsave()
        setline(2.5)
        TShape = [Point(x, -250 + 25 * cos(3 * x / 600 * π)) for x in range(-600, +600, 1000)]
        BShape = [Point(x, +250 - 25 * cos(3 * x / 600 * π)) for x in range(-600, +600, 1000)]
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
        setcolor("black")
        setfont("PT Mono", 150)
        settext("FORTUNA", Point(0, 0), halign="center", valign="bottom")
        setfont("PT Mono", 75)
        settext("STRUCTURAL AND SYSTEM", Point(0, 0), halign="center", valign="top")
        settext("RELIABILITY ANALYSIS", Point(0, +75), halign="center", valign="top")
        grestore()
    end

    # @layer begin
    #     gsave()
    #     translate(0, -180)
    #     juliacircles(20)
    #     grestore()
    # end

    setcolor("black")
    box(O, 1200, 600, 100, action=:stroke)
end 1200 600