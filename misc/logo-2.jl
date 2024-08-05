using Luxor

NumPoints = 10000

@drawsvg begin
    @layer begin
        for i in 1:-0.05:0
            setcolor(1 - i, 2 * i / 3, 1 - i / 2)

            θ = range(0, 2 * π, NumPoints)
            r = i * (450 .+ 50 * sin.(10 * θ))
            x = r .* cos.(θ)
            y = r .* sin.(θ)
            poly(Point.(x, y), action=:fill)

            rotate(π / 15)
        end
    end

    @layer begin
        setcolor("white")
        setfont("PT Mono Bold", 150)
        settext("FORTUNA", O, halign="center", valign="center")
    end
end 1000 1000