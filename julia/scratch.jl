function testy(x; name="plot")
    println("I am running $(plot)!")
    return x/2
end

function grabby(y, z = testy(y); name="plot")
    println("I am running grabby!")
    println("grabby( $(y) , $(z) )")
    return y+z
end
