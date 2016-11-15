function gauss_elim(a, b)
    m, n = size(a)
    for p = 1:n-1
        for i = p+1:n
            f = a[i, p] / a[p, p]
            a[i, p:end] = f * a[p, p:end] - a[i, p:end]
            b[i] = f * b[p] - b[i]
        end
    end
    x = zeros(typeof(b[1]), size(b))
    x[n] = b[n] / a[n, n]
    for i = n-1:-1:1
        println("Row: ", i, " ", a[i, i+1:end] * b[i+1:end]
        xx = (b[i] - a[i, i+1:end] * b[i+1:end]) / a[i,i]
        println(typeof(b[i]), " ", typeof(xx), " ", xx)
    end

    println(a)
    return x
end

a = [10.0 3 2; 3.0 8 5; 2 -3 6]
b = [30; 45.0; -10.0]
println(typeof(b))
println("Solution\n", inv(a) * b)
x = gauss_elim(a, b)
println(b)
println(x)
