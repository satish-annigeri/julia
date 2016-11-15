function f(x, y)
  return y + x
end

function euler(f, x0, y0, xn, n)
  x = linspace(x0, xn, n+1)
  y = zeros(x)
  h = x[2] - x[1]
  y[1] = y0
  for i = 1:n
    y[i+1] = y[i] + h * f(x[i], y[i])
  end
  return x, y
end

function heun(f, x0, y0, xn, n)
  x = linspace(x0, xn, n+1)
  y = zeros(x)
  h = x[2] - x[1]
  y[1] = y0
  for i = 1:n
    k1 = h * f(x[i], y[i])
    k2 = h * f(x[i+1], y[i]+k1)
    y[i+1] = y[i] + (k1 + k2) / 2
  end
  return x, y
end

function rk4(f, x0, y0, xn, n)
  x = linspace(x0, xn, n+1)
  y = zeros(x)
  h = x[2] - x[1]
  y[1] = y0
  for i = 1:n
    k1 = h * f(x[i], y[i])
    k2 = h * f(x[i]+h/2, y[i]+k1/2)
    k3 = h * f(x[i]+h/2, y[i]+k2/2)
    k4 = h * f(x[i]+h, y[i]+k3)
    y[i+1] = y[i] + (k1 + 2*k2 + 2*k3 + k4) / 6
  end
  return x, y
end

x, y = euler(f, 0, 0, 1, 5)
println(x)
println(y)

x, y = heun(f, 0, 0, 1, 5)
println(x)
println(y)

x, y = rk4(f, 0, 0, 1, 5)
println(x)
println(y)