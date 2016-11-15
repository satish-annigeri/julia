const ecy=0.002
const ecu=0.0035
const Es = 2e5
const gamma_c = 3 // 2
const gamma_s = 115 // 100

function calc_xumax(fy::Float64)
  return ecu / (ecy + ecu + fy / (gamma_s * Es))
end

function calc_xumax(fy::Vector{Float64})
  x = Vector{Number}(fy)
  for (i, ffy) in enumerate(fy)
    x[i] = calc_xumax(ffy)
  end
  return x
end

function k_a(fck=1.0, b=1.0, xu=1.0)
  return 68//189 * fck * b * xu
end

function k_x(xu=1.0)
  return 99//238 * xu
end

function k_m(fy, fck=1.0, b=1.0, d=1.0)
  xumaxbyd = calc_xumax(fy)
  return 68//189 * xumaxbyd * (1.0 - 99//238 * xumaxbyd) * fck * b * d^2
end

function calc_astlim(fck, fy, b, d)
  return 391/945 * fck / fy * calc_xumax(fy) * b * d
end

function calc_xu(fck, b, d, Mu)
  k1 = Mu / (fck * b * d^2)
  k2 = 119 / 99
  return (k2 - sqrt(k2^2 - 147/22*k1)) * d
end

function calc_ast(fck, fy, b, d, Mu)
  xu = calc_xu(fck, b, d, Mu)
  return 391 / 945 * fck / fy * xu * b
end

function calc_fc(ec, fck)
  if ec > ecu || ec < 0.0
    return 0.0
  elseif ec >= ecy
    return fck / gamma_c^2
  else
    ee = ec/ecy
    return (2 * ee - ee^2) * fck / gamma_c^2
  end
end

function calc_fs(es, fy, HYSD=true)
  es = abs(es)
  if !HYSD
    esy = fy / (gamma_s * Es)
    if es >= esy
      return fy / gamma_s
    else
      return es * Es
    end
  else
    y = [0, 0.8, 0.85, 0.9, 0.95, 0.975, 1.0] * fy / gamma_s
    x = y / Es + [0, 0.0, 0.0001, 0.0003, 0.0007, 0.001, 0.002]
    if es >= x[7]
      return y[7]
    elseif es <= x[2]
      return es * Es
    else
      i = 3
      while es > x[i]
        i += 1
      end
      x1 = x[i-1]; y1 = y[i-1]
      x2 = x[i];   y2 = y[i]
      return y1 + (y2 - y1) / (x2 - x1) * (es - x1)
    end
  end
end

function design_rect(fck, fy, b, d, Mu, dc)
  Mulim = k_m(fy, fck, b, d)
  if Mu > Mulim
    # println("Doubly reinforced section")
    xumax = calc_xumax(fy) * d
    Mu2 = Mu - Mulim
    ast1 = calc_astlim(fck, fy, b, d)
    # println("x_umax = $xumax, M_ulim = $Mulim, Mu2 = $Mu2, A_stlim = $ast1")
    ec = ecu * (1.0 - dc / xumax)
    fcc = calc_fc(ec, fck)
    fsc = calc_fs(ec, fy)
    asc = Mu2 / ((fsc - fcc) * (d - dc))
    ast2 = Mu2 / ((fy/gamma_s) * (d - dc))
    ast = ast1 + ast2
    # println("ec = $ec, fcc = $fcc, fsc = $fsc, asc = $asc, ast2 = $ast2")
    return asc, ast
  else
    println("Singly reinforced section")
    ast = calc_ast(fck, fy, b, d, Mu)
    return 0.0, ast
  end
end

function calc_bars(ast, bardia)
  n = ast / (pi * bardia^2 / 4)
  return n
end

function multof(x, m)
  return ceil(x / m) * m
end

function min_hor_dist(bardia, nom_maxsz_ca=20.0, vib=false)
  if vib
    return max(bardia, 2/3*nom_maxsz_ca)
  else
    return max(bardia, nom_maxsz_ca + 5.0)
  end
end

function min_ver_dist(bardia, nom_maxsz_ca=20.0)
  return max(bardia, 15.0, 2/3*nom_maxsz_ca)
end

function calc_layers(b, ast, bardia, endcover, nom_maxsz_ca=20.0, vib=false)
  min_hdist = min_hor_dist(bardia, nom_maxsz_ca, vib)
  nbars = multof(calc_bars(ast, bardia), 1)
  hdist = (b - 2 * endcover - nbars * bardia) / (nbars - 1)
  if hdist >= min_hdist
    return 1
  else
    nbpl = floor((b - 2 * endcover + min_hdist) / (min_hdist + bardia))
    nlayers = multof(nbars / nbpl, 1)
    return nbpl, nlayers, nbars % nbpl
  end
end



# println("e_cy = $ecy, e_cu = $ecu, gamma_c = $gamma_c")
# println("E_s = $Es, gamma_s = $gamma_s")
println("** ", calc_xumax(415.0))
fy = [250.0, 415, 500]
x_umax = calc_xumax(fy)
println("Vector x_umax ", x_umax)
println(k_a(), " ", k_x(), " ", k_m(415.0), " ", k_m(415.0, 20.0, 230.0, 415.0)/1e6)
println(calc_astlim(20.0, 415.0, 230.0, 415.0))
println(calc_xu(20.0, 230.0, 415.0, 100e6))
asc, ast = design_rect(25.0, 415.0, 230.0, 415.0, 190e6, 35.0)
println("Asc = $asc, Ast = $ast")
println(calc_bars(ast, 20), " ", multof(calc_bars(ast, 20), 1))
println(calc_layers(230.0, ast, 20, 25))
