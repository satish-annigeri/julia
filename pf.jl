function pf_calclm(n, zdbcbc)
  lm = zeros(Int64, n, 3)
  ns = size(zdbc)[1]
  for i = 1:ns
    nn = zdbc[i, 1]
    lm[nn, 1:3] = zdbc[i, 2:4]
  end
  nd = 0
  for i = 1:n
    for j = 1:3
      if lm[i, j] == 1
        lm[i, j] = 0
      elseif lm[i, j] == 0
        nd = nd + 1
        lm[i, j] = nd
      end
    end
  end
  nd1 = nd
  for i = 1:n
    for j = 1:3
      if lm[i, j] == -1
        nd += 1
        lm[i, j] = nd
      end
    end
  end
  nd2 = nd
  return lm, nd1, nd2
end

function pf_calclen(imem, xy, conn)
  nodes = conn[imem, :]
  p1 = xy[nodes[1], :]
  p2 = xy[nodes[2], :]
  dxdy = p2 - p1
  L = sqrt(sum(dxdy.^2))
  cx = dxdy[1] / L
  cy = dxdy[2] / L
  return L, cx, cy
end

function pf_calcrot(cx, cy)
  r = zeros(Float64, 6, 6)
  r[1, 1] = cx;  r[1, 2] = cy
  r[2, 1] = -cy; r[2, 2] = cx
  r[3, 3] = 1.0
  r[4:6, 4:6] = r[1:3, 1:3]
  return r
end

function pf_stiff(E, A, I, L)
  k = zeros(Float64, 6, 6)
  scm1 = E * A / L
  scm2 = 4 * E * I / L
  scm3 = 1.5 * scm2 / L
  scm4 = 2 * scm3 / L
  k[1,1] = scm1; k[1,4] = -k[1,1]
  k[2,2] = scm4; k[2,3] = scm3; k[2,5] = -scm4; k[2,6] = scm3
  k[3,3] = scm2; k[3,5] = -scm3; k[3,6] = scm2/2.0
  k[4,4] = scm1
  k[5,5] = scm4; k[5,6] = -scm3
  k[6,6] = scm2
  for i = 2:6
    for j = 1:i-1
      k[i,j] = k[j,i]
    end
  end
  return k
end

function pf_gstiff(imem, xy, conn, mprop)
  L, cx, cy = pf_calclen(imem, xy, conn)
  r = pf_calcrot(cx, cy)
  iprop = conn[imem, 3]
  E = mprop[iprop, 1]
  A = mprop[iprop, 2]
  I = mprop[iprop, 3]
  k = pf_stiff(E, A, I, L)
  K = r' * k * r
  return K
end

function pf_assemssm(imem, xy, conn, mprop, lm, ssm)
  K = pf_gstiff(imem, xy, conn, mprop)
  nj = conn[imem, 1]; nk = conn[imem, 2]
  dof = zeros(Int64, 6)
  dof[1:3] = lm[nj, 1:3]
  dof[4:6] = lm[nk, 1:3]
  for i = 1:6
    ii = dof[i]
    if ii == 0

    else
      for j = 1:6
        jj = dof[j]
        if jj == 0
        else
          ssm[ii, jj] = ssm[ii, jj] + K[i, j]
        end
      end
    end
  end
  return ssm
end

function pf_ssm(xy, conn, mprop, lm, nd)
  SSM = zeros(Float64, nd, nd)
  nmem = size(conn, 1)
  for imem = 1:nmem
    K = pf_assemssm(imem, xy, conn, mprop, lm, SSM)
  end
  return SSM
end

function pf_getdof(imem, conn, lm)
  dof = zeros(Int64, 6)
  n1 = conn[imem, 1]
  n2 = conn[imem, 2]
  dof[1:3] = lm[n1, :]
  dof[4:6] = lm[n2, :]
  return dof
end

function pf_assemloadvec_jl!(lm, jtloads, P)
  njl = size(jtloads, 1)
  for i = 1:njl
    n = jtloads[i, 1]
    dof = lm[n, :]
    for j = 1:3
      jj = dof[j]
      if jj == 0
      else
        P[jj] = P[jj] + jtloads[i, j+1]
      end
    end
  end
  return P
end

function pf_assemloadvec_ml!(iload, xy, conn, lm, memloads, P)
  imem = Int64(memloads[iload, 1])
  L, cx, cy = pf_calclen(imem, xy, conn)
  r = pf_calcrot(cx, cy)
  am = -r' * memloads[iload, 2:7]
  dof = pf_getdof(imem, conn, lm)
  for i = 1:6
    ii = dof[i]
    if ii != 0
      P[ii] = P[ii] + am[i]
    end
  end
  return P
end

function pf_memendforces(imem, xy, conn, mprop, lm, x, memloads)
  iprop = conn[imem, 3]
  E = mprop[iprop, 1]; A = mprop[iprop, 2]; I = mprop[iprop, 3]
  L, cx, cy = pf_calclen(imem, xy, conn)
  r = pf_calcrot(cx, cy)
  k = pf_stiff(E, A, I, L)
  u = zeros(Float64, 6)
  dof = pf_getdof(imem, conn, lm)
  for i = 1:6
    idof = dof[i]
    if idof != 0
      u[i] = x[idof]
    end
  end
  uu = r * u
  f = k * uu

  nml = size(memloads, 1)
  for i = 1:nml
    if memloads[i, 1] == imem
      f = f + memloads[i, 2:7]
    end
  end
  return f
end

function pf(xy, conn, zdbc, mprop, jtloads, memloads)
  n = size(xy, 1)
  lm, nd1, ndof = pf_calclm(n, zdbc)
  K = zeros(Float64, ndof, ndof)
  P = zeros(Float64, ndof)
  x = zeros(Float64, ndof)
  K = pf_ssm(xy, conn, mprop, lm, ndof)
  P = pf_assemloadvec_jl!(lm, jtloads, P)
  nml = size(memloads, 1)
  for iload = 1:nml
    P = pf_assemloadvec_ml!(iload, xy, conn, lm, memloads, P)
  end
  x = K \ P

  nmem = size(conn, 1)
  f = zeros(Float64, nmem, 6)
  println(f)
  for i = 1:nmem
    f[i, :] = pf_memendforces(i, xy, conn, mprop, lm, x, memloads)
  end
  return lm, K, P, x, f
end

xy = [0 0; 0 7; 7 7; 7 3.0]
conn = [1 2 1; 2 3 2; 4 3 1]
zdbc = [1 1 1 1; 4 1 -1 1]
nzdbc = [4 0 0.01 0]
mprop = [2e8 1e-2 3e-4; 2e8 1e-2 6e-4]
jtloads = []
memloads = [1 0 32.1 55.1 0 42.9 -73.4; 2 0 25 43.75 0 25 -43.75]
println("Coordinates")
println(xy)
println("\nMember Connectivity")
println(conn)
println("\nMaterial Properties")
println(mprop)
println("\nZero Displacement Boundary Conditions")
println(zdbc)
println("\nJoint Loads")
println(jtloads)
println("Member Loads")
println(memloads)

n = size(xy, 1)
# lm, cx, cy = pf_calclm(n, zdbc)
lm, K, P, x, f = pf(xy, conn, zdbc, mprop, jtloads, memloads)
println("\nLocation Matrix")
println(lm)
println("\nStructure Stiffness Matrix")
println(K)
println("\nStructure Load Vector")
println(P)
println("\nNode Displacements")
println(x)
println("\nMember End Forces")
for i = 1:size(conn, 1)
  println("Member $i: ", f[i, :])
end
