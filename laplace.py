import dolfin


def laplace():
    mesh = dolfin.UnitCubeMesh(3, 3, 3)
    left = dolfin.CompiledSubDomain("x[0] <= DOLFIN_EPS")
    right = dolfin.CompiledSubDomain("x[0] >= 1 - DOLFIN_EPS")
    V = dolfin.FunctionSpace(mesh, "CG", 1)

    u = dolfin.TrialFunction(V)
    v = dolfin.TestFunction(V)
    a = dolfin.dot(dolfin.grad(u), dolfin.grad(v)) * dolfin.dx
    L = v * dolfin.Constant(0) * dolfin.dx

    t = dolfin.Function(V)

    endo_bc = dolfin.DirichletBC(V, 0, left)
    epi_bc = dolfin.DirichletBC(V, 1, right)
    bcs = [endo_bc, epi_bc]
    dolfin.solve(
        a == L,
        t,
        bcs,
    )
    return t


if __name__ == "__main__":  # pragma: nocover
    t = laplace()
    dolfin.File("t.pvd") << t
