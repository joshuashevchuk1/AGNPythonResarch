from numpy import ones, copy, cos, tan, pi, linspace

# quad int params


def guassxw(N):
    a = linspace(3, 4*N-1, N)/(4*N+2)
    x = cos(pi*a+1/(8*N*N*tan(a)))
    epsilon = 1e-15
    delta = 1.0
    while delta > epsilon:
        p0 = ones(N, float)
        p1 = copy(x)
        for k in range(1, N):
            p0, p1 = p1, ((2*k+1)*x*p1-k*p0)/(k+1)
        dp = (N+1)*(p0-x*p1)/(1-x*x)
        dx = p1/dp
        x -= dx
        delta = max(abs(dx))

    w = 2*(N+1)*(N+1)/(N*N*(1-x*x)*dp*dp)

    return x, w


def guassxwb(N, a, b):
    x, w = guassxw(N)
    return 0.5*(b-a)*x+0.5*(b+a), 0.5*(b-a)*w
