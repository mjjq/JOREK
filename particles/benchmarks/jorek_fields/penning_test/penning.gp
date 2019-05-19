omega_e = 4.9
omega_b = 25.
epsilon = -1.
x0x = 10.
x0y = 0.
x0z = 0.
v0x = 50.
v0y = 0.
v0z = 20.

i = {0.,1.}

omega = sqrt(-2.*epsilon)*omega_e
omega_plus = 0.5*(omega_b + sqrt(omega_b**2 + 4.*epsilon*omega_e**2))
omega_minus = 0.5*(omega_b - sqrt(omega_b**2 + 4.*epsilon*omega_e**2))
R_minus = (omega_plus * x0x + v0y)/(omega_plus - omega_minus)
R_plus  = x0x - R_minus
T_minus = (omega_plus * x0y - v0x)/(omega_plus - omega_minus)
T_plus  = x0y - T_minus

# Calculate the result in the x-y plane in terms of w = x + iy
w(t) = (R_plus + i*T_plus)*exp(-i*omega_plus*t) + (R_minus+i*T_minus)*exp(-i*omega_minus*t)
x(t) = real(w(t))
y(t) = imag(w(t))
z(t) = x0z*cos(omega*t)+v0z*sin(omega*t)/omega
