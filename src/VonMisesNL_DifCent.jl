#  Dinamica No-Lineal: Cercha Von Mises (Barra de Green) %%%

# Parametros Estructura
rho = 7850.0 # kg/m3 (acero)
Lx  = .374/2
Lz = sqrt(.205^2-Lx^2) # m
l0 = sqrt(Lx^2+Lz^2)   # m
Lc = .240              # m
Ic = .0254*.0032^3/12  # m4
Ac = .0254*.0032       # m2
E  = 200e9             # Pa (steel)
kc = 3*E*Ic/Lc^3       # N/m
mb = l0*Ac*rho         # kg
m  = 1.4               # kg incipient buckling at 1.4
c  = 2                 # kg/s (damping)
g  = 9.81              # m/s2

# Defino Vector de fuerzas Internas: fint(u)  - u = [u1 , u2]^T
function Fint( u )
    return E * Ac * l0 * ( u[1]^2 + u[2]^2 - 2*Lx*u[1] + 2*Lz*u[2] )/2 /1e4 * [ -Lx + u[1], Lz+u[2]] + [ kc*u[1], 0 ]
end

# Defino Vector de fuerzas Externas: gravedad
function Fext( u )
    return [ 0, -(m+mb)/2*g ] # N
end

# Defino Matriz de Masa Concentrada
M = [mb 0 ; 0 (mb+m)/2]

# Defino Matriz de Amortiguamiento
C = [c/10 0 ; 0 c]

# Defino Condiciones Iniciales
t0 = 0
u0 = [0, 0]
v0 = [0, 0]
ac0 = M \ ( Fext(t0) - C*v0 - Fint(u0) )  # de ec de movimiento Mu.. + Fint(u) = ft

# Inicializacion Difrerencia Centrada
tf = 2.0
dt = .000025 # sec

a0 = 1/dt^2
a1 = 1/2/dt
a2 = 2*a0
a3 = 1/a2

numTimes = Int( tf / dt + 3 )

Us = Matrix{Float64}( undef, 2, numTimes )
Vs = Matrix{Float64}( undef, 2, numTimes )
As = Matrix{Float64}( undef, 2, numTimes )

Us[:,1] = u0 - dt*v0 + a3*ac0 # u(-dt)
Us[:,2] = u0 # u(0)

Meff = a0*M + a1*C
M2  = a0*M - a1*C
times = Vector{Float64}(undef, numTimes )
# Comienza Marcha en el Tiempo usando Diferencia Centrada
times[1] = t0-dt
times[2] = t0

k = 2
t = t0

function  epsg(u)
    return ( u[1]^2+2*Lz*u[2] - 2*Lx*u[1] + u[2]^2 )/2/1e2
end

while t < tf
    global k
    global t
print("\nk ",k)
print("t ",t)
    feff = Fext( times[ k ] ) - Fint( Us[ :, k ] ) + a2 * M * Us[ :, k ] - M2 * Us[:, k-1]
    Us[:, k+1 ] = Meff \ feff
    As[:, k+1 ] = a0*( Us[:,k+1] - 2*Us[:,k]  + Us[:,k-1] );
    Vs[:, k+1 ] = a1*( Us[:,k+1] -   Us[:,k-1])
    times[ k ] = t
    k = k + 1
    t = t + dt
end

using Plots
print(size(times))
print(size(Us))
print(size(Us[2,:]))
plot( times, Us[2,:] )

#=
subplot(3,1,1)
plot(t(1:10:end),1000*u(1,1:10:end))
xlabel('t [s]'); ylabel('u_1 [mm]');
axis([0 2 1e3*min(u(1,:))*1.1 1e3*max(u(1,:))*1.1]);
subplot(3,1,2)
plot(t(1:10:end),1000*u(2,1:10:end))
xlabel('t [s]'); ylabel('u_2 [mm]');
axis([0 2 1e3*min(u(2,:))*1.1 1e3*max(u(2,:))*1.1]);
subplot(3,1,3)
plot(t(1:10:end),E*Ac*epsg(1:10:end))
xlabel('t [s]'); ylabel('Directa [N]')
axis([0 2 E*Ac*min(epsg)*1.1 E*Ac*max(epsg)*1.1]);
=#
