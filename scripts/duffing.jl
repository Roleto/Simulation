using Pkg
# Pkg.add("Tk")
Pkg.activate(".")
using PyPlot
PyPlot.matplotlib.use("Qt5Agg")

A = 8
α = 1
δ = 0.2
ω = 1
β = 1

LONG=Int(1e5)
dt=1e−3

t=zeros(LONG)
q_pp=zeros(LONG)
q_p=zeros(LONG)
q=zeros(LONG)
q[1]=10
q_p[1]=1.0
#Simulation
for i ∈ 1:(LONG−1)
	t[i]=i*dt
	q_pp[i]=α*q[i]-δ*q_p[i]-β*q[i]^3+A*sin(ω*t[i])
	q_p[i+1]=q_p[i]+dt*q_pp[i]
	q[i+1]=q[i]+dt * q_p[i]
end
#Plotting
#Acceleration
figure(1)
grid(true)
title(L"Gyorsulás")
xlabel(L"Idő")
# xlabel(L"Idő$[s]$")
ylabel(L"Gyorsulás")
# ylabel(L"Gyorsulás$(\ddot{q})$$[\frac{m}{s^2}]$")
plot(t[1:(LONG−1)], q_pp[1:(LONG−1)])
#Velocity
figure(2)
grid(true)
title(L"Sebbeség")
xlabel(L"Idő")
# xlabel(L"Idő$[s]$")
ylabel(L"Sebbeség")
# ylabel(L"Sebbeség$(\dot{q})$$[\frac{m}{s}]$")
plot(t[1:(LONG−1)], q_p[1:(LONG−1)])
#Position
figure(3)
grid(true)
title(L"Helyzet")
# title(L"Helyzet $(q)$")
xlabel(L"Idő")
# xlabel(L"Idő$[s]$")
ylabel(L"Helyzet")
# ylabel(L"Helyzet$(q)$$[m]$")
plot(t[1:(LONG−1)], q[1:(LONG−1)])
#PhaseSpace
figure(4)
grid(true)
title(L"Fázis Tér")
xlabel(L"Helyzet")
# xlabel(L"Helyzet$(q)$$[m]$")
ylabel(L"Sebbeség$")
# ylabel(L"Sebbeség$(\dot{q})$$[\frac{m}{s}]$")
plot(q[1:(LONG−1)], q_p[1:(LONG−1)])

show()
