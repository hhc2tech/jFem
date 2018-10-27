function Bernstein(i,n,u)
    faci=factorial(i)
    facn=factorial(n)
    facn1=factorial(n-i)
    return facn*(u^i)*((1-u)^(n-i))/(faci*facn1)
end

function Rip(i,n,w,u)
    term1=w[i+1]*Bernstein(i,n,u)
    term2=0.0
    for j=0:n
        term2+=w[j+1]*Bernstein(j,n,u)
    end
    return term1/term2
end


n=8
X0=LinRange(0.0,1.0*n,n+1)
Y0=rand(n+1)
Weight=rand(n+1)*2.0

nx=100
u=LinRange(0.0,1.0,nx)
x=zeros(nx)
y=zeros(nx)

for j=1:nx
    x[j]=0.0
    y[j]=0.0
    for i=0:n
        x[j]+=Rip(i,n,Weight,u[j])*X0[i+1]
        y[j]+=Rip(i,n,Weight,u[j])*Y0[i+1]
    end
end

using Plots

plot(x,y,label="RationalBezier",dpi=500)
plot!(X0,Y0,label="ControlPoints",shape=:circle,dpi=500)
xlabel!("x")
ylabel!("y")
savefig("RationalBezier.png")
