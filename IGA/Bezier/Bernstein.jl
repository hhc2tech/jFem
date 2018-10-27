function Bernstein(i,n,u)
    faci=factorial(i)
    facn=factorial(n)
    facn1=factorial(n-i)
    return facn*(u^i)*((1-u)^(n-i))/(faci*facn1)
end

#############################
n=12
X=LinRange(0.0,10.0,n+1)
Y=rand(n+1)

nx=100
u=LinRange(0.0,1.0,nx)

x=zeros(nx)
y=zeros(nx)
for j=1:nx
    x[j]=0.0
    y[j]=0.0
    for i=0:n
        x[j]+=Bernstein(i,n,u[j])*X[i+1]
        y[j]+=Bernstein(i,n,u[j])*Y[i+1]
    end
end

using Plots

plot(x,y,label="BezierCurve",dpi=500)
plot!(X,Y,label="OriginalCurve",dpi=500)
xlabel!("X")
ylabel!("Y")
title!("BezierCurve plot")
savefig("BezierCurve1.png")
