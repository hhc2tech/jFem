function FindSpan(n,p,u,U)
    #
    if u==U[n+1+1]
        return n
    end
    low=p
    high=n+1
    mid=floor(Int,(low+high)/2.0)
    while u<U[mid+1] || u>=U[mid+2]
        if u<U[mid+1]
            high=mid
        else
            low=mid
        end
        mid=floor(Int,(low+high)/2.0)
    end
    return mid
end


function BasisFun(i,p,u,U)
    # return Ni,p-->p range from 1 to p+1
    N=zeros(1+p)
    N[1]=1.0
    left=zeros(1+p)
    right=zeros(1+p)
    for j=1:p
        left[j+1] =u-U[i+1-j+1]
        right[j+1]=U[i+j+1]-u
        saved=0.0
        for r=0:j-1
            temp=N[r+1]/(right[r+2]+left[j-r+1])
            N[r+1]=saved+right[r+2]*temp
            saved=left[j-r+1]*temp
        end
        N[j+1]=saved
    end
    return N
end


KnotVec=[0,0,0,0,1,1,1,1]
P=3
m=length(KnotVec)
n=m-P-1

X0=LinRange(0.0,1.0*n,n+1)
Y0=rand(n+1)

nx=500
xi=LinRange(0.0,KnotVec[end],nx)
x=zeros(nx)
y=zeros(nx)


for i=1:nx
    u=xi[i]
    uspan=FindSpan(n,P,u,KnotVec)
    # BasisFun(i,p,u,U)
    Nu=BasisFun(uspan,P,u,KnotVec)
    uind=uspan-P
    for j=0:P
        x[i]+=Nu[j+1]*X0[uind+j+1]
        y[i]+=Nu[j+1]*Y0[uind+j+1]
    end
end




using Plots

plot(X0,Y0,label="Control Points")
plot!(x,y,label="NURBS curve")
