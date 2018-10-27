function Ni0(i,u,knot)
    if u>=knot[i] && u<knot[i+1]
        return 1.0
    else
        return 0.0
    end
end

function Ni1(i,p,u,knot)
    term1=(u-knot[i])/(knot[i+p]-knot[i])
    term2=(knot[i+p+1]-u)/(knot[i+p+1]-knot[i+1])
    N1=Ni0(i,u,knot)
    N2=Ni0(i+1,u,knot)
    if knot[i+p]-knot[i]==0.0 && N1==0.0
        # avoid 0/0
        term1=0.0
        N1=0.0
    end
    if knot[i+p+1]-knot[i+1]==0.0 && N2==0.0
        # avoid 0/0
        term2=0.0
        N2=0.0
    end
    return term1*N1+term2*N2
end


function Nip(i,p,u,knot)
    # taken from Eq. 2.5
    if p==0
        return Ni0(i,u,knot)
    elseif p==1
        return Ni1(i,p,u,knot)
    else
        term1=(u-knot[i])/(knot[i+p]-knot[i])
        term2=(knot[i+p+1]-u)/(knot[i+p+1]-knot[i+1])
        N1=Nip(i,p-1,u,knot)
        N2=Nip(i+1,p-1,u,knot)
        if term1==0.0 && N1==0.0
            term1=0.0
            N1=0.0
        end

        if term2==0.0 && N2==0.0
            term2=0.0
            N2=0.0
        end
        return term1*N1+term2*N2
    end
end


KnotVec=[0,0,0,1,2,3,4,4,5,5,5]# length=11
P=2                            # curve order
n=length(KnotVec)-P-1          # number of control points

nx=2000
N=zeros((n,nx))
xi=LinRange(0.0,KnotVec[end],nx)

for i=1:n
    for j=1:nx
        N[i,j]=Nip(i,P,xi[j],KnotVec)
    end
end

using PGFPlots

plot( xi,N[1,:],label=L"N_{1}",dpi=500)
plot!(xi,N[2,:],label=L"N_{2}",dpi=500)
plot!(xi,N[3,:],label=L"N_{3}",dpi=500)
plot!(xi,N[4,:],label=L"N_{4}",dpi=500)
plot!(xi,N[5,:],label=L"N_{5}",dpi=500)
plot!(xi,N[6,:],label=L"N_{6}",dpi=500)
plot!(xi,N[7,:],label=L"N_{7}",dpi=500)
plot!(xi,N[8,:],label=L"N_{8}",dpi=500)
xlabel!(L"\xi")
ylabel!(L"N_{i,p}")
savefig("BSpline.png")
