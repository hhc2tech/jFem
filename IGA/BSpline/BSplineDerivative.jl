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

#############################################################
### Now calculate for derivative of BSline basis function
#############################################################
function d1Nip(i,p,u,knot)
    term1=p/(knot[i+p]-knot[i])
    N1=Nip(i,p-1,u,knot)
    term2=p/(knot[i+p+1]-knot[i+1])
    N2=Nip(i+1,p-1,u,knot)
    if knot[i+p]-knot[i]==0.0 && N1==0.0
        term1=0.0
        N1=0.0
    end
    if knot[i+p+1]-knot[i+1]==0.0 && N2==0.0
        term2=0.0
        N2=0.0
    end
    return term1*N1+term2*N2
end

function d2Nip(i,p,u,knot)
    if p<2
        return 0.0
    end

    a00=1.0
    a10=a00/(knot[i+p-1+1]-knot[i])
    if knot[i+p-1+1]-knot[i]==0.0
        a10=0.0
    end
    a20=a10/(knot[i+p-2+1]-knot[i])
    if knot[i+p-2+1]-knot[i]==0.0
        a20=0.0
    end
    ####################################
    a11=-a00/(knot[i+p+1]-knot[i+1])
    if knot[i+p+1]-knot[i+1]==0.0
        a11=0.0
    end
    a21=(a11-a10)/(knot[i+p]-knot[i+1])
    if knot[i+p]-knot[i+1]==0
        a21=0.0
    end
    a22=-a11/(knot[i+p+1]-knot[i+2])
    if knot[i+p+1]-knot[i+2]==0.0
        a22=0.0
    end

    term2=(a11-a10)/(knot[i+p]-knot[i+1])
    N2=Nip(i+1,p-2,u,knot)
    if knot[i+p]-knot[i+1]==0.0 && N2==0.0
        term2=0.0
        N2=0.0
    end

    return p*(p-1)*(a20*Nip(i,p-2,u,knot)+term2*N2+a22*Nip(i+2,p-2))
end



KnotVec=[0,0,0,0,2,4,6,8,8,8,8]# length=11
P=3                            # curve order
n=length(KnotVec)-P-1          # number of control points

nx=2000
N  =zeros((n,nx))
dN =zeros((n,nx))
d2N=zeros((n,nx))
xi=LinRange(0.0,KnotVec[end],nx)

for i=1:n
    for j=1:nx
        N[i,j]  =  Nip(i,P,xi[j],KnotVec)
        dN[i,j] =d1Nip(i,P,xi[j],KnotVec)
        #d2N[i,j]=d2Nip(i,P,xi[j],KnotVec)
    end
end



using PGFPlots

p1=plot( xi,N[1,:],label=L"N_{1}",dpi=500)
plot!(xi,N[2,:],label=L"N_{2}",dpi=500)
plot!(xi,N[3,:],label=L"N_{3}",dpi=500)
plot!(xi,N[4,:],label=L"N_{4}",dpi=500)
plot!(xi,N[5,:],label=L"N_{5}",dpi=500)
plot!(xi,N[6,:],label=L"N_{6}",dpi=500)
plot!(xi,N[7,:],label=L"N_{7}",dpi=500)
xlabel!(L"\xi")
ylabel!(L"N_{i,p}")

p2=plot( xi,dN[1,:],label=L"dN_{2}",dpi=500)
plot!(xi,dN[2,:],label=L"dN_{3}",dpi=500)
plot!(xi,dN[3,:],label=L"dN_{4}",dpi=500)
plot!(xi,dN[4,:],label=L"dN_{5}",dpi=500)
plot!(xi,dN[5,:],label=L"dN_{6}",dpi=500)
plot!(xi,dN[6,:],label=L"dN_{7}",dpi=500)
plot!(xi,dN[7,:],label=L"dN_{7}",dpi=500)
xlabel!(L"\xi")
ylabel!(L"dN_{i,p}")

# p3=plot( xi,d2N[1,:],label=L"dN_{2}",dpi=500)
# plot!(xi,d2N[2,:],label=L"dN_{3}",dpi=500)
# plot!(xi,d2N[3,:],label=L"dN_{4}",dpi=500)
# plot!(xi,d2N[4,:],label=L"dN_{5}",dpi=500)
# plot!(xi,d2N[5,:],label=L"dN_{6}",dpi=500)
# plot!(xi,d2N[6,:],label=L"dN_{7}",dpi=500)
# plot!(xi,d2N[7,:],label=L"dN_{7}",dpi=500)
# xlabel!(L"\xi")
# ylabel!(L"d2N_{i,p}")

plot(p1,p2)
