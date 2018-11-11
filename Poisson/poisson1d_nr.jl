########################################
### For 1d linear poisson equation   ###
########################################
using Printf
using LinearAlgebra

######################################
### Mesh generation                ###
######################################
function Mesh(ne,p,xmin,xmax)
    # Create 1D Lagrange mesh
    if ne<1
        println("Invalid element number(=$ne)!")
        exit()
    end
    if p<1||p>3
        println("Invalid mesh order(=$p)!")
        exit()
    end
    if abs(xmin-xmax)<1.0e-6
        print("Invalid mesh size(xmin=$xmin,xmax=$xmax)")
        exit()
    end
    nNodesPerElmt=p+1
    nElmts=ne
    nNodes=nElmts*p+1
    NodeCoords=zeros(nNodes)
    Conn=zeros(Int,(nElmts,nNodesPerElmt))
    dx=(xmax-xmin)/(nNodes-1)
    for i=1:nNodes
        NodeCoords[i]=xmin+(i-1)*dx
    end
    for e=1:nElmts
        for j=1:nNodesPerElmt
            Conn[e,j]=(e-1)*p+j
        end
    end
    return NodeCoords,Conn
end
###############################################
### 1D shape function(for N and dN)
###############################################
function Shp1D(xi,Coords)
    # For 1D shape function
    nNodes=size(Coords,1)
    if nNodes<2 || nNodes>4
        println("Invalid node number(=$nNodes) for shp1d!")
        exit()
    end
    shp=zeros(nNodes,2)
    if nNodes==2
        # Linear mesh 1+----+2
        shp[1,2]=0.5*(1.0-xi)
        shp[1,1]=-0.5

        shp[2,2]=0.5*(xi+1.0)
        shp[2,1]=0.5
    elseif nNodes==3
        # Quadratic line elements
        # 1+---2---+3
        shp[1,2]=0.5*xi*(xi-1.0)
        shp[1,1]=0.5*(2.0*xi-1.0)

        shp[2,2]=-(xi+1.0)*(xi-1.0)
        shp[2,1]=-2.0*xi

        shp[3,2]=0.5*xi*(xi+1.0)
        shp[3,1]=0.5*(2.0*xi+1.0)
    elseif nNodes==4
        # Third order mesh
        # 1+---2---3---+4
        shp[1,2]=-(3.0*xi+1.0)*(3.0*xi-1.0)*(    xi-1.0)/16.0
        shp[2,2]= (3.0*xi+3.0)*(3.0*xi-1.0)*(3.0*xi-3.0)/16.0
        shp[3,2]=-(3.0*xi+3.0)*(3.0*xi+1.0)*(3.0*xi-3.0)/16.0
        shp[4,2]= (    xi+1.0)*(3.0*xi+1.0)*(3.0*xi-1.0)/16.0
        shp[1,1]=-27.0*xi*xi/16.0+9.0*xi/8.0+ 1.0/16.0
        shp[2,1]= 81.0*xi*xi/16.0-9.0*xi/8.0-27.0/16.0
        shp[3,1]=-81.0*xi*xi/16.0-9.0*xi/8.0+27.0/16.0
        shp[4,1]= 27.0*xi*xi/16.0+9.0*xi/8.0- 1.0/16.0
    end
    dxdxi=0.0
    for i=1:nNodes
        dxdxi+=shp[i,1]*Coords[i]
    end
    DetJac=abs(dxdxi)

    if DetJac<1.0e-13
        println("Singular mesh!")
        exit()
    end
    for i=1:nNodes
        shp[i,1]=shp[i,1]/DetJac # return back to dN/dx
    end
    return shp,DetJac
end
function GaussPoint(n)
    gs=zeros(n,2)
    if n==2
        gs[1,1] =-0.577350269189625764509148780502
        gs[2,1] = 0.577350269189625764509148780502

        gs[1,2] = 1.0
        gs[2,2] = 1.0
    elseif n==3
        gs[1,2] = 0.555555555555555555555555555556
        gs[1,1] =-0.774596669241483377035853079956

        gs[2,2] = 0.888888888888888888888888888889
        gs[2,1] = 0.0

        gs[3,2] = 0.555555555555555555555555555556
        gs[3,1] = 0.774596669241483377035853079956
    end
    return gs
end
##################################################
### Form K matrix and F vector
### div(grad(phi))=c-->use newton-raphson method
### K_ij=-dR^i/dphi^j
### K=(B^j,B^i)--->B=grad(N)
### R=-(grad(phi),N^I)-(c,N^I)
### BC: phi(x=xmin)=PhiL
##      grad(phi(x=xmax))*n=0
##################################################
function FormKR(NodeCoords,Conn,U,c)
    nDofs=size(NodeCoords,1)
    nNodesPerElmt=size(Conn,2)
    K=zeros(nDofs,nDofs)
    RHS=zeros(nDofs)
    nElmts=size(Conn,1)

    localK=zeros(nNodesPerElmt,nNodesPerElmt)
    localR=zeros(nNodesPerElmt)

    ngp=3
    gs=GaussPoint(ngp)

    for e=1:nElmts
        elConn=Conn[e,:]
        Coords=NodeCoords[elConn]
        localU=U[elConn]
        #############################################
        ### Do gauss integration on local element ###
        #############################################
        localK.=0.0
        localR.=0.0
        for gpInd=1:ngp
            xi=gs[gpInd,1]
            shp,xsj=Shp1D(xi,Coords)
            JxW=xsj*gs[gpInd,2]
            # For physical quantities on each gauss point
            gradphi=0.0
            for i=1:nNodesPerElmt
                gradphi+=shp[i,1]*localU[i]
            end
            for i=1:nNodesPerElmt
                localR[i]+=-gradphi*shp[i,1]*JxW-c*shp[i,2]*JxW
                for j=1:nNodesPerElmt
                    localK[i,j]+=shp[j,1]*shp[i,1]*JxW
                end
            end
        end
        #######################################
        ### assemble local array to global
        #######################################
        for i=1:nNodesPerElmt
            iInd=elConn[i]
            RHS[iInd]+=localR[i]
            for j=1:nNodesPerElmt
                jInd=elConn[j]
                K[iInd,jInd]+=localK[i,j]
            end
        end
    end
    return K,RHS
end
################################################
### Apply boundary condition
################################################
function DirichletBC(K,RHS)
    # for the left node(1-st node), phi=PhiL
    # just ensure K11*du1=0==>du1=0, then u1=u1+du1=u1, will never be updated!
    K[1,:].=0.0
    K[:,1].=0.0
    K[1,1]=1.0
    RHS[1]=0.0
    #return K,F
end


########################################
### For the analytical solution
########################################
function Analytical(Coords,x0,x1,PhiL,c)
    A=-c*x1
    B=PhiL-0.5*c*x0^2+c*x0*x1
    nNodes=size(Coords,1)
    phi=zeros(nNodes)
    for i=1:nNodes
        phi[i]=0.5*c*Coords[i]^2+A*Coords[i]+B
    end
    return phi
end


########################################
### Mesh parameters
########################################
ne=1500;xmin=0.0;xmax=1.0;P=3

################################################
### Boundary and governing equation parameters
################################################
PhiL=1.0
c=2.5


@printf("************************************************************\n")
@printf("*** Finite element method for 1d linear poisson equation ***\n")
@printf("***   start mesh generation...                           ***\n")
NodeCoords,Conn=Mesh(ne,P,xmin,xmax)
nElmts=size(Conn,1)
nNodes=size(NodeCoords,1)
nDofs=nNodes
Phi=zeros(nDofs)
Phi[1]=PhiL
@printf("***   mesh generation finished!                          ***\n")
@printf("***   nodes number=%6d, element number=%6d         ***\n",nNodes,nElmts)

@printf("***   start Newton-Raphson iteration...                  ***\n")
MaxIters=50;tol=1.0e-9
iter=0;IsConvergent=false
@time while iter<MaxIters && !IsConvergent
    global Phi
    K,R=FormKR(NodeCoords,Conn,Phi,c)
    DirichletBC(K,R)
    dPhi=K\R

    Rnorm=norm(R)
    dUnorm=norm(dPhi)
    if iter==0
        global R0norm=Rnorm
        global dU0norm=dUnorm
    end

    Phi=Phi+dPhi
    global iter+=1
    @printf("Iter=%3d===>|R0|=%10.4e,|R|=%10.4e,|dU0|=%10.4e,|dU|=%10.4e\n",iter,R0norm,Rnorm,dU0norm,dUnorm)
    if Rnorm<tol
        IsConvergent=true
        break
    end
end


@printf("***   N-R iteration finished!                            ***\n")
@printf("***   start filling analytical solution...               ***\n")
@time  sol=Analytical(NodeCoords,xmin,xmax,PhiL,c)
@printf("***   analytical solution is filled!                     ***\n")
err=norm(Phi-sol)
@printf("Abs Error=%12.6e,Rel Error=%12.6e\n",err,err/norm(sol))
@printf("************************************************************\n")

using Plots
plot(NodeCoords,Phi,label="FEM")
plot!(NodeCoords,sol,label="Analytical")
xlabel!("ϕ")
ylabel!("x")
savefig("res.png")
