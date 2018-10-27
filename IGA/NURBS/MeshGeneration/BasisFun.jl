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
