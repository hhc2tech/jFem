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
