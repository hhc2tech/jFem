nx=150;ny=180

x=LinRange(0,1,nx)
y=LinRange(0,2,ny)
f=zeros((nx,ny))

for i=1:nx
    for j=1:ny
        f[i,j]=sin(x[i])*cos(y[j])
    end
end


using Plots

# contour(x,y,f,fill=true)
heatmap(x,y,f,aspect_ratio=1,fill=true,dpi=800)
