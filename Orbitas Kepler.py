def F(x,y,xd,yd):
    F1=xd
    F2=yd
    F3=-x/(x**2+y**2)**(3/2)
    F4=-y/(x**2+y**2)**(3/2)
    return F1,F2,F3,F4

T=1
N=10
deltat=T/N

x=1
y=0
xd=0
yd=1

for i in range (N):
    F1,F2,F3,F4=F(x,y,xd,yd)
    x=x+deltat*F1
    y=y+deltat*F2
    xd=xd+deltat*F3
    yd=yd+deltat*F4
    print(x,y,xd,yd)