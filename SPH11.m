dx=1/1000;  %����x^2-x+1(0<x<1)
h=dx*5;
ad=1.25/h;
X=0:dx:1;
Y=linspace(0,0,1001);
L=linspace(0,0,1001);
for i=50:950
    x=X(i);
    y=0;y1=0;
    for j=1:1000
        x1=0+j*dx;
        R=(x-x1)/h;
        m=-1;
        if R<0
            m=1;
            R=-R;
        end
        if R<=1
            y=y+ad*(1+3*R )*(1-R)^3*(x1^2-x1+1)*dx;
            y1=y1-ad*(1-R)^2*(-12*R)*(x1^2-x1+1)*m*dx/h;
        end
    end
        Y(i)=y;
        L(i)=y1;
   
end
figure
plot(X,Y,X,L)