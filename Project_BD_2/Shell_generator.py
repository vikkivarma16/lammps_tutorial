from math import sin, acos, cos
li=[]
v=float(input("Put the radios of spherical icosahedron: "))
b=(2.0*v)/(((2.618033**2.0)+(1.618033**2.0))**(1.0/2.0))
a=1.618033*b
x1=(a+b)/2.0
y1=a/2.0
z1=0.0
li.append([x1,y1,z1])
x2=x1
y2=(-1*y1)
z2=0.0
li.append([x2,y2,z2])
x3=(-1*x1)
y3=y1
z3=0.0
li.append([x3,y3,z3])
x4=(-1*x1)
y4=(-1*y1)
z4=0.0
li.append([x4,y4,z4])
y1=(a+b)/2.0
z1=a/2.0
x1=0.0
li.append([x1,y1,z1])
y2=y1
z2=(-1*z1)
x2=0.0
li.append([x2,y2,z2])
y3=(-1*y1)
z3=z1
x3=0.0
li.append([x3,y3,z3])
y4=(-1*y1)
z4=(-1*z1)
x4=0.0
li.append([x4,y4,z4])
z1=(a+b)/2.0
x1=a/2.0
y1=0.0
li.append([x1,y1,z1])
z2=z1
x2=(-1*x1)
y2=0.0
li.append([x2,y2,z2])
z3=(-1*z1)
x3=x1
y3=0.0
li.append([x3,y3,z3])
z4=(-1*z1)
x4=(-1*x1)
y4=0.0
li.append([x4,y4,z4])
print("The coordinate of the principle vertices is given by:")
print(li)
edge=v/(sin((2.0*3.141)/5.0))
print("The length of the edge is given by:")
print(edge)
ite=0          
drive=[]
fu=open("Icosacordi.txt","w")
radon=0.0
for i in range(len(li)):
    x=li[i][0]
    y=li[i][1]
    z=li[i][2]
    ite=ite+1
    tite=str(ite)
    tx=str(x)
    ty=str(y)
    tz=str(z)
    fu.write(tx)
    fu.write("   ")
    fu.write(ty)
    fu.write("   ")
    fu.write(tz)
    fu.write("\n")
    drive.append([x,y,z])
    don=((li[i][0]**2.0)+(li[i][1]**2.0)+(li[i][2]**2.0))**0.5
    radon=radon+don
radon=radon/len(li)
redon=radon
print("The calculated value of radius is given by:")
print(radon)
mi=[]
for i in range(len(li)):
    for j in range((i+1),len(li)):
        dsi=0
        for n in range(3):
            dsi=((li[i][n]-li[j][n])**2.0)+dsi
        di=dsi**0.5
        dmi=[]
        if di<=(edge+(edge/3.0)):
            dmi.append(li[i])
            dmi.append(li[j])
            mi.append(dmi)
print("The coordinate of the edges is given by:")

leng=len(mi)
print("The total no of edges is given by:")
print(leng)
bum=[]
for i in range(leng):
    for j in range((i+1),leng):
        if mi[i][0]==mi[j][0]:
            for k in range((i+1),leng):
                gum=[]
                if mi[k][0]==mi[j][1] and mi[k][1]==mi[i][1]:
                    gum.append(mi[i][0])
                    gum.append(mi[j][1])
                    gum.append(mi[i][1])
                    bum.append(gum)
                elif mi[k][1]==mi[j][1] and mi[k][0]==mi[i][1]:
                    gum.append(mi[i][0])
                    gum.append(mi[j][1])
                    gum.append(mi[i][1])
                    bum.append(gum)
        elif mi[i][0]==mi[j][1]:
            for k in range((i+1),leng):
                gum=[]
                if mi[j][0]==mi[k][0] and mi[k][1]==mi[i][1]:
                    gum.append(mi[i][0])
                    gum.append(mi[j][0])
                    gum.append(mi[i][1])
                    bum.append(gum)
                elif mi[j][0]==mi[k][1] and mi[k][0]==mi[i][1]:
                    gum.append(mi[i][0])
                    gum.append(mi[j][0])
                    gum.append(mi[i][1])
                    bum.append(gum)
print("The faces is given by:")
print(bum)
print("The number of faces is given by")
print(len(bum))   
que=float(input("Put the no of slice for the edge of icosahedron: "))
for i in range(len(mi)):
    x1=mi[i][1][0]
    y1=mi[i][1][1]
    z1=mi[i][1][2]
    x2=mi[i][0][0]
    y2=mi[i][0][1]
    z2=mi[i][0][2]
    a=(y1*z2)-(y2*z1)
    if abs(a)<=0.000000001:
        a=abs(a)
    b=(x2*z1)-(x1*z2)
    if abs(b)<=0.000000001:
        b=abs(b)
    c=(y2*x1)-(y1*x2)
    if abs(c)<=0.000000001:
        c=abs(c)
    beta=(acos(((x1*x2)+(y1*y2)+(z1*z2))/(radon**2.0)))/que
    gi=beta
    qq=int(que-1.0)
    if a==0.0 and b==0.0:
        for m in range(qq):
                    alpha=cos(beta)*(redon**2.0)
                    zr1=0.0
                    zr2=0.0
                    p=((y1/x1)*(y1/x1))+1.0
                    q=((-2.0*alpha*y1)/(x1**2.0))
                    r=((alpha/x1)**2.0)-(redon**2.0)
                    ant=((q**2.0)-(4.0*p*r))**0.5
                    yr1=(ant-q)/(2.0*p)
                    yr2=(ant+q)/(-2.0*p)
                    xr1=(alpha/x1)-((y1/x1)*yr1)
                    xr2=(alpha/x1)-((y1/x1)*yr2)
                    di1=(((xr1-x2)**2.0)+((yr1-y2)**2.0)+((yr1-y2)**2.0))**0.5
                    di2=(((xr2-x2)**2.0)+((yr2-y2)**2.0)+((yr2-y2)**2.0))**0.5
                    if di1<=di2:
                        x=xr1
                        y=yr1
                        z=zr1
                    else:
                        x=xr2
                        y=yr2
                        z=zr2
                    ite=ite+1.0
                    tite=str(ite)
                    tx=str(x)
                    ty=str(y)
                    tz=str(z)
                    fu.write(tx)
                    fu.write("   ")
                    fu.write(ty)
                    fu.write("   ")
                    fu.write(tz)
                    fu.write("\n")
                    drive.append([x,y,z]) 
                    beta=beta+gi
    elif b==0.0 and c==0.0:
        for m in range(qq):
                    alpha=cos(beta)*(redon**2.0)
                    xr1=0.0
                    xr2=0.0
                    p=((y1/z1)**2.0)+1.0
                    q=((-2.0*alpha*y1)/(z1**2.0))
                    r=((alpha/z1)**2.0)-(redon**2.0)
                    ant=((q**2.0)-(4.0*p*r))**0.5
                    yr1=(ant-q)/(2.0*p)
                    yr2=(ant+q)/(-2.0*p)
                    zr1=(alpha/z1)-((y1/z1)*yr1)
                    zr2=(alpha/z1)-((y1/z1)*yr2)
                    di1=(((xr1-x2)**2.0)+((yr1-y2)**2.0)+((yr1-y2)**2.0))**0.5
                    di2=(((xr2-x2)**2.0)+((yr2-y2)**2.0)+((yr2-y2)**2.0))**0.5
                    if di1<=di2:
                        x=xr1
                        y=yr1
                        z=zr1
                    else:
                        x=xr2
                        y=yr2
                        z=zr2
                    ite=ite+1
                    tite=str(ite)
                    tx=str(x)
                    ty=str(y)
                    tz=str(z)
                    fu.write(tx)
                    fu.write("   ")
                    fu.write(ty)
                    fu.write("   ")
                    fu.write(tz)
                    fu.write("\n")
                    drive.append([x,y,z]) 
                    beta=beta+gi
    elif a==0.0 and c==0.0:
        for m in range(qq):
                    alpha=cos(beta)*(redon**2.0)
                    yr1=0.0
                    yr2=0.0
                    p=((x1/z1)**2.0)+1.0
                    q=((-2.0*alpha*x1)/(z1**2.0))
                    r=((alpha/z1)**2.0)-(redon**2.0)
                    ant=((q**2.0)-(4.0*p*r))**0.5
                    xr1=(ant-q)/(2.0*p)
                    xr2=(ant+q)/(-2.0*p)
                    zr1=(alpha/z1)-((x1/z1)*xr1)
                    zr2=(alpha/z1)-((x1/z1)*xr2)
                    di1=(((xr1-x2)**2.0)+((yr1-y2)**2.0)+((yr1-y2)**2.0))**0.5
                    di2=(((xr2-x2)**2.0)+((yr2-y2)**2.0)+((yr2-y2)**2.0))**0.5
                    if di1<=di2:
                        x=xr1
                        y=yr1
                        z=zr1
                    else:
                        x=xr2
                        y=yr2
                        z=zr2
                    ite=ite+1
                    tite=str(ite)
                    tx=str(x)
                    ty=str(y)
                    tz=str(z)
                    fu.write(tx)
                    fu.write("   ")
                    fu.write(ty)
                    fu.write("   ")
                    fu.write(tz)
                    fu.write("\n")
                    drive.append([x,y,z]) 
                    beta=beta+gi
    elif x1==0.0:
        for m in range(qq):
                    alpha=cos(beta)*(redon**2.0)
                    ca=(alpha/z1)
                    cc=(y1/z1)
                    cb=((c*cc)-b)/a
                    cd=(c*ca)/a
                    p=(cc**2.0)+(cb**2.0)+1.0
                    q=-2.0*((ca*cc)+(cb*cd))
                    r=(ca**2.0)+(cd**2.0)-(redon**2.0)
                    ant=((q**2.0)-(4.0*p*r))**0.5
                    yr1=(ant-q)/(2.0*p)
                    yr2=(ant+q)/(-2.0*p)
                    zr1=ca-(cc*yr1)
                    xr1=(cb*yr1)-cd
                    zr2=ca-(cc*yr2)
                    xr2=(cb*yr2)-cd
                    di1=(((xr1-x2)**2.0)+((yr1-y2)**2.0)+((yr1-y2)**2.0))**0.5
                    di2=(((xr2-x2)**2.0)+((yr2-y2)**2.0)+((yr2-y2)**2.0))**0.5
                    if di1<=di2:
                        x=xr1
                        y=yr1
                        z=zr1
                    else:
                        x=xr2
                        y=yr2
                        z=zr2
                    drive.append([x,y,z])
                    ite=ite+1
                    tite=str(ite)
                    tx=str(x)
                    ty=str(y)
                    tz=str(z)
                    fu.write(tx)
                    fu.write("   ")
                    fu.write(ty)
                    fu.write("   ")
                    fu.write(tz)
                    fu.write("\n")
                    beta=beta+gi                     
    else:
        for m in range(qq):
                    alpha=cos(beta)*(redon**2.0)
                    ca=(c-((a*z1)/x1))/(((y1*a)/x1)-b)
                    cb=((alpha*a)/x1)/(((y1*a)/x1)-b)
                    cc=(alpha-(y1*cb))/x1
                    cd=((y1*ca)+z1)/x1
                    p=(ca**2.0)+(cd**2.0)+1.0
                    q=2.0*((ca*cb)-(cc*cd))
                    r=((cc**2.0)+(cb**2.0))-(redon**2.0)
                    ant=((q**2.0)-(4.0*p*r))**0.5
                    zr1=(ant-q)/(2.0*p)
                    zr2=(ant+q)/(-2.0*p)
                    yr1=(ca*zr1)+cb
                    xr1=cc-(cd*zr1)
                    yr2=(ca*zr2)+cb
                    xr2=cc-(cd*zr2)
                    di1=(((xr1-x2)**2.0)+((yr1-y2)**2.0)+((yr1-y2)**2.0))**0.5
                    di2=(((xr2-x2)**2.0)+((yr2-y2)**2.0)+((yr2-y2)**2.0))**0.5
                    if di1<=di2:
                        x=xr1
                        y=yr1
                        z=zr1
                    else:
                        x=xr2
                        y=yr2
                        z=zr2
                    drive.append([x,y,z])
                    ite=ite+1
                    tite=str(ite)
                    tx=str(x)
                    ty=str(y)
                    tz=str(z)
                    fu.write(tx)
                    fu.write("   ")
                    fu.write(ty)
                    fu.write("   ")
                    fu.write(tz)
                    fu.write("\n")
                    beta=beta+gi
mi=[]
for i in range(len(bum)):
  
    tum=[]
    chum=[]
    x1=bum[i][2][0]
    y1=bum[i][2][1]
    z1=bum[i][2][2]
    x2=bum[i][0][0]
    y2=bum[i][0][1]
    z2=bum[i][0][2]
    a=(y1*z2)-(y2*z1)
    if abs(a)<=0.000000001:
        a=abs(a)
    b=(x2*z1)-(x1*z2)
    if abs(b)<=0.000000001:
        b=abs(b)
    c=(y2*x1)-(y1*x2)
    if abs(c)<=0.000000001:
        c=abs(c)
    beta=(acos(((x1*x2)+(y1*y2)+(z1*z2))/(redon**2.0)))/que
    gi=beta
    qq=int(que-1.0)
    if a==0.0 and b==0.0:
        for m in range(qq):
                    alpha=cos(beta)*(redon**2.0)
                    zr1=0.0
                    zr2=0.0
                    p=((y1/x1)*(y1/x1))+1.0
                    q=((-2.0*alpha*y1)/(x1**2.0))
                    r=((alpha/x1)**2.0)-(redon**2.0)
                    ant=((q**2.0)-(4.0*p*r))**0.5
                    yr1=(ant-q)/(2.0*p)
                    yr2=(ant+q)/(-2.0*p)
                    xr1=(alpha/x1)-((y1/x1)*yr1)
                    xr2=(alpha/x1)-((y1/x1)*yr2)
                    di1=(((xr1-x2)**2.0)+((yr1-y2)**2.0)+((yr1-y2)**2.0))**0.5
                    di2=(((xr2-x2)**2.0)+((yr2-y2)**2.0)+((yr2-y2)**2.0))**0.5
                    if di1<=di2:
                        x=xr1
                        y=yr1
                        z=zr1
                    else:
                        x=xr2
                        y=yr2
                        z=zr2
                    tum.append([x,y,z]) 
                    beta=beta+gi
    elif b==0.0 and c==0.0:
        for m in range(qq):
                    alpha=cos(beta)*(redon**2.0)
                    xr1=0.0
                    xr2=0.0
                    p=((y1/z1)**2.0)+1.0
                    q=((-2.0*alpha*y1)/(z1**2.0))
                    r=((alpha/z1)**2.0)-(redon**2.0)
                    ant=((q**2.0)-(4.0*p*r))**0.5
                    yr1=(ant-q)/(2.0*p)
                    yr2=(ant+q)/(-2.0*p)
                    zr1=(alpha/z1)-((y1/z1)*yr1)
                    zr2=(alpha/z1)-((y1/z1)*yr2)
                    di1=(((xr1-x2)**2.0)+((yr1-y2)**2.0)+((yr1-y2)**2.0))**0.5
                    di2=(((xr2-x2)**2.0)+((yr2-y2)**2.0)+((yr2-y2)**2.0))**0.5
                    if di1<=di2:
                        x=xr1
                        y=yr1
                        z=zr1
                    else:
                        x=xr2
                        y=yr2
                        z=zr2
                    tum.append([x,y,z]) 
                    beta=beta+gi
    elif a==0.0 and c==0.0:
        for m in range(qq):
                    alpha=cos(beta)*(redon**2.0)
                    yr1=0.0
                    yr2=0.0
                    p=((x1/z1)**2.0)+1.0
                    q=((-2.0*alpha*x1)/(z1**2.0))
                    r=((alpha/z1)**2.0)-(redon**2.0)
                    ant=((q**2.0)-(4.0*p*r))**0.5
                    xr1=(ant-q)/(2.0*p)
                    xr2=(ant+q)/(-2.0*p)
                    zr1=(alpha/z1)-((x1/z1)*xr1)
                    zr2=(alpha/z1)-((x1/z1)*xr2)
                    di1=(((xr1-x2)**2.0)+((yr1-y2)**2.0)+((yr1-y2)**2.0))**0.5
                    di2=(((xr2-x2)**2.0)+((yr2-y2)**2.0)+((yr2-y2)**2.0))**0.5
                    if di1<=di2:
                        x=xr1
                        y=yr1
                        z=zr1
                    else:
                        x=xr2
                        y=yr2
                        z=zr2
                    tum.append([x,y,z]) 
                    beta=beta+gi
    elif x1==0.0:
        for m in range(qq):
                    alpha=cos(beta)*(redon**2.0)
                    ca=(alpha/z1)
                    cc=(y1/z1)
                    cb=((c*cc)-b)/a
                    cd=(c*ca)/a
                    p=(cc**2.0)+(cb**2.0)+1.0
                    q=-2.0*((ca*cc)+(cb*cd))
                    r=(ca**2.0)+(cd**2.0)-(redon**2.0)
                    ant=((q**2.0)-(4.0*p*r))**0.5
                    yr1=(ant-q)/(2.0*p)
                    yr2=(ant+q)/(-2.0*p)
                    zr1=ca-(cc*yr1)
                    xr1=(cb*yr1)-cd
                    zr2=ca-(cc*yr2)
                    xr2=(cb*yr2)-cd
                    di1=(((xr1-x2)**2.0)+((yr1-y2)**2.0)+((yr1-y2)**2.0))**0.5
                    di2=(((xr2-x2)**2.0)+((yr2-y2)**2.0)+((yr2-y2)**2.0))**0.5
                    if di1<=di2:
                        x=xr1
                        y=yr1
                        z=zr1
                    else:
                        x=xr2
                        y=yr2
                        z=zr2
                    tum.append([x,y,z])
                    beta=beta+gi                     
    else:
        for m in range(qq):
                    alpha=cos(beta)*(redon**2.0)
                    ca=(c-((a*z1)/x1))/(((y1*a)/x1)-b)
                    cb=((alpha*a)/x1)/(((y1*a)/x1)-b)
                    cc=(alpha-(y1*cb))/x1
                    cd=((y1*ca)+z1)/x1
                    p=(ca**2.0)+(cd**2.0)+1.0
                    q=2.0*((ca*cb)-(cc*cd))
                    r=((cc**2.0)+(cb**2.0))-(redon**2.0)
                    ant=((q**2.0)-(4.0*p*r))**0.5
                    zr1=(ant-q)/(2.0*p)
                    zr2=(ant+q)/(-2.0*p)
                    yr1=(ca*zr1)+cb
                    xr1=cc-(cd*zr1)
                    yr2=(ca*zr2)+cb
                    xr2=cc-(cd*zr2)
                    di1=(((xr1-x2)**2.0)+((yr1-y2)**2.0)+((yr1-y2)**2.0))**0.5
                    di2=(((xr2-x2)**2.0)+((yr2-y2)**2.0)+((yr2-y2)**2.0))**0.5
                    if di1<=di2:
                        x=xr1
                        y=yr1
                        z=zr1
                    else:
                        x=xr2
                        y=yr2
                        z=zr2
                    tum.append([x,y,z])
                    beta=beta+gi
    x1=bum[i][2][0]
    y1=bum[i][2][1]
    z1=bum[i][2][2]
    x2=bum[i][1][0]
    y2=bum[i][1][1]
    z2=bum[i][1][2]
    a=(y1*z2)-(y2*z1)
    if abs(a)<=0.000000001:
        a=abs(a)
    b=(x2*z1)-(x1*z2)
    if abs(b)<=0.000000001:
        b=abs(b)
    c=(y2*x1)-(y1*x2)
    if abs(c)<=0.000000001:
        c=abs(c)
    beta=(acos(((x1*x2)+(y1*y2)+(z1*z2))/(redon**2.0)))/que
    gi=beta
    if a==0.0 and b==0.0:
        for m in range(qq):
                    alpha=cos(beta)*(redon**2.0)
                    zr1=0.0
                    zr2=0.0
                    p=((y1/x1)**2.0)+1.0
                    q=((-2.0*alpha*y1)/(x1**2.0))
                    r=((alpha/x1)**2.0)-(redon**2.0)
                    ant=((q**2.0)-(4.0*p*r))**0.5
                    yr1=(ant-q)/(2.0*p)
                    yr2=(ant+q)/(-2.0*p)
                    xr1=(alpha/x1)-((y1/x1)*yr1)
                    xr2=(alpha/x1)-((y1/x1)*yr2)
                    di1=(((xr1-x2)**2.0)+((yr1-y2)**2.0)+((yr1-y2)**2.0))**0.5
                    di2=(((xr2-x2)**2.0)+((yr2-y2)**2.0)+((yr2-y2)**2.0))**0.5
                    if di1<=di2:
                        x=xr1
                        y=yr1
                        z=zr1
                    else:
                        x=xr2
                        y=yr2
                        z=zr2
                    chum.append([x,y,z]) 
                    beta=beta+gi
    elif b==0.0 and c==0.0:
        for m in range(qq):
                    alpha=cos(beta)*(redon**2.0)
                    xr1=0.0
                    xr2=0.0
                    p=((y1/z1)**2.0)+1.0
                    q=((-2.0*alpha*y1)/(z1**2.0))
                    r=((alpha/z1)**2.0)-(redon**2.0)
                    ant=((q**2.0)-(4.0*p*r))**0.5
                    yr1=(ant-q)/(2.0*p)
                    yr2=(ant+q)/(-2.0*p)
                    zr1=(alpha/z1)-((y1/z1)*yr1)
                    zr2=(alpha/z1)-((y1/z1)*yr2)
                    di1=(((xr1-x2)**2.0)+((yr1-y2)**2.0)+((yr1-y2)**2.0))**0.5
                    di2=(((xr2-x2)**2.0)+((yr2-y2)**2.0)+((yr2-y2)**2.0))**0.5
                    if di1<=di2:
                        x=xr1
                        y=yr1
                        z=zr1
                    else:
                        x=xr2
                        y=yr2
                        z=zr2
                    chum.append([x,y,z]) 
                    beta=beta+gi
    elif a==0.0 and c==0.0:
        for m in range(qq):
                    alpha=cos(beta)*(redon**2.0)
                    yr1=0.0
                    yr2=0.0
                    p=((x1/z1)**2.0)+1.0
                    q=((-2.0*alpha*x1)/(z1**2.0))
                    r=((alpha/z1)**2.0)-(redon**2.0)
                    ant=((q**2.0)-(4.0*p*r))**0.5
                    xr1=(ant-q)/(2.0*p)
                    xr2=(ant+q)/(-2.0*p)
                    zr1=(alpha/z1)-((x1/z1)*xr1)
                    zr2=(alpha/z1)-((x1/z1)*xr2)
                    di1=(((xr1-x2)**2.0)+((yr1-y2)**2.0)+((yr1-y2)**2.0))**0.5
                    di2=(((xr2-x2)**2.0)+((yr2-y2)**2.0)+((yr2-y2)**2.0))**0.5
                    if di1<=di2:
                        x=xr1
                        y=yr1
                        z=zr1
                    else:
                        x=xr2
                        y=yr2
                        z=zr2
                    chum.append([x,y,z]) 
                    beta=beta+gi
    elif x1==0.0:
        for m in range(qq):
                    alpha=cos(beta)*(redon**2.0)
                    ca=(alpha/z1)
                    cc=(y1/z1)
                    cb=((c*cc)-b)/a
                    cd=(c*ca)/a
                    p=(cc**2.0)+(cb**2.0)+1.0
                    q=-2.0*((ca*cc)+(cb*cd))
                    r=(ca**2.0)+(cd**2.0)-(redon**2.0)
                    ant=((q**2.0)-(4.0*p*r))**0.5
                    yr1=(ant-q)/(2.0*p)
                    yr2=(ant+q)/(-2.0*p)
                    zr1=ca-(cc*yr1)
                    xr1=(cb*yr1)-cd
                    zr2=ca-(cc*yr2)
                    xr2=(cb*yr2)-cd
                    di1=(((xr1-x2)**2.0)+((yr1-y2)**2.0)+((yr1-y2)**2.0))**0.5
                    di2=(((xr2-x2)**2.0)+((yr2-y2)**2.0)+((yr2-y2)**2.0))**0.5
                    if di1<=di2:
                        x=xr1
                        y=yr1
                        z=zr1
                    else:
                        x=xr2
                        y=yr2
                        z=zr2
                    chum.append([x,y,z])
                    beta=beta+gi
    else:
        for m in range(qq):
                    alpha=cos(beta)*(redon**2.0)
                    ca=(c-((a*z1)/x1))/(((y1*a)/x1)-b)
                    cb=((alpha*a)/x1)/(((y1*a)/x1)-b)
                    cc=(alpha-(y1*cb))/x1
                    cd=((y1*ca)+z1)/x1
                    p=(ca**2.0)+(cd**2.0)+1.0
                    q=2.0*((ca*cb)-(cc*cd))
                    r=((cc**2.0)+(cb**2.0))-(redon**2.0)
                    ant=((q**2.0)-(4.0*p*r))**0.5
                    zr1=(ant-q)/(2.0*p)
                    zr2=(ant+q)/(-2.0*p)
                    yr1=(ca*zr1)+cb
                    xr1=cc-(cd*zr1)
                    yr2=(ca*zr2)+cb
                    xr2=cc-(cd*zr2)
                    di1=(((xr1-x2)**2.0)+((yr1-y2)**2.0)+((yr1-y2)**2.0))**0.5
                    di2=(((xr2-x2)**2.0)+((yr2-y2)**2.0)+((yr2-y2)**2.0))**0.5
                    if di1<=di2:
                        x=xr1
                        y=yr1
                        z=zr1
                    else:
                        x=xr2
                        y=yr2
                        z=zr2
                    chum.append([x,y,z]) 
                    beta=beta+gi   
    divi=1.0
    for k in range((len(chum))-1): 
        j=k+1
        divi=divi+1.0
        x1=tum[j][0]
        if abs(x1)<=0.0000000001:
                    x1=0.0
        y1=tum[j][1]
        if abs(y1)<=0.0000000001:
                    y1=0.0
        z1=tum[j][2]
        if abs(z1)<=0.0000000001:
                    z1=0.0
        x2=chum[j][0]
        if abs(x2)<=0.0000000001:
                    x2=0.0
        y2=chum[j][1]
        if abs(y2)<=0.0000000001:
                    y2=0.0
        z2=chum[j][2]
        if abs(z2)<=0.0000000001:
                    z2=0.0
        a=(y1*z2)-(y2*z1)
        if abs(a)<=0.000000001:
            a=abs(a)
        b=(x2*z1)-(x1*z2)
        if abs(b)<=0.000000001:
            b=abs(b)
        c=(y2*x1)-(y1*x2)
        if abs(c)<=0.000000001:
            c=abs(c)
        beta=(acos(((x1*x2)+(y1*y2)+(z1*z2))/(redon**2.0)))/divi
        n=int((divi-1.0))
        gi=beta
        if a==0.0 and b==0.0:
            for m in range(n):
                    alpha=cos(beta)*(redon**2.0)
                    zr1=0.0
                    zr2=0.0
                    p=((y1/x1)**2.0)+1.0
                    q=((-2.0*alpha*y1)/(x1**2.0))
                    r=((alpha/x1)**2.0)-(redon**2.0)
                    ant=((q**2.0)-(4.0*p*r))**0.5
                    yr1=(ant-q)/(2.0*p)
                    yr2=(ant+q)/(-2.0*p)
                    xr1=(alpha/x1)-((y1/x1)*yr1)
                    xr2=(alpha/x1)-((y1/x1)*yr2)
                    di1=(((xr1-x2)**2.0)+((yr1-y2)**2.0)+((yr1-y2)**2.0))**0.5
                    di2=(((xr2-x2)**2.0)+((yr2-y2)**2.0)+((yr2-y2)**2.0))**0.5
                    if di1<=di2:
                        x=xr1
                        y=yr1
                        z=zr1
                    else:
                        x=xr2
                        y=yr2
                        z=zr2
                    ite=ite+1
                    tite=str(ite)
                    tx=str(x)
                    ty=str(y)
                    tz=str(z)
                    fu.write(tx)
                    fu.write("   ")
                    fu.write(ty)
                    fu.write("   ")
                    fu.write(tz)
                    fu.write("\n")
                    drive.append([x,y,z]) 
                    beta=beta+gi
        elif b==0.0 and c==0.0:
            for m in range(n):
                    alpha=cos(beta)*(redon**2.0)
                    xr1=0.0
                    xr2=0.0
                    p=((y1/z1)**2.0)+1.0
                    q=((-2.0*alpha*y1)/(z1**2.0))
                    r=((alpha/z1)**2.0)-(redon**2.0)
                    ant=((q**2.0)-(4.0*p*r))**0.5
                    yr1=(ant-q)/(2.0*p)
                    yr2=(ant+q)/(-2.0*p)
                    zr1=(alpha/z1)-((y1/z1)*yr1)
                    zr2=(alpha/z1)-((y1/z1)*yr2)
                    di1=(((xr1-x2)**2.0)+((yr1-y2)**2.0)+((yr1-y2)**2.0))**0.5
                    di2=(((xr2-x2)**2.0)+((yr2-y2)**2.0)+((yr2-y2)**2.0))**0.5
                    if di1<=di2:
                        x=xr1
                        y=yr1
                        z=zr1
                    else:
                        x=xr2
                        y=yr2
                        z=zr2
                    ite=ite+1
                    tite=str(ite)
                    tx=str(x)
                    ty=str(y)
                    tz=str(z)
                    fu.write(tx)
                    fu.write("   ")
                    fu.write(ty)
                    fu.write("   ")
                    fu.write(tz)
                    fu.write("\n")
                    drive.append([x,y,z]) 
                    beta=beta+gi
        elif a==0.0 and c==0.0:
            for m in range(n):
                    alpha=cos(beta)*(redon**2.0)
                    yr1=0.0
                    yr2=0.0
                    p=((x1/z1)**2.0)+1.0
                    q=((-2.0*alpha*x1)/(z1**2.0))
                    r=((alpha/z1)**2.0)-(redon**2.0)
                    ant=((q**2.0)-(4.0*p*r))**0.5
                    xr1=(ant-q)/(2.0*p)
                    xr2=(ant+q)/(-2.0*p)
                    zr1=(alpha/z1)-((x1/z1)*xr1)
                    zr2=(alpha/z1)-((x1/z1)*xr2)
                    di1=(((xr1-x2)**2.0)+((yr1-y2)**2.0)+((yr1-y2)**2.0))**0.5
                    di2=(((xr2-x2)**2.0)+((yr2-y2)**2.0)+((yr2-y2)**2.0))**0.5
                    if di1<=di2:
                        x=xr1
                        y=yr1
                        z=zr1
                    else:
                        x=xr2
                        y=yr2
                        z=zr2
                    ite=ite+1
                    tite=str(ite)
                    tx=str(x)
                    ty=str(y)
                    tz=str(z)
                    fu.write(tx)
                    fu.write("   ")
                    fu.write(ty)
                    fu.write("   ")
                    fu.write(tz)
                    fu.write("\n")
                    drive.append([x,y,z]) 
                    beta=beta+gi
        elif x1==0.0:
            for m in range(n):
                    alpha=cos(beta)*(redon**2.0)
                    ca=(alpha/z1)
                    cc=(y1/z1)
                    cb=((c*cc)-b)/a
                    cd=(c*ca)/a
                    p=(cc**2.0)+(cb**2.0)+1.0
                    q=-2.0*((ca*cc)+(cb*cd))
                    r=(ca**2.0)+(cd**2.0)-(redon**2.0)
                    ant=((q**2.0)-(4.0*p*r))**0.5
                    yr1=(ant-q)/(2.0*p)
                    yr2=(ant+q)/(-2.0*p)
                    zr1=ca-(cc*yr1)
                    xr1=(cb*yr1)-cd
                    zr2=ca-(cc*yr2)
                    xr2=(cb*yr2)-cd
                    di1=(((xr1-x2)**2.0)+((yr1-y2)**2.0)+((yr1-y2)**2.0))**0.5
                    di2=(((xr2-x2)**2.0)+((yr2-y2)**2.0)+((yr2-y2)**2.0))**0.5
                    if di1<=di2:
                        x=xr1
                        y=yr1
                        z=zr1
                    else:
                        x=xr2
                        y=yr2
                        z=zr2
                    drive.append([x,y,z])
                    ite=ite+1
                    tite=str(ite)
                    tx=str(x)
                    ty=str(y)
                    tz=str(z)
                    fu.write(tx)
                    fu.write("   ")
                    fu.write(ty)
                    fu.write("   ")
                    fu.write(tz)
                    fu.write("\n")
                    beta=beta+gi            
        else:
            for m in range(n):
                    alpha=cos(beta)*(redon**2.0)
                    ca=(c-((a*z1)/x1))/(((y1*a)/x1)-b)
                    cb=((alpha*a)/x1)/(((y1*a)/x1)-b)
                    cc=(alpha-(y1*cb))/x1
                    cd=((y1*ca)+z1)/x1
                    p=(ca**2.0)+(cd**2.0)+1.0
                    q=2.0*((ca*cb)-(cc*cd))
                    r=((cc**2.0)+(cb**2.0))-(redon**2.0)
                    ant=((q**2.0)-(4.0*p*r))**0.5
                    zr1=(ant-q)/(2.0*p)
                    zr2=(ant+q)/(-2.0*p)
                    yr1=(ca*zr1)+cb
                    xr1=cc-(cd*zr1)
                    yr2=(ca*zr2)+cb
                    xr2=cc-(cd*zr2)
                    di1=(((xr1-x2)**2.0)+((yr1-y2)**2.0)+((yr1-y2)**2.0))**0.5
                    di2=(((xr2-x2)**2.0)+((yr2-y2)**2.0)+((yr2-y2)**2.0))**0.5
                    if di1<=di2:
                        x=xr1
                        y=yr1
                        z=zr1
                    else:
                        x=xr2
                        y=yr2
                        z=zr2
                    drive.append([x,y,z])
                    ite=ite+1
                    tite=str(ite)
                    tx=str(x)
                    ty=str(y)
                    tz=str(z)
                    fu.write(tx)
                    fu.write("   ")
                    fu.write(ty)
                    fu.write("   ")
                    fu.write(tz)
                    fu.write("\n") 
                    beta=beta+gi
fu.close()





oradius = float(input("Please enter the opening radius of the spherical icosahedron (should be less than the radius of the shells): "))

z_boundary = (v * v - oradius * oradius) ** 0.5  # Compute z-boundary

su = open("shell_template.mol", "w")  
fu = open("Icosacordi.txt", "r")  

drive1 = []  

atom_id = 1  # Start atom numbering from 1
valid_particles = []  # List to store particles satisfying the condition

line = fu.readline()  

while line:  
    ju = line.split()  
    chi = float(ju[0])  
    khi = float(ju[1])  
    fri = float(ju[2])  

    # Check if the particle satisfies the z-boundary condition
    if fri <= z_boundary:
        valid_particles.append((atom_id, chi, khi, fri))  # Store valid particle data
        drive1.append([chi, khi, fri])  # Store for counting  
        atom_id += 1  # Increment only for valid particles  

    line = fu.readline()  

fu.close()  

# Writing to the output file
su.write(f"\n{len(valid_particles)} atoms\n\n")  # Write valid particle count
su.write("Coords\n\n")  # Section title

for atom in valid_particles:
    su.write(f"{atom[0]} {atom[1]:.4f} {atom[2]:.4f} {atom[3]:.4f}\n")  

su.write("\nTypes\n\n")  # Section title for atom types  

for atom in valid_particles:
    atom_type = 1  # Assigning type 1 to all atoms
    su.write(f"{atom[0]} {atom_type}\n")  

su.close()  

print(len(drive1), "valid atoms written.")  


    
            
            
