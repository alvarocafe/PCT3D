
function gera_CCFace(ELEM,CondCont)
faces=unique(ELEM[:,5])
nfaces=size(faces,1)
CCFace=zeros(nfaces,3)
CCFace[:,1]=faces;
CCFace[:,2]=1.
CCFace[:,3]=0.
ncondcont=size(CondCont,1)
for i=1:ncondcont
  iface=Int32(CondCont[i,1])
  tipo=CondCont[i,2]
  valor=CondCont[i,3]
  CCFace[iface,2:3] = [tipo,valor]
end
return CCFace
end

function gera_vars(ELEM,CCFace,NOS_GEO)

# Descri��o: Gera a matriz CDC a partir das condi��es de contorno das
#            faces.
# Autor:     Gustavo Gontijo
#
# �ltima modifica��o: 24/03/2014 - 13h21m



# Gera a matriz de CDC
# CDC = [n�mero do elemento, tipo da CDC, valor da CDC no n� 1,...
#                               valor da CDC no n� 2, valor da CDC no n� 3]
nelem = size(ELEM,1)
CDC = zeros(nelem,3)
NOS=zeros(nelem,4)
for i=1:nelem
    CDC[i,1] = i
    CDC[i,2] = CCFace[ELEM[i,5],2];
    CDC[i,3] = CCFace[ELEM[i,5],3];
    noselem = ELEM[i,2:4]
    X1=NOS_GEO[noselem[1],2]
    Y1=NOS_GEO[noselem[1],3]
    Z1=NOS_GEO[noselem[1],4]
    X2=NOS_GEO[noselem[2],2]
    Y2=NOS_GEO[noselem[2],3]
    Z2=NOS_GEO[noselem[2],4]
    X3=NOS_GEO[noselem[3],2]
    Y3=NOS_GEO[noselem[3],3]
    Z3=NOS_GEO[noselem[3],4]
    XM=(X1+X2+X3)/3;
    YM=(Y1+Y2+Y3)/3;
    ZM=(Z1+Z2+Z3)/3;
    NOS[i,1]=i
    NOS[i,2]=XM
    NOS[i,3]=YM
    NOS[i,4]=ZM
end
return CDC,NOS
end

function monta_matrizvetor(NOS,NOS_GEO,ELEM,k,CDC)

nelem=size(ELEM,1); # Número de elementos (número de linhas da
#  matriz ELEM)
nnos=nelem; # Número de nós
A=zeros(nnos,nnos); # Inicialização da matriz G
b=zeros(nnos); # Inicialização da matriz H
for pc=1:nelem # Laço sobre os elementos
    tipoCDC = CDC[pc,2]; # Tipo da condição de contorno
    valorCDC = CDC[pc,3]; # Tipo da condição de contorno
    nos = ELEM[pc,2:4]
    X1=NOS_GEO[nos[1],2]
    Y1=NOS_GEO[nos[1],3]
    Z1=NOS_GEO[nos[1],4]
    X2=NOS_GEO[nos[2],2]
    Y2=NOS_GEO[nos[2],3]
    Z2=NOS_GEO[nos[2],4]
    X3=NOS_GEO[nos[3],2]
    Y3=NOS_GEO[nos[3],3]
    Z3=NOS_GEO[nos[3],4]
    l12=(X2-X1);
    l23=(X3-X2);
    m12=(Y2-Y1);
    m23=(Y3-Y2);
    n12=(Z2-Z1);
    n23=(Z3-Z2);
    vet1=[l12 ;m12 ;n12]/norm([l12; m12 ;n12])
    vet3=cross(vet1,[l23 ;m23; n23])
    vet3=vet3/norm(vet3) # Vetor normal unitário ao elemento
    vet2=cross(vet3,vet1)
    for pf=1:nnos # Laço sobre os pontos fonte
        XS=NOS[pf,2] # coordenada x do terceiro nó do quadrilatero desgenerado
        YS=NOS[pf,3] # coordenada y do terceiro nó do quadrilatero desgenerado
        ZS=NOS[pf,4] # coordenada z do terceiro nó do quadrilatero desgenerado
        gg=0; # Inicializa o somatorio de g
        hh=0; # Inicializa o somatorio de h
        d=dot(([X1;Y1;Z1]-[XS;YS;ZS]),vet3)/(norm(vet3)^2)
        o=d*vet3+[XS;YS;ZS] # Coordenadas da origem do sistema local
        arestas=[1 2
        2  3
        3 1] # Número dos nós das 3 arestas do elemento
        NOSElem=[1  X1  Y1  Z1 # Coordenada dos nós que compõe o elemento
            2 X2 Y2 Z2
           3 X3 Y3 Z3]
        zS=abs(dot(vet3,[XS;YS;ZS])-dot(vet3,o)) # Coordenada z do ponto fonte
        for ii=1:3
            p1=arestas[ii,1]
            p2=arestas[ii,2]
            x1=dot(vet1,[NOSElem[p1,2];NOSElem[p1,3];NOSElem[p1,4]])-dot(vet1,o)
            y1=dot(vet2,[NOSElem[p1,2];NOSElem[p1,3];NOSElem[p1,4]])-dot(vet2,o)
            x2=dot(vet1,[NOSElem[p2,2];NOSElem[p2,3];NOSElem[p2,4]])-dot(vet1,o)
            y2=dot(vet2,[NOSElem[p2,2];NOSElem[p2,3];NOSElem[p2,4]])-dot(vet2,o)
            # Comprimento da aresta
            # Cossenos diretores da aresta
            t12=sqrt((x2-x1)^2+(y2-y1)^2)
            l12=(x2-x1)/t12
            m12=(y2-y1)/t12
            # Distância da origem a aresta
            d12=l12*y1-m12*x1
            # Distância do ponto O até os vértices
            t1=l12*x1+m12*y1
            t2=l12*x2+m12*y2
            # Limites de integração (apêndice B do artigo japonês)

            if(d12!=0)
                a1=t1/d12;
                a2=t2/d12;
                h=1/(4*pi)*(-atan2(sqrt(zS^2+(a2^2+1)*d12^2),(a2*zS))+
                    atan2(sqrt(zS^2+(a1^2+1)*d12^2),(a1*zS))+atan2(1,a2)-atan2(1,a1))
                g=-1/(4*pi*k)*(-zS*atan2(sqrt(zS^2+(a2^2+1)*d12^2),(a2*zS))+
                    zS*atan2(sqrt(zS^2+(a1^2+1)*d12^2),(a1*zS))+
                    d12*asinh((a2*d12)/sqrt(zS^2+d12^2))-d12*asinh((a1*d12)/sqrt(zS^2+d12^2))+
                    (atan2(1,a2)-atan2(1,a1))*zS)
            else
                h=0;
                g=0;
            end

            hh=hh+h;
            gg=gg+g;
        end
        if pf == pc
            # Integração singular
            hh=-1/2;
        end
        if tipoCDC == 0
            A[pf,pc] = -gg; # Os elementos de G vão para a matriz A
            b[pf,1] = b[pf,1] - hh*valorCDC# Os elementos de H v�o para o vetor b
        else # O fluxo é conhecido
            A[pf,pc] = + hh; # Os elementos de H vão para a matriz A
            b[pf,1] = b[pf,1] + gg*valorCDC# Os elementos de G vão para o vetor b
        end
    end
end
return A,b
end

function monta_Teq(CDC,x)
# Separa fluxo e temperatura

# ncdc = n�mero de linhas da matriz CDC
# T = vetor que cont�m as temperaturas nos n�s
# q = vetor que cont�m o fluxo nos n�s

ncdc = size(CDC,1)
T=zeros(ncdc)
q=zeros(ncdc)
for i=1:ncdc # La�o sobre as condi��es de contorno
    tipoCDC=CDC[i,2]; # Tipo da condi��o de contorno
    valorCDC=CDC[i,3]; # Valor da condi��o de contorno
    valorcalculado=x[i]; # Valor que antes era desconhecido
    if tipoCDC == 1 # Fluxo � conhecido
        T[i] = valorcalculado; # A temperatura � o valor calculado
        q[i] = valorCDC; # O fluxo � a condi�ao de contorno
    else # A temperatura � conhecida
        T[i] = valorCDC; # A temperatura � a condi�ao de contorno
        q[i] = valorcalculado; # O fluxo � o valor calculado
    end
end
return T,q
end


function mostra_resultados(XYZ,tri,T)
npontos=2
xil = linspace(0, 1, npontos)
eta = linspace(0, 1, npontos)
x=zeros(3)
y=zeros(3)
z=zeros(3)
X=zeros(npontos,npontos)
Y=zeros(npontos,npontos)
Z=zeros(npontos,npontos)
Tmin=minimum(T)
Tmax=maximum(T)
cmp=cm.ScalarMappable(col.Normalize(Tmin,Tmax),cm.jet)
for elem =1:size(tri,1)
    no1=tri[elem,2]
    no2=tri[elem,3]
    no3=tri[elem,4]
    x[1]=XYZ[no1,2]
    y[1]=XYZ[no1,3]
    z[1]=XYZ[no1,4]
    x[2]=XYZ[no2,2]
    y[2]=XYZ[no2,3]
    z[2]=XYZ[no2,4]
    x[3]=XYZ[no3,2]
    y[3]=XYZ[no3,3]
    z[3]=XYZ[no3,4]
    for i = 1: npontos
        for j = 1:npontos
            xi=(1-eta[j])*xil[i];
            N1=xi
            N2=eta[j]
            N3=1-xi-eta[j]
            X[i,j] = N1*x[1]+N2*x[2]+N3*x[3]
            Y[i,j] = N1*y[1]+N2*y[2]+N3*y[3]
            Z[i,j] = N1*z[1]+N2*z[2]+N3*z[3]
        end
    end
    co=cmp[:to_rgba](T[elem])
    plt.plot_surface(X, Y, Z,rstride=2,edgecolors="k", color=co,cstride=2, alpha=1, linewidth=0.25)
end
cmp[:set_array]([Tmin,Tmax])
plt.colorbar(cmp)
xlabel("X")
ylabel("Y")
title("Temperatura")
plt.show()
return
end

function mostra_resultados2(XYZ,tri,T)
nelem=size(tri,1)
x=zeros(3)
y=zeros(3)
z=zeros(3)
zc=zeros(nelem,1)
pc=[zeros(3,3)];
triang=zeros(3,3)
for elem =1:nelem
    no1=tri[elem,2]
    no2=tri[elem,3]
    no3=tri[elem,4]
    x[1]=XYZ[no1,2]
    y[1]=XYZ[no1,3]
    z[1]=XYZ[no1,4]
    x[2]=XYZ[no2,2]
    y[2]=XYZ[no2,3]
    z[2]=XYZ[no2,4]
    x[3]=XYZ[no3,2]
    y[3]=XYZ[no3,3]
    z[3]=XYZ[no3,4]
    triang=[[x[1] y[1] z[1]
    x[2] y[2] z[2]
    x[3] y[3] z[3]]]
    append!(pc,triang)
end
fig = plt.figure()
ax = mp.Axes3D(fig)
q = ar.Poly3DCollection(pc[2:end], linewidths=1,edgecolor="k")
ax[:add_collection3d](q)
m = cm.ScalarMappable(cmap=cm.jet)
b=m[:to_rgba](T[1:nelem])
q[:set_facecolor](b[:,1:3])
m[:set_array]([minimum(T),maximum(T)])
m[:set_clim](vmin=minimum(T),vmax=maximum(T))
plt.colorbar(m, orientation="vertical",shrink=0.9)
ax[:set_xlabel]("x")
ax[:set_ylabel]("y")
ax[:set_zlabel]("z")
xmin=minimum(XYZ[:,2])
ymin=minimum(XYZ[:,3])
zmin=minimum(XYZ[:,4])
xmax=maximum(XYZ[:,2])
ymax=maximum(XYZ[:,3])
zmax=maximum(XYZ[:,4])
deltax=xmax-xmin
deltay=ymax-ymin
deltaz=zmax-zmin
deltamax=maximum([deltax,deltay,deltaz])
ax[:set_xlim3d](xmin-0.2*deltamax,xmin+1.2*deltamax)
ax[:set_ylim3d](ymin-0.2*deltamax,ymin+1.2*deltamax)
ax[:set_zlim3d](zmin-0.2*deltamax,zmin+1.2*deltamax)
ax[:view_init](elev=18., azim=43.)
ax[:axis]("off")
plt.show()
return
end
