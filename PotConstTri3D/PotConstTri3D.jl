# Programa de elementos de contorno aplicado a problemas de condu��o de
# calor tri-dimensional sem gera��o de calor
# Tipo de elementos: Triangulares e constantes
# Brasília, março de 2017
# Por: Éder Lima de Albuquerque

# cd C:/Eder/Work/Julia/PotConstTri3D
include("bem_functions.jl")
include("dad.jl")
include("lermsh.jl")
using PyPlot
using PyCall
# @pyimport matplotlib.colors as col
@pyimport matplotlib.cm as cm
@pyimport mpl_toolkits.mplot3d as mp
@pyimport mpl_toolkits.mplot3d.art3d as ar

plt=PyPlot

PyPlot.close("all")
# NOS_GEO,ELEM,CCFace,k=dad1() # Arquivo de entrada de dados
# dadGID # Arquivo de entrada de dados
arquivo,CondCont,k=dad5()

NOS_GEO,ELEM,elemint,elemtipo=lermsh(arquivo,3)
CCFace=gera_CCFace(ELEM,CondCont)
nelem=size(ELEM,1)
println("Número de nós: $nelem")

CDC,NOS = gera_vars(ELEM,CCFace,NOS_GEO);   # Gera a matriz de condições de contorno
# Mostra a geometria do problema

println("1. Montando as matrizes")
tic()
A, b = monta_matrizvetor(NOS, NOS_GEO, ELEM,k,CDC); # Monta as matrizes H e G
tempo1=toq();


println("2. Resolvendo o sistema")
tic()
x = A\b; # Calcula as vari�veis desconhecidas
tempo2=toq();

println("3. Reordenando as variaveis")
tic()
T,q= monta_Teq(CDC,x); # Monta T e q
tempo3=toq();# Mostra o mapa de cor
println("4. Criando o mapa de cor")
tic()
mostra_resultados2(NOS_GEO,ELEM,T)
tempo4=toq()

println("Tempo para montar as matrizes = $tempo1")
println("Tempo para resolver o sistema = $tempo2")
println("Tempo para reordenar as variaveis = $tempo3")
println("Tempo para plotar os resultados = $tempo4")
# salva(arquivo,T,elemtipo)
