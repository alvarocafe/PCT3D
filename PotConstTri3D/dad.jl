function dad1()
# Programa: dad1.m

arquivo="cubo.msh"
k = 1.; # Condutividade t�rmica do material

# Condições de contorno das faces:
# CondCont=[número da face, tipo da CDC, valor da CDC]
# Tipo da CDC = 0 => A temperatura é conhecida
# Tipo da CDC = 1 => O fluxo é conhecido
# O número da face é obtido no GMSH
CondCont=[5 0. 0.
6 0. 1.]

return arquivo,CondCont,k
end

function dad2()
    #name of the file
    #file = 'PlacaFuro';
    arquivo= "viga.msh" # Condutividade t�rmica do material

    # Condições de contorno das faces:
    # CondCont=[número da face, tipo da CDC, valor da CDC]
    # Tipo da CDC = 0 => A temperatura é conhecida
    # Tipo da CDC = 1 => O fluxo é conhecido
    # O número da face é obtido no GMSH
    CondCont=[1.  0. 0.
     2. 0. 1.]
    k=1
    return arquivo,CondCont,k
end
function dad3()
    #name of the file
    #file = 'PlacaFuro';
    arquivo= "Placa_furo.msh" # Condutividade t�rmica do material

    # Condições de contorno das faces:
    # CondCont=[número da face, tipo da CDC, valor da CDC]
    # Tipo da CDC = 0 => A temperatura é conhecida
    # Tipo da CDC = 1 => O fluxo é conhecido
    # O número da face é obtido no GMSH
    CondCont=[4. 0. 0.
    2. 0. 1.]
    k=1
    return arquivo,CondCont,k
end
function dad4()
    #name of the file
    #file = 'PlacaFuro';
    arquivo= "L.msh" # Condutividade t�rmica do material

    # Condições de contorno das faces:
    # CondCont=[número da face, tipo da CDC, valor da CDC]
    # Tipo da CDC = 0 => A temperatura é conhecida
    # Tipo da CDC = 1 => O fluxo é conhecido
    # O número da face é obtido no GMSH
    CondCont=[1. 0. 0.
       6. 0. 1.]
    k=1
    return arquivo,CondCont,k
end

function dad5()
    #name of the file
    #file = 'PlacaFuro';
    arquivo= "Tutorial_v4.msh" # Condutividade t�rmica do material

    # Condições de contorno das faces:
    # CondCont=[número da face, tipo da CDC, valor da CDC]
    # Tipo da CDC = 0 => A temperatura é conhecida
    # Tipo da CDC = 1 => O fluxo é conhecido
    # O número da face é obtido no GMSH
    CondCont=[13. 0. 0.
        6. 0. 1.]
    k=1
    return arquivo,CondCont,k

end

function dad6()
    #name of the file
    #file = 'PlacaFuro';
    arquivo= "Polia.msh" # Condutividade t�rmica do material

    # Condições de contorno das faces:
    # CondCont=[número da face, tipo da CDC, valor da CDC]
    # Tipo da CDC = 0 => A temperatura é conhecida
    # Tipo da CDC = 1 => O fluxo é conhecido
    # O número da face é obtido no GMSH
    CondCont=[32. 0. 220.
        3. 0. 0.
        7. 0. 0.]
    k=1
    return arquivo,CondCont,k
end
