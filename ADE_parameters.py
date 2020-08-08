#parameters
N =
lext =

fsick = 
fdead =
fPI_sick =
fPI_dead =
fADE_sick =
fADE_dead =
fIstar_sick =
fIstar_dead =
fLstar_sick =
fLstar_dead =
fUplus_sick = f_sick / (f_sick + (1 - f_sick) * fUplus_I)

fR_S =
fR_E =
fR_P =
fR_I =
fR_L =

fNI_S =
fPI_S =
fADE_S =

fNI_E =
fPI_E =
fADE_E =

fNI_P =
fPI_P =
fADE_P =

fUplus_I = 
fNI_I =
fPI_I =
fADE_I =
fItilde_I =

fUplus_L =
fNI_L =
fPI_L =
fADE_L =
fItilde_L =
fLtilde_L =

alpha =
nu =

epsilon =
varphi =
gamma =
delta =

phome =
fiso =

Qmax = 
tiso1 = 
tiso2 =

a =
tROmax =
R0bar =
cP =  
cI =  
cL = 
pCont = 

def f_parameter(subscript, upperscript):
    if subscript == 'sick':
        if upperscript in ['', 'V', 'NI']:
            return fsick
        elif upperscript == 'PI'
            return fPI_sick
        elif upperscript == 'ADE'
            return fADE_sick
        elif upperscript == 'U+'
            return fUplus_sick
        elif upperscript == 'Istar'
            return fIstar_sick
        elif upperscript == 'Lstar'
            return fLstar_sick
        else:
            print("Error: enter a correct values for the subscript and upperscript")
    elif subscript == 'dead':
        if upperscript == '':
            return fdead
        elif upperscript == 'PI'
            return fPI_dead
        elif upperscript == 'ADE'
            return fADE_dead
        elif upperscript == 'U+'
            return fUplus_dead
        elif upperscript == 'Istar'
            return fIstar_dead
        elif upperscript == 'Lstar'
            return fLstar_dead
        else:
            print("Error: enter a correct values for the subscript and upperscript")
    elif subscript == 'S'
        if upperscript == 'R':
            return fR_S
        elif upperscript == 'NI':
            return fNI_S
        elif upperscript == 'PI':
            return fPI_S
        elif upperscript == 'ADE':
            return fADE_S
        else:
            print("Error: enter a correct values for the subscript and upperscript")
    elif subscript == 'E'
        if upperscript == 'R':
            return fR_E
        elif upperscript == 'NI':
            return fNI_E
        elif upperscript == 'PI':
            return fPI_E
        elif upperscript == 'ADE':
            return fADE_E
        else:
            print("Error: enter a correct values for the subscript and upperscript")
    elif subscript == 'P'
        if upperscript == 'R':
            return fR_P
        elif upperscript == 'NI':
            return fNI_P
        elif upperscript == 'PI':
            return fPI_P
        elif upperscript == 'ADE':
            return fADE_P
        else:
            print("Error: enter a correct values for the subscript and upperscript")
    elif subscript == 'I'
        if upperscript == 'R':
            return fR_I
        elif upperscript == 'NI':
            return fNI_I
        elif upperscript == 'PI':
            return fPI_I
        elif upperscript == 'ADE':
            return fADE_I
        elif upperscript == 'Uplus':
            return fUplus_I
        elif upperscript == 'Itilde':
            return fItilde_I
        else:
            print("Error: enter a correct values for the subscript and upperscript")
    elif subscript == 'L'
        if upperscript == 'R':
            return fR_L
        elif upperscript == 'NI':
            return fNI_L
        elif upperscript == 'PI':
            return fPI_L
        elif upperscript == 'ADE':
            return fADE_L
        elif upperscript == 'Uplus':
            return fUplus_L
        elif upperscript == 'Itilde':
            return fItilde_L
        elif upperscript == 'Ltilde':
            return fLtilde_L           
        else:
            print("Error: enter a correct values for the subscript and upperscript")
            
            





