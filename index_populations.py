"""
This Code  Build a map which associate to a population name of the model  an 
index, in order to create a list that will record the dynamics of the populations.

We first created a dictionary named "index" which will has as keys the names of 
populations and as  vallues the associated index.

The illustration of the model can be found at 

We decided to go row wise to attribute the index corresponding to eah population.

To do so, we first created 

compartments: a list that cointains 3-uplets. each 3-uplets are of
the name of the stage [S, E, P, I, L, D, R] of the disease, the number different 
possibility to be in that stage and number of Erlang stage associated with it. 

upperscript: a list cointaining all the upperscripts used to descript the different
populations. this is list is strategically sorted because some stages don't 
end up with all the different upperscripts. We sorted this list by the frequency
of use of the upperscript. for instance, in our model the upperscript 'NI' 
is almost used in all the stages.  

output desired ==> index : dictionary
"""
from ADE_parameters import *

compartments = [('S', 5, 0), ('E', 5, NE), ('P', 5, NP), ('I', 9, NI), \
            ('L', 12, NL), ('D', 0, 0), ('R', 0, 0)]
upperscripts = ['U', 'NV', 'V', 'NI', 'PI', 'ADE', 'U-', 'U+', 'IV', 'Itilde',\
            'Istar', 'LV', 'Ltilde', 'Lstar']

index = dict()
ind = 0
for i in compartments:
    notation = i[0]
    # in case, there is no upperscripts we just update the index 
    if i[1] == 0 :
        index[notation] = ind
        ind += 1
    else:
        for j in upperscripts[:i[1]+1]:
            # Remenber that the compartments I and L don't have the upperscript (U)
            # we are ignoring to concatenate U to the compartments I and L
            if (i[0] in ['I','L']) and (j in ['U']):
                continue
            notation = i[0] + '^' + j
            # in case, there is no subscripts we just update the index 
            if i[2] == 0 :
                index[notation] = ind
                ind += 1
            else:
                for k in list(range(1,i[2]+1)) :
                    notation = i[0] + '^' + j + '_' + str(k)
                    index[notation] = ind
                    ind += 1
print(index)
