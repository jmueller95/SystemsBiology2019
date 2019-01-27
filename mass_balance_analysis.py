# SBMLSqueezer was run with the default parameters except:
# - Set "Ueberschreibe vorhandene..." to True, because model already contained (dummy) kinetic laws
# - Set "Maximale Anzahl an Edukten" to 10 (maximum number of educts might be bigger than 3 and it can't hurt)

import tellurium as te
import libsbml
import numpy as np
import pandas as pd
from sympy import Matrix
from pprint import PrettyPrinter

asansm_sbml = "fbc/iPAE1146_Amino_sugar_and_nucleotide_sugar_metabolism_fbc_squeezed.xml"
pm_sbml = "fbc/iPAE1146_Pyrimidine_metabolism_fbc_squeezed.xml"
lb_sbml = "fbc/iPAE1146_Lipopolysaccharide_biosynthesis_fbc_squeezed.xml"
combined_sbml = "fbc/iPAE1146_Combined_Subsystems_fbc_squeezed.xml"


def get_reaction_names(listOfReactions, reactionIDs):
    return [listOfReactions.get(id).getName() for id in reactionIDs]


def get_metabolite_names(listOfMetabolites, metaboliteIDs):
    return [listOfMetabolites.get(id).getName() for id in metaboliteIDs]


def nullspace(M):
    # Calculate nullspace
    # Copied and modified from https://tellurium.readthedocs.io/en/latest/_modules/tellurium/utils/matrix.html#nullspace
    atol = 1e-13
    rtol = 0
    u, s, vh = np.linalg.svd(M)
    tol = max(atol, rtol * s[0])
    nnz = (s >= tol).sum()
    ns = vh[nnz:].conj().T
    return ns


def find_conservative_relations(S_matrix):
    S_matrix_transpose = np.transpose(S_matrix)
    S_matrix_left_null_space = Matrix(S_matrix_transpose).nullspace()

    number_of_conservative_relations = len(S_matrix_left_null_space)
    conservative_relationship_dict = {}

    for i in range(number_of_conservative_relations):
        counter = 0
        conservative_relationship_dict[i] = ''
        for metabolite in S_matrix.rownames:
            if S_matrix_left_null_space[i][counter] != 0:
                conservative_relationship_dict[i] += ' +('
                conservative_relationship_dict[i] += str(int(S_matrix_left_null_space[i][counter]))
                conservative_relationship_dict[i] += ')"'
                conservative_relationship_dict[i] += metabolite
                conservative_relationship_dict[i] += '"'
            counter += 1
        conservative_relationship_dict[i] += ' = 0'
    return conservative_relationship_dict


def write_stoichiometric_matrix_to_file(M, filepath):
    M_df = pd.DataFrame(data=M, index=M.rownames, columns=M.colnames, dtype=int)
    with open(filepath, "w") as outfile:
        outfile.write(M_df.to_csv(sep="\t"))


def write_left_nullspace_to_file(M, filepath, species_names):  # pass rownames of corresponding S as species_names
    M_df = pd.DataFrame([[round(x, 4) for x in line] for line in M], columns=species_names)
    with open(filepath, "w") as outfile:
        outfile.write(M_df.to_csv(sep="\t"))


def write_eigenvalues_to_file(model, filepath, species_names):
    M_df = pd.DataFrame([model.getFullEigenValues()], columns=species_names)
    with open(filepath, "w") as outfile:
        outfile.write(M_df.to_csv(sep="\t"))


# We start with the smallest subsystem: Amino Sugar and Nucleotide Sugar Metabolism
asansm_libsbml_doc = libsbml.readSBML(asansm_sbml)
asansm_libsbml = asansm_libsbml_doc.getModel()
asansm_model = te.loadSBMLModel(asansm_sbml)
S_asansm = asansm_model.getFullStoichiometryMatrix()
S_asansm.colnames = get_reaction_names(asansm_libsbml.getListOfReactions(), S_asansm.colnames)
S_asansm.rownames = get_metabolite_names(asansm_libsbml.getListOfSpecies(), S_asansm.rownames)
asansm_right_nullspace = Matrix(S_asansm).nullspace()
asansm_left_nullspace = Matrix(S_asansm.T).nullspace()
asansm_rank = Matrix(S_asansm).rank()

# Next subsystem: Pyrimidine Metabolism
# This one is a little special because Matrix().nullspace computes an empty right nullspace here, but there is one!
# So we do it the pedestrian way...
pm_libsbml_doc = libsbml.readSBML(pm_sbml)
pm_libsbml = pm_libsbml_doc.getModel()
pm_model = te.loadSBMLModel(pm_sbml)
S_pm = pm_model.getFullStoichiometryMatrix()
S_pm.colnames = get_reaction_names(pm_libsbml.getListOfReactions(), S_pm.colnames)
S_pm.rownames = get_metabolite_names(pm_libsbml.getListOfSpecies(), S_pm.rownames)
pm_right_nullspace_numpy = nullspace(S_pm)
pm_right_nullspace_numpy = [round(elem / min(pm_right_nullspace_numpy)) for elem in pm_right_nullspace_numpy]
pm_right_nullspace_sympy = Matrix(S_pm).nullspace()
pm_left_nullspace = Matrix(S_pm.T).nullspace()
pm_rank_numpy = np.linalg.matrix_rank(S_pm)
pm_rank_sympy = Matrix(S_pm).rank()

# Third one: Lipopolysaccharide Biosynthesis
lb_libsbml_doc = libsbml.readSBML(lb_sbml)
lb_libsbml = lb_libsbml_doc.getModel()
lb_model = te.loadSBMLModel(lb_sbml)
S_lb = lb_model.getFullStoichiometryMatrix()
S_lb.colnames = get_reaction_names(lb_libsbml.getListOfReactions(), S_lb.colnames)
S_lb.rownames = get_metabolite_names(lb_libsbml.getListOfSpecies(), S_lb.rownames)
lb_right_nullspace = Matrix(S_lb).nullspace()
lb_left_nullspace = Matrix(S_lb.T).nullspace()
lb_rank = Matrix(S_lb).rank()

# Combined one
combined_libsbml_doc = libsbml.readSBML(combined_sbml)
combined_libsbml = combined_libsbml_doc.getModel()
combined_model = te.loadSBMLModel(combined_sbml)
S_combined = combined_model.getFullStoichiometryMatrix()
S_combined.colnames = get_reaction_names(combined_libsbml.getListOfReactions(), S_combined.colnames)
S_combined.rownames = get_metabolite_names(combined_libsbml.getListOfSpecies(), S_combined.rownames)
combined_right_nullspace = Matrix(S_combined).nullspace()
combined_left_nullspace = Matrix(S_combined.T).nullspace()
combined_rank = Matrix(S_combined).rank()


#Write all the tables to files 
write_stoichiometric_matrix_to_file(S_asansm, "matrices/S_asansm.tsv")
write_left_nullspace_to_file(asansm_left_nullspace, "matrices/asansm_left_nullspace.tsv", S_asansm.rownames)
write_stoichiometric_matrix_to_file(S_pm, "matrices/S_pm.tsv")
write_left_nullspace_to_file(pm_left_nullspace, "matrices/pm_left_nullspace.tsv", S_pm.rownames)
# Writing to file is also a little different in PM subsystem due to the different calculation of the nullspace
with open("matrices/pm_nullspace.tsv", "w") as outfile:
    pm_nullspace_df = pd.DataFrame([pm_right_nullspace_numpy], columns=S_pm.colnames)
    outfile.write(pm_nullspace_df.to_csv(sep="\t"))
write_stoichiometric_matrix_to_file(S_lb, "matrices/S_lb.tsv")
write_left_nullspace_to_file(lb_left_nullspace, "matrices/lb_left_nullspace.tsv", S_lb.rownames)
write_stoichiometric_matrix_to_file(S_combined, "matrices/S_combined.tsv")
write_left_nullspace_to_file(combined_left_nullspace, "matrices/combined_left_nullspace.tsv", S_combined.rownames)
write_eigenvalues_to_file(asansm_model, "matrices/asansm_eigenvalues.tsv", S_asansm.rownames)
write_eigenvalues_to_file(lb_model, "matrices/lb_eigenvalues.tsv", S_lb.rownames)
write_eigenvalues_to_file(pm_model, "matrices/pm_eigenvalues.tsv", S_pm.rownames)
write_eigenvalues_to_file(combined_model, "matrices/combined_eigenvalues.tsv", S_combined.rownames)

#Print ranks and conservation relations
pp = PrettyPrinter()

print("The stoichiometric matrix of Amino Sugar and Nucleotide Sugar Metabolism has rank " + str(
    asansm_rank) + ".")
print("The stoichiometric matrix of Pyrimidine Metabolism has rank " + str(pm_rank_numpy) + ".")
print("The stoichiometric matrix of Lipopolysaccharide Biosynthesis has rank " + str(lb_rank) + ".")
print("\nThe conservative relations in the Amino Sugar and Nucleotide Sugar Metabolism subsystem are: \n")
pp.pprint(find_conservative_relations(S_asansm))

print("\nThe conservative relations in the Pyrimidine Metabolism subsystem are: \n")
pp.pprint(find_conservative_relations(S_pm))

print("\nThe conservative relations in the Lipopolysaccharide Biosynthesis subsystem are: \n")
pp.pprint(find_conservative_relations(S_lb))

print("\nThe conservative relations in the combined subsystem are: \n")
pp.pprint(find_conservative_relations(S_combined))
S_asansm_latex = asansm_model.getFullStoichiometryMatrix()
S_asansm_latex_df = pd.DataFrame(data=S_asansm_latex, index=S_asansm_latex.rownames, columns=S_asansm_latex.colnames, dtype=int)
with open("matrices/S_asansm_latex.txt", "w") as outfile:
    outfile.write(S_asansm_latex_df.to_latex(bold_rows=True))

S_asansm_latex = asansm_model.getFullStoichiometryMatrix()
asansm_left_null_latex_df = pd.DataFrame(data=[[elem for elem in line] for line in asansm_left_nullspace],
                                         columns=S_asansm_latex.rownames, dtype=int)
with open("matrices/asansm_left_nullspace_latex.txt","w") as outfile:
    outfile.write(asansm_left_null_latex_df.to_latex(bold_rows=True))
# Todo next: Redo with _fbc'd models, repaint pathway in pm map (mark 2 UDPs)
