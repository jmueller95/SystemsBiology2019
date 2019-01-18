import tellurium as te
import libsbml


def get_reaction_names(listOfReactions, reactionIDs):
    return [listOfReactions.get(id).getName() for id in reactionIDs]


def get_metabolite_names(listOfMetabolites, metaboliteIDs):
    return [listOfMetabolites.get(id).getName() for id in metaboliteIDs]


# We start with the smallest subsystem: Amino Sugar and Nucleotide Sugar Metabolism
asansm_sbml = "xml/iPAE1146_Amino_sugar_and_nucleotide_sugar_metabolism.xml"
asansm_libsbml_doc = libsbml.readSBML(asansm_sbml)
asansm_libsbml = asansm_libsbml_doc.getModel()
asansm_model = te.loadSBMLModel(asansm_sbml)
S_asansm = asansm_model.getFullStoichiometryMatrix()
S_asansm.colnames = get_reaction_names(asansm_libsbml.getListOfReactions(), S_asansm.colnames)
S_asansm.rownames = get_metabolite_names(asansm_libsbml.getListOfSpecies(), S_asansm.rownames)

# Next subsystem: Pyrimidine Metabolism
pm_sbml = "xml/iPAE1146_Pyrimidine_metabolism.xml"
pm_libsbml_doc = libsbml.readSBML(pm_sbml)
pm_libsbml = pm_libsbml_doc.getModel()
pm_model = te.loadSBMLModel(pm_sbml)
S_pm = pm_model.getFullStoichiometryMatrix()
S_pm.colnames = get_reaction_names(pm_libsbml.getListOfReactions(), S_pm.colnames)
S_pm.rownames = get_metabolite_names(pm_libsbml.getListOfSpecies(), S_pm.rownames)

# Third one: Lipopolysaccharide Biosynthesis
lb_sbml = "xml/iPAE1146_Lipopolysaccharide_biosynthesis.xml"
lb_libsbml_doc = libsbml.readSBML(lb_sbml)
lb_libsbml = lb_libsbml_doc.getModel()
lb_model = te.loadSBMLModel(lb_sbml)
S_lb = lb_model.getFullStoichiometryMatrix()
S_lb.colnames = get_reaction_names(lb_libsbml.getListOfReactions(), S_lb.colnames)
S_lb.rownames = get_metabolite_names(lb_libsbml.getListOfSpecies(), S_lb.rownames)

# Todo next: Get Fluxes v (to formulate dx=Sv)
