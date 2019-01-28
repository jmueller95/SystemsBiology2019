import libsbml
from collections import defaultdict

"""CODE FROM MEPHENOR START"""
doc = libsbml.readSBMLFromFile("iPAE1146_with_groups.xml")
model = doc.getModel()

groups_plugin = model.getPlugin("groups")
groups = groups_plugin.getListOfGroups()

subsystem_balance = dict()
neighbors = dict()
for group in groups:
    net_balance = defaultdict(int)
    reactant_balance = defaultdict(int)
    product_balance = defaultdict(int)
    for member in group.getListOfMembers():
        r_id = member.getIdRef()
        reaction = model.getReaction(r_id)
        for reactant in reaction.getListOfReactants():
            s_id = reactant.getSpecies()
            name = model.getSpecies(s_id).getName()
            stoichiometry = reactant.getStoichiometry()
            reactant_balance[name] += stoichiometry
        for product in reaction.getListOfProducts():
            s_id = product.getSpecies()
            name = model.getSpecies(s_id).getName()
            stoichiometry = product.getStoichiometry()
            product_balance[name] += stoichiometry
        keys = set(product_balance.keys()).union(reactant_balance.keys())
        for key in keys:
            net_balance[key] += product_balance[key] - reactant_balance[key]
        to_remove = list()
        for key in net_balance.keys():
            if net_balance[key] == 0:
                to_remove.append(key)
        for key in to_remove:
            net_balance.pop(key)
    subsystem_balance[group.name] = net_balance
"""CODE FROM MEPHENOR END"""


def add_Boundary_Reactions(model, balance_dict):
    for species_name, rate in balance_dict.items():
        species = [species for species in model.getListOfSpecies() if species.getName() == species_name][0]
        species.setBoundaryCondition(True)
        rate_rule = model.createRateRule()
        rate_rule.setVariable(species.getId())
        rate_rule.setFormula(str(rate))


asansm_sbml = "fbc/iPAE1146_Amino_sugar_and_nucleotide_sugar_metabolism_fbc_squeezed.xml"
pm_sbml = "fbc/iPAE1146_Pyrimidine_metabolism_fbc_squeezed.xml"
lb_sbml = "fbc/iPAE1146_Lipopolysaccharide_biosynthesis_fbc_squeezed.xml"
combined_sbml = "fbc/iPAE1146_Combined_Subsystems_fbc_squeezed.xml"

asansm_libsbml_doc = libsbml.readSBML(asansm_sbml)
asansm_libsbml = asansm_libsbml_doc.getModel()
add_Boundary_Reactions(asansm_libsbml,
                       subsystem_balance["Amino sugar and nucleotide sugar metabolism"])
libsbml.writeSBMLToFile(asansm_libsbml_doc,
                        "xml/iPAE1146_Amino_sugar_and_nucleotide_sugar_metabolism_fbc_squeezed_with_boundaries.xml")

pm_libsbml_doc = libsbml.readSBML(pm_sbml)
pm_libsbml = pm_libsbml_doc.getModel()
add_Boundary_Reactions(pm_libsbml,
                       subsystem_balance["Pyrimidine metabolism"])
libsbml.writeSBMLToFile(pm_libsbml_doc,
                        "xml/iPAE1146_Pyrimidine_metabolism_fbc_squeezed_with_boundaries.xml")

lb_libsbml_doc = libsbml.readSBML(lb_sbml)
lb_libsbml = lb_libsbml_doc.getModel()
add_Boundary_Reactions(lb_libsbml,
                       subsystem_balance["Lipopolysaccharide biosynthesis"])
libsbml.writeSBMLToFile(lb_libsbml_doc,
                        "xml/iPAE1146_Lipopolysaccharide_biosynthesis_fbc_squeezed_with_boundaries.xml")

combined_balance = {}
for balance_dict in [subsystem_balance["Amino sugar and nucleotide sugar metabolism"],
                     subsystem_balance["Pyrimidine metabolism"],
                     subsystem_balance["Lipopolysaccharide biosynthesis"]]:
    for key, value in balance_dict.items():
        if key in combined_balance:
            combined_balance[key] += value
        else:
            combined_balance[key] = value

combined_libsbml_doc = libsbml.readSBML(combined_sbml)
combined_libsbml = combined_libsbml_doc.getModel()
add_Boundary_Reactions(combined_libsbml,
                       combined_balance)
libsbml.writeSBMLToFile(combined_libsbml_doc,
                        "xml/iPAE1146_Combined_Subsystems_fbc_squeezed_with_boundaries.xml")
