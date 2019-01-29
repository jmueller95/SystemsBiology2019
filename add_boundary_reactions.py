import libsbml
from collections import defaultdict


def main():
    """CODE FROM MEPHENOR START"""
    doc = libsbml.readSBMLFromFile("iPAE1146_with_groups.xml")
    model = doc.getModel()

    groups_plugin = model.getPlugin("groups")
    chosen_groups = ["Lipopolysaccharide biosynthesis",
                     "Amino sugar and nucleotide sugar metabolism",
                     "Pyrimidine metabolism"]

    net_fluxes = defaultdict(list)
    for group in groups_plugin.getListOfGroups():
        if group.getName() not in chosen_groups:
            continue
        for member in group.getListOfMembers():
            r_id = member.getIdRef()
            reaction = model.getReaction(r_id)
            sum_reactants = 0
            sum_products = 0
            for reactant in reaction.getListOfReactants():
                sum_reactants += reactant.getStoichiometry()
            for product in reaction.getListOfProducts():
                sum_products += product.getStoichiometry()
            net = sum_products - sum_reactants
            if net != 0:
                net_fluxes[group.getName()].append(reaction.getId())
    """CODE FROM MEPHENOR END"""

    asansm_sbml = "fbc/iPAE1146_Amino_sugar_and_nucleotide_sugar_metabolism_fbc_squeezed.xml"
    pm_sbml = "fbc/iPAE1146_Pyrimidine_metabolism_fbc_squeezed.xml"
    lb_sbml = "fbc/iPAE1146_Lipopolysaccharide_biosynthesis_fbc_squeezed.xml"
    combined_sbml = "fbc/iPAE1146_Combined_Subsystems_fbc_squeezed.xml"

    asansm_libsbml_doc = libsbml.readSBML(asansm_sbml)
    asansm_libsbml = asansm_libsbml_doc.getModel()
    add_boundary_reactions(asansm_libsbml, net_fluxes["Amino sugar and nucleotide sugar metabolism"])
    libsbml.writeSBMLToFile(asansm_libsbml_doc,
                            "xml/iPAE1146_Amino_sugar_and_nucleotide_sugar_metabolism_fbc_squeezed_with_boundaries.xml")

    pm_libsbml_doc = libsbml.readSBML(pm_sbml)
    pm_libsbml = pm_libsbml_doc.getModel()
    add_boundary_reactions(pm_libsbml, net_fluxes["Pyrimidine metabolism"])
    libsbml.writeSBMLToFile(pm_libsbml_doc, "xml/iPAE1146_Pyrimidine_metabolism_fbc_squeezed_with_boundaries.xml")

    lb_libsbml_doc = libsbml.readSBML(lb_sbml)
    lb_libsbml = lb_libsbml_doc.getModel()
    add_boundary_reactions(lb_libsbml, net_fluxes["Lipopolysaccharide biosynthesis"])
    libsbml.writeSBMLToFile(lb_libsbml_doc,
                            "xml/iPAE1146_Lipopolysaccharide_biosynthesis_fbc_squeezed_with_boundaries.xml")

    combined_balance = defaultdict(int)
    for balance_dict in [net_fluxes["Amino sugar and nucleotide sugar metabolism"],
                         net_fluxes["Pyrimidine metabolism"],
                         net_fluxes["Lipopolysaccharide biosynthesis"]]:
        for key, value in balance_dict.items():
            combined_balance[key] += value

    combined_libsbml_doc = libsbml.readSBML(combined_sbml)
    combined_libsbml = combined_libsbml_doc.getModel()
    add_boundary_reactions(combined_libsbml, combined_balance)
    libsbml.writeSBMLToFile(combined_libsbml_doc, "xml/iPAE1146_Combined_Subsystems_fbc_squeezed_with_boundaries.xml")


def add_boundary_reactions(model, boundary_reactions):
    species_occurrence = defaultdict(int)
    for reaction in model.getListOfReactions():
        for reactant in reaction.getListOfReactants():
            species_occurrence[reactant.getId()] += 1
        for product in reaction.getListOfProducts():
            species_occurrence[product.getId()] += 1

    for reaction in boundary_reactions:
        for reactant in reaction.getListOfReactants():
            if not species_occurrence[reactant.getId()] > 1:
                species = model.getSpecies(reactant.getId())
                species.setBoundaryCondition(True)
        for product in reaction.getListOfProducts():
            if not species_occurrence[product.getId()] > 1:
                species = model.getSpecies(product.getId())
                species.setBoundaryCondition(True)


if __name__ == "__main__":
    main()
