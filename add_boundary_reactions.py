#!/usr/bin/env python
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
    add_boundary_metabolites(asansm_libsbml, net_fluxes["Amino sugar and nucleotide sugar metabolism"])
    libsbml.writeSBMLToFile(asansm_libsbml_doc,
                            "xml/iPAE1146_Amino_sugar_and_nucleotide_sugar_metabolism_fbc_squeezed_with_boundaries.xml")

    pm_libsbml_doc = libsbml.readSBML(pm_sbml)
    pm_libsbml = pm_libsbml_doc.getModel()
    add_boundary_metabolites(pm_libsbml, net_fluxes["Pyrimidine metabolism"])
    libsbml.writeSBMLToFile(pm_libsbml_doc, "xml/iPAE1146_Pyrimidine_metabolism_fbc_squeezed_with_boundaries.xml")

    lb_libsbml_doc = libsbml.readSBML(lb_sbml)
    lb_libsbml = lb_libsbml_doc.getModel()
    add_boundary_metabolites(lb_libsbml, net_fluxes["Lipopolysaccharide biosynthesis"])
    libsbml.writeSBMLToFile(lb_libsbml_doc,
                            "xml/iPAE1146_Lipopolysaccharide_biosynthesis_fbc_squeezed_with_boundaries.xml")

    combined_libsbml_doc = libsbml.readSBML(combined_sbml)
    combined_libsbml = combined_libsbml_doc.getModel()
    for boundary_reactions in [net_fluxes["Amino sugar and nucleotide sugar metabolism"],
                               net_fluxes["Pyrimidine metabolism"],
                               net_fluxes["Lipopolysaccharide biosynthesis"]]:
        add_boundary_metabolites(combined_libsbml, boundary_reactions)
    libsbml.writeSBMLToFile(combined_libsbml_doc, "xml/iPAE1146_Combined_Subsystems_fbc_squeezed_with_boundaries.xml")


def add_boundary_metabolites(model, boundary_reactions):
    reactants = defaultdict(int)
    products = defaultdict(int)
    for reaction in model.getListOfReactions():
        for reactant in reaction.getListOfReactants():
            reactants[reactant.getSpecies()] += reactant.getStoichiometry()
        for product in reaction.getListOfProducts():
            products[product.getSpecies()] += product.getStoichiometry()

    for r_id in boundary_reactions:
        reaction = model.getReaction(r_id)
        for reactant in reaction.getListOfReactants():
            if reactants[reactant.getSpecies()] != products[reactant.getSpecies()]:
                species = model.getSpecies(reactant.getSpecies())
                species.setBoundaryCondition(True)
        for product in reaction.getListOfProducts():
            if reactants[product.getSpecies()] != products[product.getSpecies()]:
                species = model.getSpecies(product.getSpecies())
                species.setBoundaryCondition(True)


if __name__ == "__main__":
    main()
