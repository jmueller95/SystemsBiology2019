#!/usr/bin/env python
import libsbml
import re

from collections import defaultdict

SPECIES_PATTERN = re.compile("M_\\w*_c")


def main():
    doc = libsbml.readSBMLFromFile("iPAE1146_with_groups.xml")
    model = doc.getModel()

    groups_plugin = model.getPlugin("groups")
    groups = groups_plugin.getListOfGroups()

    subsystem_balance = dict()
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

    chosen_groups = ["Amino sugar and nucleotide sugar metabolism",
                     "Lipopolysaccharide biosynthesis",
                     "Pyrimidine metabolism"]

    subsystem_balance["Combined"] = subsystem_balance[chosen_groups[0]]
    for group_name in chosen_groups[1:]:
        balance = subsystem_balance[group_name]
        for key, value in balance.items():
            subsystem_balance["Combined"][key] += value

    subsystem_metabolites = dict()
    for key in subsystem_balance.keys():
        subsystem_metabolites[key] = list(subsystem_balance[key].keys())

    subsystem_metabolites["Combined"] = list(set(
        subsystem_balance[chosen_groups[0]].keys()).
        union(subsystem_balance[chosen_groups[1]].keys()).
        union(subsystem_balance[chosen_groups[2]].keys()))

    candidates = defaultdict(list)

    for this in chosen_groups:
        for other in subsystem_metabolites.keys():
            if this == other or this == "Combined" and other in chosen_groups:
                continue
            metabolites = set(subsystem_metabolites[this]).intersection(subsystem_metabolites[other])
            if metabolites:
                for metabolite in metabolites:
                    is_part_sink = subsystem_balance[this][metabolite] > 0 and subsystem_balance[other][metabolite] < 0
                    is_part_demand = subsystem_balance[this][metabolite] < 0 and subsystem_balance[other][metabolite] > 0
                    if is_part_sink or is_part_demand:
                        candidates[this].append(metabolite)

    add_boundary_species(candidates)


def add_boundary_species(candidates):
    asansm_sbml = "fbc/iPAE1146_Amino_sugar_and_nucleotide_sugar_metabolism_fbc_squeezed.xml"
    pm_sbml = "fbc/iPAE1146_Pyrimidine_metabolism_fbc_squeezed.xml"
    lb_sbml = "fbc/iPAE1146_Lipopolysaccharide_biosynthesis_fbc_squeezed.xml"
    combined_sbml = "fbc/iPAE1146_Combined_Subsystems_fbc_squeezed.xml"

    asansm_libsbml_doc = libsbml.readSBML(asansm_sbml)
    asansm_libsbml = asansm_libsbml_doc.getModel()
    set_boundary_species(asansm_libsbml, candidates["Amino sugar and nucleotide sugar metabolism"])
    libsbml.writeSBMLToFile(asansm_libsbml_doc, "xml/iPAE1146_Amino_sugar_and_nucleotide_sugar_metabolism_fbc_squeezed_with_boundaries.xml")

    pm_libsbml_doc = libsbml.readSBML(pm_sbml)
    pm_libsbml = pm_libsbml_doc.getModel()
    set_boundary_species(pm_libsbml, candidates["Pyrimidine metabolism"])
    libsbml.writeSBMLToFile(pm_libsbml_doc, "xml/iPAE1146_Pyrimidine_metabolism_fbc_squeezed_with_boundaries.xml")

    lb_libsbml_doc = libsbml.readSBML(lb_sbml)
    lb_libsbml = lb_libsbml_doc.getModel()
    set_boundary_species(lb_libsbml, candidates["Lipopolysaccharide biosynthesis"])
    libsbml.writeSBMLToFile(lb_libsbml_doc, "xml/iPAE1146_Lipopolysaccharide_biosynthesis_fbc_squeezed_with_boundaries.xml")

    combined_libsbml_doc = libsbml.readSBML(combined_sbml)
    combined_libsbml = combined_libsbml_doc.getModel()
    boundary_species = set()
    for candidate_list in candidates.values():
        for species in candidate_list:
            boundary_species.add(species)
    set_boundary_species(combined_libsbml, candidates["Combined"])
    libsbml.writeSBMLToFile(combined_libsbml_doc, "xml/iPAE1146_Combined_Subsystems_fbc_squeezed_with_boundaries.xml")


def set_boundary_species(model, candidates):
    species_occurrence = defaultdict(int)
    species_reaction_map = defaultdict(list)
    for reaction in model.getListOfReactions():
        kinetic_law = reaction.getKineticLaw()
        for l_param in kinetic_law.getListOfLocalParameters():
            model.addParameter(l_param)
            param = model.getParameter(l_param.getId())
            param.setConstant(False)
        for reactant in reaction.getListOfReactants():
            species = model.getSpecies(reactant.getSpecies())
            species_occurrence[species.getName()] += 1
            species_reaction_map[species.getName()].append(reaction.getId())
        for product in reaction.getListOfProducts():
            species = model.getSpecies(product.getSpecies())
            species_occurrence[species.getName()] += 1
            species_reaction_map[species.getName()].append(reaction.getId())

    to_remove = set()
    for key in species_occurrence.keys():
        if species_occurrence[key] == 1 or key in candidates:
            to_remove.add(key)

    for key in to_remove:
        species_occurrence.pop(key)

    # build filling/draining fluxes for species possibly exchanged with other subsystems
    non_boundary_species = [species for species in species_occurrence.keys()]
    for species in model.getListOfSpecies():
        if species.getName() not in non_boundary_species:
            for r_id in species_reaction_map[species.getName()]:
                reaction = model.getReaction(r_id)
                if species.getId() in reaction.getListOfReactants() or reaction.getReversible():
                    # add flux for filling
                    r_flux = model.createReaction()
                    r_flux.setId(r_id + "_" + species.getId() + "_fill")
                    r_flux.addProduct(species)
                    r_flux.setSBOTerm(627)
                    # discontinued in level 3 version 2, but needed in level 3 version 1
                    r_flux.setFast(False)
                    r_flux.setReversible(False)
                    kinetic_law = r_flux.createKineticLaw()
                    l_param = kinetic_law.createLocalParameter()
                    l_param.setId("FILL_RATE")
                    l_param.setValue(1)
                    l_param.setConstant(False)
                    l_param.setUnits("mmol_per_gDW_per_hr")
                    kinetic_law.setFormula("FILL_RATE")
                if species.getId() in reaction.getListOfProducts() or reaction.getReversible():
                    # add flux for emptying
                    r_flux = model.createReaction()
                    r_flux.setId(r_id + "_" + species.getId() + "_drain")
                    r_flux.addReactant(species)
                    r_flux.setSBOTerm(627)
                    # discontinued in level 3 version 2, but needed in level 3 version 1
                    r_flux.setFast(False)
                    r_flux.setReversible(False)
                    kinetic_law = r_flux.createKineticLaw()
                    kl = reaction.getKineticLaw()
                    formula = kl.getFormula()
                    modifiers = set()
                    for modifier in re.findall(SPECIES_PATTERN, formula):
                        modifiers.add(modifier)
                    for modifier in modifiers:
                        if modifier != species.getId():
                            modifier_reference = r_flux.createModifier()
                            modifier_reference.setSpecies(modifier)
                    kinetic_law.setFormula(formula)

    # merge filling reactions, as they do not depend on a kinetic law
    for species in model.getListOfSpecies():
        species.setInitialAmount(0)
        if species.getName() not in non_boundary_species:
            reactions = [reaction.getId() for reaction in model.getListOfReactions()
                         if species.getId() in reaction.getId() and "fill" in reaction.getId()]
            if len(reactions) > 1:
                for r_id in reactions:
                    model.removeReaction(r_id)
                r_flux = model.createReaction()
                r_flux.setId("R_fill_" + species.getId())
                r_flux.addProduct(species)
                # discontinued in level 3 version 2, but needed in level 3 version 1
                r_flux.setFast(False)
                r_flux.setReversible(False)
                kinetic_law = r_flux.createKineticLaw()
                l_param = kinetic_law.createLocalParameter()
                l_param.setId("FILL_RATE")
                l_param.setValue(1)
                l_param.setConstant(False)
                l_param.setUnits("mmol_per_gDW_per_hr")
                kinetic_law.setFormula("FILL_RATE")


if __name__ == "__main__":
    main()
