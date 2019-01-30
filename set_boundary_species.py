#!/usr/bin/env python
import libsbml
import operator

from collections import defaultdict, OrderedDict
from pprint import pprint

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
        
    subsystem_metabolites = dict()
    for key in subsystem_balance.keys():
        subsystem_metabolites[key] = list(subsystem_balance[key].keys())

    
    candidates = defaultdict(list)
    chosen_groups = ["Lipopolysaccharide biosynthesis",
                     "Amino sugar and nucleotide sugar metabolism", 
                     "Pyrimidine metabolism"]
    for this in chosen_groups:
        for other in subsystem_metabolites.keys():
            if this == other:
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
    set_boundary_species(combined_libsbml, list(boundary_species))
    libsbml.writeSBMLToFile(combined_libsbml_doc, "xml/iPAE1146_Combined_Subsystems_fbc_squeezed_with_boundaries.xml")


def set_boundary_species(model, candidates):
    species_occurrence = defaultdict(int)
    for reaction in model.getListOfReactions():
        for reactant in reaction.getListOfReactants():
            species = model.getSpecies(reactant.getSpecies())
            species_occurrence[species.getName()] += 1
        for product in reaction.getListOfProducts():
            species = model.getSpecies(product.getSpecies())
            species_occurrence[species.getName()] += 1
    
    to_remove = set()
    for key in species_occurrence.keys():
         if species_occurrence[key] == 1 or key in candidates:
                to_remove.add(key)
    
    for key in to_remove:
        species_occurrence.pop(key)
        
    non_boundary_species = [species for species in species_occurrence.keys()]
    for species in model.getListOfSpecies():
        if species.getName() not in non_boundary_species:
            species.setBoundaryCondition(True)
            
if __name__ == "__main__":
    main()
