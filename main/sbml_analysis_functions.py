# Most of this is obsolete now that we can use the groups package
from libsbml import *
import csv
import cobra

sbmlFile = "../iPAE1146_with_groups.xml"
reader = SBMLReader()
document = reader.readSBML(sbmlFile)
libsbml_model = document.getModel()
import xml.etree.ElementTree as ET

"""Map of reaction ids to subsystems"""
reactionIds_to_subsystems = {reaction.getId(): tag[tag.find(":") + 2:tag.find("<")] for reaction in
                             libsbml_model.getListOfReactions()
                             for tag in reaction.getNotesString().split("<p>")
                             if tag.startswith("SUBSYSTEM")}

"""Map of metabolite id to metabolite name, generated from the csv file"""
metabolites_dict = {}
with open('../iPae1146_metabolites.csv', 'rb') as csvfile:
    metabolites_reader = csv.reader(csvfile, delimiter=",", quotechar="\"")
    for row in metabolites_reader:
        id = "M_" + row[0][:-3] + "_c"
        name = row[1]
        metabolites_dict[id] = name


def reaction_ids_of_subsystem(subsystem_name):
    """This function returns a list of reaction ids for a given subsystem (given as a string)"""
    # Get list of reactions associated with the subsystem
    return [reaction for reaction in reactionIds_to_subsystems.keys() if
            reactionIds_to_subsystems[reaction] == subsystem_name]


def reactions_of_subsystem(subsystem_name):
    """This function returns a list of reactions for a given subsystem (given as a string)"""
    model_groups = libsbml_model.getPlugin("groups").getListOfGroups()
    subsystem_group = [group for group in model_groups if group.getName() == subsystem_name][0]
    return [member.getReferencedElement() for member in subsystem_group.getListOfMembers()]


def create_cobra_json(sbml_name):
    cobra.io.save_json_model(cobra.io.read_sbml_model(sbml_name), sbml_name.split(".")[0] + ".json")
