from libsbml import *
import csv
sbmlFile = "../iPAE1146_with_groups.xml"
reader = SBMLReader()
document = reader.readSBML(sbmlFile)
libsbml_model = document.getModel()

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

"""This function returns a list of reactions for a given subsystem (given as a string)"""
def reactions_of_subsystem(subsystem_name):

    # Get list of reactions associated with the subsystem
    return [reaction for reaction in reactionIds_to_subsystems.keys() if
            reactionIds_to_subsystems[reaction] == subsystem_name]
