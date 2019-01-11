from libsbml import *
import csv

sbmlFile = "../iPAE1146.xml"
# Code for later use (convert sbml to json)
# cobra_model = cobra.io.read_sbml_model(path2sbml)
# cobra.io.save_json_model(cobra_model, "iPAE1146.json")

reader = SBMLReader()
document = reader.readSBML(sbmlFile)
libsbml_model = document.getModel()

# Map all the reactions in the model to their subsystems
reactionIds_to_subsystems = {reaction.getId(): tag[tag.find(":") + 2:tag.find("<")] for reaction in
                             libsbml_model.getListOfReactions()
                             for tag in reaction.getNotesString().split("<p>")
                             if tag.startswith("SUBSYSTEM")}

# Get list of reactions associated with lipopolysaccharide synthesis
lipopolysaccharide_reactions = [reaction for reaction in reactionIds_to_subsystems.keys()
                                if reactionIds_to_subsystems[reaction] == 'Lipopolysaccharide biosynthesis']

# Get set of metabolites that occur in lipopolysaccharide synthesis
lipopolysaccharide_metabolites = set()
for reaction in [libsbml_model.getReaction(id) for id in lipopolysaccharide_reactions]:
    reagent_list = reaction.getListOfReactants().clone()
    reagent_list.appendFrom(reaction.getListOfProducts())
    for reagent in reagent_list:
        if reagent not in lipopolysaccharide_metabolites:
            lipopolysaccharide_metabolites.add(reagent.getSpecies())

# Use csv file to generate metabolite id to metabolite name
metabolites_dict = {}
with open('../iPae1146_metabolites.csv', 'rb') as csvfile:
    metabolites_reader = csv.reader(csvfile, delimiter=",", quotechar="\"")
    for row in metabolites_reader:
        id = "M_" + row[0][:-3] + "_c"
        name = row[1]
        metabolites_dict[id] = name

# Map metabolites of lipopolysaccharide biosynthesis to all subsystems they occur in
metabolites_to_subsystems = {}
for metabolite_id in lipopolysaccharide_metabolites:
    for reaction in libsbml_model.getListOfReactions():
        reagent_list = reaction.getListOfReactants().clone()
        reagent_list.appendFrom(reaction.getListOfProducts())
        reagent_list_string = [reagent.getSpecies() for reagent in reagent_list]
        if metabolite_id in reagent_list_string:
            if metabolites_dict[metabolite_id] in metabolites_to_subsystems:
                metabolites_to_subsystems[metabolites_dict[metabolite_id]].add(
                    reactionIds_to_subsystems[reaction.getId()])
            else:
                metabolites_to_subsystems[metabolites_dict[metabolite_id]] = {
                reactionIds_to_subsystems[reaction.getId()]}

