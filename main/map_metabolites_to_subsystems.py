from libsbml import *
from collections import Counter
from sbml_analysis_functions import reaction_ids_of_subsystem, reactionIds_to_subsystems, metabolites_dict

sbmlFile = "../iPAE1146_with_groups.xml"
# Code for later use (convert sbml to json)
# cobra_model = cobra.io.read_sbml_model(sbmlFile)
# cobra.io.save_json_model(cobra_model, "iPAE1146.json")

reader = SBMLReader()
document = reader.readSBML(sbmlFile)
libsbml_model = document.getModel()

# Get list of reactions associated with lipopolysaccharide synthesis
lipopolysaccharide_reactions = reaction_ids_of_subsystem('Lipopolysaccharide biosynthesis')

# Get set of metabolites that occur in lipopolysaccharide synthesis
lipopolysaccharide_metabolites = set()
for reaction in [libsbml_model.getReaction(id) for id in lipopolysaccharide_reactions]:
    reagent_list = reaction.getListOfReactants().clone()
    reagent_list.appendFrom(reaction.getListOfProducts())
    for reagent in reagent_list:
        if reagent not in lipopolysaccharide_metabolites:
            lipopolysaccharide_metabolites.add(reagent.getSpecies())
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

# Count occurrences of subsystems
counter = Counter([value for list in metabolites_to_subsystems.values() for value in list])
most_common_subsystems = counter.most_common(5)
reaction_counts = {subsystem + " (" + str(counter[subsystem]) + ")": len(reaction_ids_of_subsystem(subsystem)) for
                   subsystem in [tuple[0] for tuple in most_common_subsystems]}
"""Top three would be Lipopolysaccharide synthesis, pyrimidine metabolism and fatty acid biosynthesis,
however fatty acid biosynthesis has a lot of reactions so I suggest we take pyrimidine metabolism instead"""
