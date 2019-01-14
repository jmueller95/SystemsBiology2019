import json
from pprint import pprint

"""
remove quotation marks for using the print functions
"""


# Open JSON file
jsonFile = "../iPAE1146.json"

with open(jsonFile) as f:
    data = json.load(f)

# Select data from JSON file
metabolites = data["metabolites"]
reactions = data["reactions"]
genes = data["genes"]

# Choose the name of the subsystem to be analyzed
subsystem_name = 'Lipopolysaccharide biosynthesis'

# Find all the reaction in the chosen subsystem
lipopolysaccharide_reactions = []
for reaction in reactions:
    if reaction["subsystem"] == subsystem_name:
        lipopolysaccharide_reactions.append(reaction)
"""
print("lipopolysaccharide_reactions")
pprint(lipopolysaccharide_reactions)
"""

# Find all metabolites ids in the chosen subsystem
lipopolysaccharide_metabolites_ids = list()
for reaction in lipopolysaccharide_reactions:
    for metabolite in reaction['metabolites'].keys():
        if metabolite not in lipopolysaccharide_metabolites_ids:
            lipopolysaccharide_metabolites_ids.append(metabolite)

# Find all information regarding the metabolites in the chosen subsystem
lipopolysaccharide_metabolites = []
for metabolite in metabolites:
    if metabolite['id'] in lipopolysaccharide_metabolites_ids:
        lipopolysaccharide_metabolites.append(metabolite)
"""
print("lipopolysaccharide_metabolites")
pprint(lipopolysaccharide_metabolites)
print('Number of lipopolysaccharide_metabolites: ', len(lipopolysaccharide_metabolites))
print('Number of lipopolysaccharide_reactions: ', len(lipopolysaccharide_reactions))
"""


# Associate metabolites with the reactions of the chosen subsystem in which they are present
metabolite_in_reaction = {}
for metabolite in lipopolysaccharide_metabolites_ids:
    for reaction in lipopolysaccharide_reactions:
        if metabolite in reaction['metabolites']:
            metabolite_in_reaction.setdefault(metabolite, []).append(reaction['id'])

print('metabolites in the same reactions: ')
pprint(metabolite_in_reaction)
