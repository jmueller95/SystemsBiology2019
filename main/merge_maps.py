# This script merges the three Escher maps into one (by merging their Json files together)
import json

with open("../maps/iPAE1146_Amino_sugar_and_nucleotide_sugar_metabolism.json") as asansm_file:
    asansm_json = json.load(asansm_file)

with open("../maps/iPAE1146_Pyrimidine_metabolism.json") as pm_file:
    pm_json = json.load(pm_file)

with open("../maps/iPAE1146_Lipopolysaccharide_biosynthesis.json")as lb_file:
    lb_json = json.load(lb_file)

# The merged json will be a list of two dictionaries (we call them 'header' and 'content')

# We hard-code the header
merged_header = {"map_name": "merged_map",
                 "map_id": "sysbio2019",
                 "map_description": "Merging of three maps",
                 "homepage": "https://escher.github.io",
                 "schema": "https://escher.github.io/escher/jsonschema/1-0-0#"}


def rename_nodes_and_add_offset(subsystem_prefix, reactions, nodes, offset_x, offset_y):
    # Reactions
    for reaction in reactions.values():
        # Offset 'label_x' and 'label_y'
        reaction['label_x'] += offset_x
        reaction['label_y'] += offset_y
        for segment in reaction['segments'].values():
            # Prefix all "from/to_node_ids"
            segment['from_node_id'] = subsystem_prefix + segment['from_node_id']
            segment['to_node_id'] = subsystem_prefix + segment['to_node_id']
            # If 'b1/b2' exist, offset them
            if segment['b1'] is not None:
                segment['b1']['x'] += offset_x
                segment['b1']['y'] += offset_y
            if segment['b2'] is not None:
                segment['b2']['x'] += offset_x
                segment['b2']['y'] += offset_y
        # Prefix all keys of the segments
        reaction['segments'] = {subsystem_prefix + key: value for key, value in reaction['segments'].items()}

    # Nodes
    for node in nodes.values():
        # Offset labels and coordinates
        if 'label_x' in node:
            node['label_x'] += offset_x
        if 'label_y' in node:
            node['label_y'] += offset_y
        node['x'] += offset_x
        node['y'] += offset_y

    # Finally, rename the keys of the reactions and nodes themselves and return them
    return {subsystem_prefix + key: value for key, value in reactions.items()}, \
           {subsystem_prefix + key: value for key, value in nodes.items()}


# The content is a dict with four keys: reactions, nodes, canvas, text_labels
def merge_dicts(x, y, z):
    result = x.copy()
    result.update(y)
    result.update(z)
    return result


# Approach: We leave lipopolysaccharide metabolism where it is
# and move pyrimidine to the bottom and amino/nucleotide sugar to the right (and a little to the bottom)
lb_offset_x = 0
lb_offset_y = 0
pm_offset_x = 0
pm_offset_y = 1250
asansm_offset_x = 4500
asansm_offset_y = 750

lb_reactions, lb_nodes = rename_nodes_and_add_offset(
    "lb_", lb_json[1]['reactions'], lb_json[1]['nodes'], lb_offset_x, lb_offset_y)
pm_reactions, pm_nodes = rename_nodes_and_add_offset(
    "pm_", pm_json[1]['reactions'], pm_json[1]['nodes'], pm_offset_x, pm_offset_y)
asansm_reactions, asansm_nodes = rename_nodes_and_add_offset(
    "asansm_", asansm_json[1]['reactions'], asansm_json[1]['nodes'], asansm_offset_x, asansm_offset_y)

# Reactions need to be offset
reactions_merged = merge_dicts(asansm_reactions, pm_reactions, lb_reactions)
# Nodes need to be offset in the same way
nodes_merged = merge_dicts(asansm_nodes, pm_nodes, lb_nodes)
# Canvas needs to be extended based on maximum coordinates
canvas_merged = {'x': -1700, 'y': -1500, 'height': 5000, 'width': 8750}
# Text labels can simply be merged
text_labels_merged = merge_dicts(asansm_json[1]['text_labels'], pm_json[1]['text_labels'], lb_json[1]['text_labels'])
# We also add some custom text labels to label our submaps
text_labels_merged['lb_label'] = {"x": 1000, "y": -1250, "text": "Lipopolysaccharide Biosynthesis"}
text_labels_merged['pm_label'] = {"x": 1000.0, "y": 1400.0, "text": "Pyrimidine Metabolism"}
text_labels_merged['asansm_label'] = {"x": 4700, "y": 500,
                                      "text": "Amino Sugar And Nucleotide Sugar Metabolism"}

merged_content = {"reactions": reactions_merged,
                  "nodes": nodes_merged,
                  "canvas": canvas_merged,
                  "text_labels": text_labels_merged}

merged_json = [merged_header, merged_content]
merged_json_string = json.dumps(merged_json)
with open("../maps/combined_map.json", "w") as outfile:
    outfile.write(merged_json_string)
