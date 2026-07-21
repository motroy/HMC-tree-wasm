'''
Here we collect functions that edit either the json or the svg outputs of grapetree (https://achtman-lab.github.io/GrapeTree/MSTree_holder.html)
'''
import json
import os
import re
import matplotlib
import dan_auxiliary_functions as daf

def split_svg_file(svg_file):
    with open(svg_file) as file:
        lines = [line.rstrip() for line in file]
    svg_line = ' '.join(lines)
    by_line = svg_line.split("<")
    by_line = ["<"+i for i in by_line if len(i)>1]
    return by_line

def print_svg_lines(svg_lines,output_file):
    with open(output_file, "w") as f:
        for i in svg_lines:
            f.write(i+"\n")

def get_json_data(json_file_name):
    with open(json_file_name) as json_file:
        json_data = json.load(json_file)
    return json_data

def modify_new_samples_in_json(json_data,clone_map):
    for i in json_data['metadata'].keys():
        for j in clone_map.keys():
            if str(j) in str(i):
                json_data['metadata'][i]['new_samples'] = clone_map[j]
    json_data['initial_category'] = 'new_samples'
    json_data['layout_data']['display_category'] = 'new_samples'
    return json_data

def output_json_data(json_data,output_file):
    with open(output_file, "w") as outfile:
        json.dump(json_data, outfile)

def find_pie_charted_nodes(json_file):
    json_data = get_json_data(json_file)
    names_dictionary = {}
    for i in json_data['metadata'].keys():
        if json_data['metadata'][i]['ID'] != json_data['metadata'][i]['__Node']:
            shown = json_data['metadata'][i]['__Node']
            to_add = json_data['metadata'][i]['ID']
            if shown not in names_dictionary.keys():
                names_dictionary[shown] = []
            # Here we fix the '\n' in a middle of a node name
            if '\n' in to_add:
                to_add = to_add.replace('\n', '')
            names_dictionary[shown].append(to_add)
    return names_dictionary

# if index%2 is odd, we reduce the translate factor (text will be more up)
# if index%2 is even, we add to the translate factor (text will be more down)
# spacing factor is multiplied every two samples (0,1 --> 1 , 2,3 --> 2, etc) - for index%2==0, factor is 0.5x+1, for index%2==1, factor is 0.5x+0.5
def find_translate_factor(index,spacing_factor,start_point):
    # just to make things tidy
    sign = 0
    multiplying_factor = 0
    if index%2 == 0:
        sign = 1
        multiplying_factor = 0.5*index + 1
    if index%2==1:
        sign = -1
        multiplying_factor = 0.5*index + 0.5
    return spacing_factor*sign*multiplying_factor+start_point

# for each key in the multi_sample_nodes:
# Looks for the lines:
# <text class="node-group-number" dy=".71em" text-anchor="middle" font-size="12" font-family="sans-serif" transform="translate(0,-4)">|sampleName|
# </text>
# e.g.:
# <text class="node-group-number" dy=".71em" text-anchor="middle" font-size="12" font-family="sans-serif" transform="translate(0,-4)">J3472212
# </text>
# and adds after them the other samples, one up by spacing faactor, one down by spacing factor.
def add_multi_sample_nodes(svg_lines,multi_sample_nodes,spacing_factor):
    for i in multi_sample_nodes.keys():
        current_key_line = '<text class="node-group-number" dy=".71em" text-anchor="middle" font-size="12" font-family="sans-serif" transform="translate(0,-4)">'+i
        current_key_index = svg_lines.index(current_key_line)
        # If the total number of samples in the node is even (i.e. the number of samples in the dictionary value is odd), we'll rearrange the locations, splitting the total space using the spacing_factor.
        if len(multi_sample_nodes[i])%2==1:
            key_translate_factor = -4 - 0.5 * spacing_factor
            svg_lines[current_key_index] = '<text class="node-group-number" dy=".71em" text-anchor="middle" font-size="12" font-family="sans-serif" transform="translate(0,'+str(key_translate_factor)+')">'+i

        # If the total number of samples in the node is odd (i.e. the list in the dictionary value is even), the locations are symmetrical around the (0,-4) location of the original key
        for index,j in enumerate(multi_sample_nodes[i]):
            to_insert = current_key_index+2
            if len(multi_sample_nodes[i]) % 2 == 0:
                translate_factor = find_translate_factor(index,spacing_factor,-4)
            else:
                translate_factor = find_translate_factor(index, spacing_factor,-4 - 0.5 * spacing_factor)
            svg_lines.insert(to_insert,'<text class="node-group-number" dy=".71em" text-anchor="middle" font-size="12" font-family="sans-serif" transform="translate(0,'+str(translate_factor)+')">'+j)
            svg_lines.insert(to_insert+1,'</text>')
            to_insert+=2
    return svg_lines

def get_loci_string(cgmlst,x,y,font_size):
    with open(cgmlst) as f:
        first_line = f.readline().strip('\n')
    title = '#Loci: '+str(len(first_line.split())-1)
    returned = '<text x="'+str(x)+'" y="'+str(y)+'" font-size="'+str(font_size)+'">'+title+'</text>'
    return returned

# Looks for lines that start with <g class="node mst-element fixed" id="
# and forms a dictionary of node:line index
def get_mst_element_line_index(svg_lines):
    returned = {}
    for index,value in enumerate(svg_lines):
        if value.startswith('<g class="node mst-element fixed" id="'):
            returned[value.split('"')[3]] = index
    return returned

def get_clone_color_map(clone_color_map_file):
    returned = {}
    with open(clone_color_map_file) as file:
        lines = [line.rstrip() for line in file]
    for i in lines:
        after_split = i.split(',')
        for index,j in enumerate(after_split[1:]):
            returned[j] = (after_split[0],int(index))
    return returned

def add_pattern_to_svg(svg_lines,patterns_file):
    patterns_lines = daf.lines_from_file(patterns_file)
    svg_lines[1:1] = patterns_lines
    return svg_lines

# Here we change the "fill" property of the path segment, which is found in the line_index+1 for any particular node
def color_nodes_by_map(svg_lines,sample_clone_map,clone_color_map,node_mst_element_line_index):
    for i in node_mst_element_line_index.keys():
        current_clone = sample_clone_map[i]
        current_color = clone_color_map[current_clone][0]
        current_index = node_mst_element_line_index[i]+1
        hashtag_index = svg_lines[current_index].index('#')
        svg_lines[current_index] = svg_lines[current_index][:hashtag_index]+current_color+svg_lines[current_index][hashtag_index+7:]
    return svg_lines

# Here we:
# use the URL of the pattern files to change the fill property of the path segment, which is found in the line_index+1 for any particular node
# Example: <circle cx="8" cy="10" r="16" style="stroke-width: 0.5; stroke: black; fill: url(#NC);">
def color_nodes_by_pattern(svg_lines,sample_clone_map,clone_color_map,node_mst_element_line_index):
    for i in node_mst_element_line_index.keys():
        current_clone = sample_clone_map[i]
        current_color = 'fill="url(#'+clone_color_map[current_clone][0]+')" '
        current_index = node_mst_element_line_index[i] + 1
        fill_index = svg_lines[current_index].index('fill=')
        svg_lines[current_index] = svg_lines[current_index][:fill_index] + current_color + svg_lines[current_index][fill_index + 15:]
    return svg_lines

# Here we add second (third, etc) circles around nodes
def add_second_layer(svg_lines, sample_clone_map, clone_color_map, node_mst_element_line_index):
    for i in node_mst_element_line_index.keys():
        current_clone = sample_clone_map[i]
        if clone_color_map[current_clone][1]>0:
            current_index = node_mst_element_line_index[i]
            original_circle_size = float(svg_lines[current_index + 1].split(',')[1].split('A')[0])
            for j in range(clone_color_map[current_clone][1], 0, -1): # We might need more circles to accomodate for a third+ layer
                tmp_fill = svg_lines[current_index + 1].replace(str(original_circle_size),str(10 * (j+1) + original_circle_size))
                tmp_fill = tmp_fill[0:tmp_fill.index('fill')] + 'fill="none" style="stroke: black; stroke-width: 5px"/>'
                svg_lines.insert(current_index + 3, tmp_fill)
        node_mst_element_line_index = get_mst_element_line_index(svg_lines)
    return svg_lines

def get_batch_mark_map(sample_batch_map,batches_marks):
    returned = {}
    batches = list(set(sample_batch_map.values()))
    batches.sort()
    for i,v in enumerate(batches):
        returned[v] = batches_marks[i]
    return returned

# Here we add signs near nodes with new batch samples. Multiple batches are allowed
def mark_batches(svg_lines,sample_batch_map,multi_sample_nodes,node_mst_element_line_index,batch_mark_map):
    batches = list(set(sample_batch_map.values()))
    batches.sort()
    sample_batch_mark = {}
    for x in sample_batch_map.keys():
        sample_batch_mark[x] = batch_mark_map[sample_batch_map[x]]
    shown_node_marks = {x:[] for x in node_mst_element_line_index.keys()}
    # Now we check for each node, whether its sample is in a mapped batch, or if the other nodes in the node are in a mapped batch
    for i in node_mst_element_line_index.keys():
        if i in sample_batch_map.keys():
            shown_node_marks[i].append(sample_batch_mark[i])
        if i in multi_sample_nodes.keys():
            for j in multi_sample_nodes[i]:
                if j in sample_batch_map.keys():
                    shown_node_marks[i].append(sample_batch_mark[j])
    for i in shown_node_marks.keys():
        if len(shown_node_marks[i]) > 0:
            shown_node_marks[i] = list(set(shown_node_marks[i]))
            shown_node_marks[i].sort()
            marks_line = '<text dx="-1em" dy="-1em" font-size="65" font-weight="bold">'+''.join(shown_node_marks[i])+'</text>'
            # the markings are inserted right before the closing </g> of the node
            index_to_insert = svg_lines[node_mst_element_line_index[i]:].index('</g>')+node_mst_element_line_index[i]
            svg_lines.insert(index_to_insert,marks_line)
            node_mst_element_line_index = get_mst_element_line_index(svg_lines)
    return svg_lines

def add_batches(returned,batch_mark_map,locations,sample_batch_map):
    batches_marks = ['*', '%', '@', '$', '&', '|']
    if batch_mark_map!=None:
        sample_batch_dictionary = daf.dictionary_from_csv_file(sample_batch_map)
        batches = list(set(list(sample_batch_dictionary.values())))
        batches.sort()
        if len(locations) == 0:
            current_y = 0
        else:
            current_y = 70 * (len(locations)+1) - 20
        for i, v in enumerate(batches):
            returned.append('<g class="legend-item" transform="translate(0,' + str(current_y) + ')">')
            returned.append('<text x="0" y="9" dx="0.10em" dy=".45em" font-family="Arial" style="text-anchor: start;" font-size="20">' + v)
            returned.append('</text>')
            returned.append('<text x="150" y="9" dy=".35em" font-family="Arial" style="text-anchor: start;" font-size="20" font-weight="bold">' + batches_marks[i])
            returned.append('</text>')
            returned.append('</g>')
            current_y+=70
    return returned

def parse_legend_table(legend_table):
    lines = daf.lines_from_file(legend_table)
    header_row = lines[0].split("\t")
    header_column = lines[1].split("\t")
    locations = []
    for i in range(2,len(lines)):
        locations.append(lines[i].split("\t"))
    return(header_row,header_column,locations)

def get_default_locations(pattern_clone_legend_color_map):
    file_to_parse = pattern_clone_legend_color_map
    lines = daf.lines_from_file(file_to_parse)
    returned = []
    for i in lines:
        returned.append(i.split(',')[1:])
    return returned

def get_legend_header(items,row_or_column):
    returned = []
    for i,v in enumerate(items):
        if row_or_column == 'column':
            y_pos = 70*i-3
            returned.append('<text x="30" y="'+str(y_pos)+'" style="text-anchor: start;" font-size="20">'+v)
            returned.append("</text>")
        if row_or_column == 'row':
            x_pos = 125*(i+1)+30
            returned.append('<text x="'+str(x_pos)+'" y="0" style="text-anchor: start;" font-size="20">' + v)
            returned.append("</text>")
    return returned

def hex_to_rgb(hex):
    current_color = matplotlib.colors.to_rgb(hex)
    rgb_color = [round(x * 255) for x in current_color]
    returned = 'rgb(' + str(rgb_color[0]) + ', ' + str(rgb_color[1]) + ', ' + str(rgb_color[2]) + ')'
    return returned

def get_legend_shape_dictionary(pattern_clone_legend_color_map,clone_counter,legend_text):
    parsed_file = pattern_clone_legend_color_map
    returned = {}
    lines = daf.lines_from_file(parsed_file)
    for i in lines:
        after_split = i.split(',')
        fill = '="url(#'+after_split[0]+')"'
        for index,value in enumerate(after_split[1:]):
            returned[value] = ['<circle cx="8" cy="10" r="8" style="stroke-width: 1.5px; stroke: black;" fill'+fill+'>']
            returned[value].append('</circle>')
            # Adding concentric circles for extra layers, if necessary
            for j in range(index,0,-1):
                returned[value].append('<circle cx="8" cy="10" r="'+str(2.5*j+8)+'" style="stroke-width: 1.5px; stroke: black; fill:none;">')
                returned[value].append('</circle>')
            if legend_text==True:
                returned[value].append('<text x="20" y="9" dy=".35em" font-family="Arial" style="text-anchor: start;">'+value+'['+str(clone_counter[value])+']')
            else:
                returned[value].append('<text x="20" y="9" dy=".35em" font-family="Arial" style="text-anchor: start;">[' + str(clone_counter[value]) + ']')
            returned[value].append('</text>')
    return returned

def get_legend_shape_lines(legend_shape_dictionary,locations):
    returned = []
    for i,v in enumerate(locations):
        for i2,v2 in enumerate(v):
            if len(v2)>0:
                xpos = 125*i2+160
                ypos = 70*(i+1)-20
                returned.append('<g class="legend-item" transform="translate('+str(xpos)+','+str(ypos)+')">')
                returned.extend(legend_shape_dictionary[v2])
                returned.append('</g>')
    return returned

#def get_legend_lines_with_table(legend_table,pattern_clone_legend_color_map,clone_counter,batch_mark_map, legend_text, legend_x, legend_y):
def get_legend_lines_with_table(arguments,clone_counter,batch_mark_map):
    legend_table = arguments['legend_table']
    pattern_clone_legend_color_map = arguments['pattern_clone_legend_color_map']
    legend_text = arguments['legend_text']
    legend_x = arguments['legend_x']
    legend_y = arguments['legend_y']
    returned = ['<g class="legend" transform="translate('+str(legend_x)+','+str(legend_y)+')">']
    if legend_table=='NA' and pattern_clone_legend_color_map=='NA':
        return add_batches(returned,batch_mark_map,[])
    header_row = []
    header_column = []
    locations = get_default_locations(pattern_clone_legend_color_map) # Each clone on its own line. If there are second layers, they are in the line as the original layer.
    if legend_table!='NA':
        header_row,header_column,locations = parse_legend_table(legend_table)
    if len(header_row)>0:
        returned.extend(get_legend_header(header_row,'row'))
    if len(header_column)>0:
        returned.extend(get_legend_header(header_column,'column'))
    legend_shape_dictionary = get_legend_shape_dictionary(pattern_clone_legend_color_map,clone_counter,legend_text)
    returned.extend(get_legend_shape_lines(legend_shape_dictionary,locations))
    returned = add_batches(returned, batch_mark_map,locations,arguments['sample_batch_map'])
    returned.append('</g>')
    return returned

def remove_samples_text(svg_lines):
    returned = []
    for index,value in enumerate(svg_lines):
        sample_text = value.startswith('<text class="node-group-number"')
        end_sample_text = value=='</text>' and svg_lines[index-1].startswith('<text class="node-group-number"')
        if not (sample_text or end_sample_text):
            returned.append(value)
    return returned

# keeps one line if two consecutive lines contain </g>
def clean_svg_lines(svg_lines):
    indexes_to_remove = []
    for i in range(1,len(svg_lines)):
        if svg_lines[i]=="</g>" and svg_lines[i-1]=="</g>":
            indexes_to_remove.append(i)
    returned = [v for i,v in enumerate(svg_lines) if i not in indexes_to_remove]
    return returned