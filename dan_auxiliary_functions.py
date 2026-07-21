'''
This file includes some dumb-but-useful functions that are very often used
'''

# Assumes it has a header line
def dictionary_from_csv_file(input_file, has_header=True):
    lines = lines_from_file(input_file)
    returned = {}
    start_index = 1
    if has_header == False:
        start_index = 0
    for i in lines[start_index:]:
        after_split = i.split(',')
        returned[after_split[0]] = after_split[1]
    return returned
def lines_from_file(file_name):
    with open(file_name) as file:
        lines = [line.rstrip() for line in file]
    return lines
def list_to_file(file_name,input_list):
    f = open(file_name, "w")
    for i in input_list:
        f.write(i + "\n")
    f.close()
    print("Created " + file_name + " with " + str(len(input_list)) + " records")