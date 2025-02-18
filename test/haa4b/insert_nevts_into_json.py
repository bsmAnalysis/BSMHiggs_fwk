import json

# Set the year 
year = "2017"  

# Read the input JSON file as text
input_json_file = open("samples{}.json".format(year), "r")
input_json_lines = input_json_file.readlines()
input_json_file.close()

# Read the summary totals file
input_totals_file = open("summary-{}-totals.txt".format(year), "r")
input_totals_lines = input_totals_file.readlines()
input_totals_file.close()

# Create a dictionary from the totals file
totals_dict = dict()
for line in input_totals_lines:
    line = line.rstrip("\n")
    chunks = line.split()
    totals_dict[str(chunks[0])] = chunks[1]  # Store as string

print("\n\n")
for key in totals_dict.keys():
    print(" key = {:70}  val = {:15}".format(key, totals_dict[key]))
print("\n\n")

# Process JSON lines and modify the "nevts" value
updated_json_lines = []
inside_entry = False
dtag_value = None

for line in input_json_lines:
    stripped_line = line.strip()

    # Detect "dtag" and store the value
    if stripped_line.startswith('"dtag"'):
        chunks = stripped_line.split(":")
        dtag_value = chunks[1].strip().replace('"', '').replace(',', '')

        print("  found dtag ---{}---".format(dtag_value))

    # Detect "nevts" and replace it
    elif stripped_line.startswith('"nevts"') and dtag_value:
        if "Data" in dtag_value:
            new_nevts_line = '                    "nevts": 1,\n'
        elif dtag_value in totals_dict:
            new_nevts_line = '                    "nevts": {},\n'.format(totals_dict[dtag_value])
        else:
            print("*** Warning: dtag {} found in JSON but not in totals file. Skipping update...".format(dtag_value))
            new_nevts_line = line  # Keep the original line if dtag is not found

        # Replace the existing line
        updated_json_lines.append(new_nevts_line)
        continue  # Skip adding the original "nevts" line

    # Add the current line to output
    updated_json_lines.append(line)

# Write the updated JSON content to a new file
output_file = open("samples{}_with_nevts.json".format(year), "w")
output_file.writelines(updated_json_lines)
output_file.close()

print("Updated JSON file has been saved as 'samples{}_with_nevts.json'.".format(year))
