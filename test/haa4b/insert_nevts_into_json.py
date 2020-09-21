

input_json_file = open( "samples2018.json", "r" )

output_file     = open( "samples2018_with_nevts.json", "w" )

input_json_lines = input_json_file.readlines()


input_totals_file = open( "summary-2018-totals.txt", "r" )

input_totals_lines = input_totals_file.readlines()

totals_dict = dict()
for line in input_totals_lines:
   line = line.rstrip("\n")
   chunks = line.split()
   #print(" line ---{}---  key = {}, val = {}".format(line,chunks[0], chunks[1]) )
   totals_dict[str(chunks[0])] = chunks[1]

print("\n\n")
for key in totals_dict.keys():
   print( " key = {:70}  val = {:15}".format(key,totals_dict[key] ) )

print("\n\n")



for line in input_json_lines:
   output_file.write(line)
   if "dtag" in line:
      chunks = line.split(":")
      dtag = chunks[1]
      dtag = dtag.replace(' ','')
      dtag = dtag.replace('"','')
      dtag = dtag.replace(',','')
      dtag = dtag.replace('\n','')
      print("  found dtag ---{}---".format(dtag) )
      if "Data" in dtag:
         output_file.write('                    "nevts": 1,\n')
      else:
         if dtag not in totals_dict.keys():
            print("\n\n *** found dtag {} in json file but not in totals file.\n\n")
            exit()
         else:
            output_file.write('                    "nevts": {},\n'.format( totals_dict[dtag] ) )




