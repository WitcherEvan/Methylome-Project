import sys
import math

# ------------ Default values ------------ #
radius = 1000.0
startPlant = ''
maxPop = 150

if (len(sys.argv) > 1): # then we have parameters
    # Sort out the given parameters.
    index = 1
    while ( index <= len(sys.argv)-1 ):
        arg = sys.argv[index]

# ------------ Radius in km ------------ #
        if arg == 'kmRadius' and index+1 <= len(sys.argv)-1:
            try:
                radius = round(math.fabs(float(sys.argv[index+1])),2)
                index = index + 1
            except ValueError:
                print('Error: \'kmRadius\' must be followed by a number value.')
                sys.exit()

# ------------ Max population ------------ #
        elif arg == 'popSize' and index+1 <= len(sys.argv)-1:
            try:
                maxPop = round(math.fabs(float(sys.argv[index+1])))
                index = index + 1
            except ValueError:
                print('Error: \'popSize\' must be followed by a number value.')
                sys.exit()

# ------------ Identify the starting plant ------------ #
        elif arg == 'ecotypeID' and index+1 <= len(sys.argv)-1:
            eID = sys.argv[index+1]
            with open('masterList.csv', 'r') as master:
                for line in master:
                    if eID == line.split(',')[0]:
                        startPlant = line
                        break # from reading file
            index = index + 1

# ------------ Identify the starting plant ------------ #
        elif arg == 'ecotypeName' and index+1 <= len(sys.argv)-1:
            eName = sys.argv[index+1]
            with open('masterList.csv', 'r') as master:
                for line in master:
                    if eName == line.split(',')[1]:
                        startPlant = line
                        break # from reading file
            index = index + 1

# ------------ Identify the starting plant ------------ #
        elif arg == 'GSMnumber' and index+1 <= len(sys.argv)-1:
            eGSM = sys.argv[index+1]
            with open('masterList.csv', 'r') as master:
                for line in master:
                     # -1 loops around to get last element
                    if eGSM == line.split(',')[-1].split('_')[0]:
                        startPlant = line
                        break # from reading file
            index = index + 1

        index = index + 1 # end of iteration

if startPlant == '':
    print('A starting methylome is required.')
    sys.exit()

# ------------ Create dictionary ------------ #
web = {}
startLat = float(startPlant.split(',')[4])
startLon = float(startPlant.split(',')[5])
# convert coords to Radians
p1_lat = math.radians(startLat)
p1_lon = math.radians(startLon)

with open('masterList.csv', 'r') as master:
    next(master) # skip headers which cause typeError
    for line in master:
        lat = line.split(',')[4]
        if lat == '':
            continue # to next line, because coords are missing.
        lat = float(lat)
        lon = line.split(',')[5]
        if lon == '':
            continue # to next line, because coords are missing.
        lon = float(lon)

        # convert coords to Radians
        p2_lat = math.radians(lat)
        p2_lon = math.radians(lon)

        # haversine formula described here: https://www.movable-type.co.uk/scripts/latlong.html
        # a = sin²(Δφ/2) + cos φ1 ⋅ cos φ2 ⋅ sin²(Δλ/2)
        a = math.pow( math.sin(p1_lat - p2_lat), 2) + math.cos(p1_lat) * math.cos(p2_lat) * math.pow( math.sin(p1_lon - p2_lon) ,2)
        # c = atan2( √a, √(1−a) )     there was a *2 in there, but computed distances were 2*actual distances.
        c = math.atan2( math.sqrt(a), math.sqrt(1-a))
        distance = 6371 * c # mean radius of Earth = 6371 km

        if distance < radius: #
            GSMfile = line.split(',')[-1].strip('\n')
            web[GSMfile] = distance

# to get closest n plants, we are essentially sorting the dictionary. Room to improve this
# If there was a magic way to sort a dict based on values, that'd be great.

# two lists, so I can use indexes, and keep the values in parallel. Initialized with safe values
Distances = []
Distances.insert(0,22000) # circumference of Earth ~40'000km so no distance will be greater than this.
GSM_list = []
GSM_list.insert(0,'not_of_this_world')

# Here we order/sort the plants and their distances:
for GSMfile in web:
    d = web.get(GSMfile)
    idx = 0
    while idx <= len(Distances):
        if d < Distances[idx]:
            Distances.insert(idx,d)
            GSM_list.insert(idx,GSMfile)
            break # while loop, not web-loop
        else:
            idx = idx+1

Distances.pop() # to remove inital 'not_a_file' and large distance
GSM_list.pop()  # pop() removes last element, which a sufficiently large distance guarantees.

# Create a Batch file to list the methylomes and distances:
description = './Batch/Pop_max_'+str(maxPop)+'_radius_'+str(radius)+'_on_'+startPlant.split(',')[-1].split('_')[0]
with open(description, 'w') as result:
    idx = 0
    while idx < maxPop and idx < len(Distances):# Grab up to the specified maximum of closest methylomes:
        result.write( GSM_list[idx] +','+ str(round(Distances[idx],3)) +'\n')
        idx = idx+1
