import cantera as ct
    
gas = ct.Solution('MFuller.cti')
gas.TPX = 1690, 1000000, 'iC4H9O:0.0001,AR:0.9999'
 
# Define a reactor, let it react until the temperature reaches 1800 K:
r = ct.IdealGasReactor(gas)
sim = ct.ReactorNet([r])
sim.advance(2E-3)   # advance simulation to 2 ms
 
# Define the element to follow in the reaction path diagram:
element = 'C'
 
# Initiate the reaction path diagram:
diagram = ct.ReactionPathDiagram(gas, element)
 
# Options for cantera:
diagram.flow_type = 'NetFlow' # or OneWayFlow
diagram.threshold = 0.01
# diagram.scale = -1
diagram.font = 'Arial'
diagram.dot_options = 'dpi=100'
diagram.show_details = True
 
dot_text = diagram.get_dot() # write the dot-text 
'''
# Let's open the just created dot file and make so adjustments before generating the image
# The dot_file is opened and read, the adjustements are saved in the modified_dot_file:
for line in dot_text:
    # Remove the line with the label:
    if line.startswith(' label'):
        line = ""
    # Change lightness to zero, erase style command and replace setlinewidth with penwidth:
    line = line.replace(', 0.9"', ', 0.0"')
    line = line.replace('style="setlinewidth(', 'penwidth=')
    line = line.replace(')", arrowsize', ', arrowsize')
    # Find the lines with a color statement:
    if line.find('color="0.7, ') != -1:
        # Find the position of the saturation value:
        start=line.find('color="0.7, ')
        end=line.find(', 0.0"')
        saturation=float(line[start+12:end])
        if saturation > 1:
            # The highest values are set to 0 (black):
            saturationnew = 0
        else:
            # All values between 0 and 1 are inverted:
            saturationnew=round(abs(saturation-1),2)
    # Switch positions between saturation and lightness:
    line = line.replace(', 0.0"', '"')
    line = line.replace('0.7, ', '0.7, 0.0, ')
    # Replace saturation value with new calculated one:
    try:
        line = line.replace(str(saturation), str(saturationnew))
    except NameError:
        pass
'''
for species, molfrac in zip(gas.species_names, gas.X):
    print('{:15s}{:.2e}'.format(species, molfrac))

plot_dot(dot_text)