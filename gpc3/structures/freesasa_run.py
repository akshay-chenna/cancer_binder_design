import freesasa

structure = freesasa.Structure("correct_uniprot_gpc3_59-477.pdb")
result = freesasa.calc(structure)
residue_areas = result.residueAreas()

test = residue_areas['A'].items()

relative_total = []
for a, b in test:
	relative_total.append(b.relativeTotal)

relative_apolar = []
for a, b in test:
	relative_apolar.append(b.relativeApolar)

relative_polar = []
for a, b in test:
	relative_polar.append(b.relativePolar)

relative_total_str = [str(item) for item in relative_total]
with open('relative_total_sasa.txt','w') as file:
	file.write('\n'.join(relative_total_str))

relative_apolar_str = [str(item) for item in relative_apolar]
with open('relative_apolar_sasa.txt','w') as file:
	file.write('\n'.join(relative_apolar_str))

relative_polar_str = [str(item) for item in relative_polar]
with open('relative_polar_sasa.txt','w') as file:
	file.write('\n'.join(relative_polar_str))
