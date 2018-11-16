nathan_tax_file = '/home/dkoslicki/Data/MiCOPMinHash/Test/fnames_to_info.tsv'
output_file = '/home/dkoslicki/Data/MiCOPMinHash/Test/fnames_to_info_repophlan_format.tsv'
file_names = []
org_names_with_taxid = []
tax_ids = []
tax_path_names_with_taxid = []

prefixes = ['k', 'p', 'c', 'o', 'f', 'g', 's', 't']
with open(nathan_tax_file, 'r') as fid:
	for line in fid.readlines():
		line = line.strip()
		file_name, tax_path_names, tax_path_ids = line.split('\t')
		file_names.append(file_name)
		tax_path_names_split = tax_path_names.split('|')
		tax_path_ids_split = tax_path_ids.split('|')
		for index in reversed(range(len(tax_path_ids_split))):
			if tax_path_ids_split[index] and tax_path_names_split[index]:
				bottom_tax_id = tax_path_ids_split[index]
				bottom_tax_name = '_'.join(tax_path_names_split[index].split(' '))
				break
		org_name_with_taxid = bottom_tax_id + '_' + bottom_tax_name
		org_names_with_taxid.append(org_name_with_taxid)
		tax_ids.append(bottom_tax_id)
		tax_path_name_with_taxid = ''
		for index in range(len(tax_path_ids_split)):
			tax_id = tax_path_ids_split[index]
			name = '_'.join(tax_path_names_split[index].split(' '))
			prefix = prefixes[index]
			if tax_id and name:
				tax_path_name_with_taxid += prefix + '__' + tax_id + '_' + name
			tax_path_name_with_taxid += '|'
		tax_path_name_with_taxid = tax_path_name_with_taxid[:-2]  # remove the last pipe and strain
		tax_path_names_with_taxid.append(tax_path_name_with_taxid)

with open(output_file, 'w') as fid:
	for index in range(len(file_names)):
		file_name = file_names[index]
		org_name_with_taxid = org_names_with_taxid[index]
		tax_id = tax_ids[index]
		tax_path_name_with_taxid = tax_path_names_with_taxid[index]
		fid.write('\t'.join([file_name, org_name_with_taxid, tax_id, tax_path_name_with_taxid]))
		fid.write('\n')











