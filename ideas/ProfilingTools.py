# This is a collection of scripts that will allow manipulation of CAMI profiling files
import sys
import copy
import os
import numpy as np
import timeit


class Profile(object):
	def __init__(self, input_file_name=None):
		# Initialize file name (if appropriate)
		self.input_file_name = input_file_name
		self._data = dict()
		# Stick in the root node just to make sure everything is consistent
		self._data["-1"] = dict()
		self._data["-1"]["rank"] = None
		self._data["-1"]["tax_path"] = list()
		self._data["-1"]["tax_path_sn"] = list()
		self._data["-1"]["abundance"] = 0
		self._data["-1"]["descendants"] = list()
		self._data["-1"]["branch_length"] = 0
		self._header = list()
		self._tax_id_pos = None
		self._rank_pos = None
		self._tax_path_pos = None
		self._tax_path_sn_pos = None
		self._abundance_pos = None
		self._eps = .0000000000000001  # This is to act like zero, ignore any lines with abundance below this quantity
		self._all_keys = ["-1"]
		self._merged_flag = False

		if self.input_file_name:
			if not os.path.exists(self.input_file_name):
				print("Input file %s does not exist" % self.input_file_name)
				raise Exception
			else:
				self.parse_file()

	def parse_file(self):
		input_file_name = self.input_file_name
		_data = self._data
		_header = self._header
		_all_keys = self._all_keys
		# all_tax_ids = set()
		with open(input_file_name, 'r') as read_handler:
			for line in read_handler:
				line = line.rstrip()
				if len(line) == 0:
					continue  # skip blank lines
				if line[0] == '@' and line[1] == '@':
					headers = line.strip().split()
					for header_iter in range(len(headers)):
						header = headers[header_iter]
						header = header.replace('@', '')
						if header == 'TAXID':
							tax_id_pos = header_iter
							self._tax_id_pos = tax_id_pos
						elif header == 'RANK':
							rank_pos = header_iter
							self._rank_pos = rank_pos
						elif header == 'TAXPATH':
							tax_path_pos = header_iter
							self._tax_path_pos = tax_path_pos
						elif header == 'TAXPATHSN' or header == "TAXPATH_SN":
							tax_path_sn_pos = header_iter
							self._tax_path_sn_pos = tax_path_sn_pos
						elif header == 'PERCENTAGE':
							abundance_pos = header_iter
							self._abundance_pos = abundance_pos
				if line[0] in ['@', '#']:
					_header.append(line)  # store data and move on
					continue
				if not all([isinstance(x, int) for x in [tax_id_pos, tax_path_pos, abundance_pos]]):
					print("Appears the headers TAXID, TAXPATH, and PERCENTAGE are missing from the "
											"header (should start with line @@)")
					sys.exit(2)
				temp_split = line.split('\t')
				tax_id = temp_split[tax_id_pos].strip()
				_all_keys.append(tax_id)
				tax_path = temp_split[tax_path_pos].strip().split("|")  # this will be a list, join up late
				# all_tax_ids.update(tax_path)
				abundance = float(temp_split[abundance_pos].strip())
				if isinstance(rank_pos, int):  # might not be present
					rank = temp_split[rank_pos].strip()
				if isinstance(tax_path_sn_pos, int):  # might not be present
					tax_path_sn = temp_split[tax_path_sn_pos].strip().split("|")  # this will be a list, join up later
				if tax_id in _data:  # If this tax_id is already present, add the abundance. NOT CHECKING FOR CONSISTENCY WITH PATH
					_data[tax_id]["abundance"] += abundance
					_data[tax_id]["tax_path"] = tax_path
					if isinstance(rank_pos, int):  # might not be present
						_data[tax_id]["rank"] = rank
					if isinstance(tax_path_sn_pos, int):  # might not be present
						_data[tax_id]["tax_path_sn"] = tax_path_sn
					# Find the ancestor
					if len(tax_path) <= 1:
						_data[tax_id]["ancestor"] = "-1"  # no ancestor, it's a root
						_data[tax_id]["branch_length"] = 1
						ancestor = "-1"
					else:
						ancestor = tax_path[-2]
						_data[tax_id]["branch_length"] = 1
						i = -3
						while ancestor is "" or ancestor == tax_id:  # if it's a blank or repeated, go up until finding ancestor
							ancestor = tax_path[i]
							_data[tax_id]["branch_length"] += 1
							i -= 1
						_data[tax_id]["ancestor"] = ancestor
				else:  # Otherwise populate the data
					_data[tax_id] = dict()
					_data[tax_id]["abundance"] = abundance
					_data[tax_id]["tax_path"] = tax_path
					if isinstance(rank_pos, int):  # might not be present
						_data[tax_id]["rank"] = rank
					if isinstance(tax_path_sn_pos, int):  # might not be present
						_data[tax_id]["tax_path_sn"] = tax_path_sn
					# Find the ancestor
					if len(tax_path) <= 1:
						_data[tax_id]["ancestor"] = "-1"  # no ancestor, it's a root
						_data[tax_id]["branch_length"] = 1
						ancestor = "-1"
					else:
						ancestor = tax_path[-2]
						_data[tax_id]["branch_length"] = 1
						i = -3
						while ancestor is "" or ancestor == tax_id:  # if it's a blank or repeated, go up until finding ancestor
							ancestor = tax_path[i]
							_data[tax_id]["branch_length"] += 1
							i -= 1
							if i+len(tax_path)<0:  # no ancestor available, manually set to root (-1) 
								ancestor="-1"
								break
						_data[tax_id]["ancestor"] = ancestor
				# Create a placeholder descendant key initialized to [], just so each tax_id has a descendant key associated to it
				if "descendants" not in _data[tax_id]:  # if this tax_id doesn't have a descendant list,
					_data[tax_id]["descendants"] = list()  # initialize to empty list
				# add the descendants
				if ancestor in _data:  # see if the ancestor is in the data so we can add this entry as a descendant
					if "descendants" not in _data[ancestor]:  # if it's not present, create the descendant list
						_data[ancestor]["descendants"] = list()
					_data[ancestor]["descendants"].append(tax_id)  # since ancestor is an ancestor, add this descendant to it
				else:  # if it's not already in the data, create the entry
					_data[ancestor] = dict()
					_data[ancestor]["descendants"] = list()
					_data[ancestor]["descendants"].append(tax_id)
		# only fix if need be
		# if all_tax_ids.intersection(_all_keys):
		self._delete_missing()  # make sure there aren't any missing internal nodes

	def _delete_missing(self):
		# This is really only useful for MetaPhlAn profiles, which due to the infrequently updated taxonomy, contain
		# many missing internal nodes, just delete the tax_id (make them "") and adjust branch lengths accordingly
		_data = self._data
		_good_keys = self._all_keys
		current_keys = _data.keys()
		bad_keys = []
		# Remove the bad keys from the dictionary
		for key in current_keys:  # go through the current keys
			if key not in _good_keys:  # if it's not a good key
				bad_keys.append(key)  # store the bad key
				del _data[key]  # delete the entry
		# Get the bad keys in the intermediate places
		for key in _good_keys:
			tax_path = _data[key]["tax_path"]
			for tax_id in tax_path:  # look in full taxpath
				if tax_id not in _good_keys:  # if it's a bad key
					if tax_id not in bad_keys:  # and it's not already in the list
						bad_keys.append(tax_id)  # add it to the list
		# Don't regard the blank tax_id
		if '' in bad_keys:
			bad_keys.pop(bad_keys.index(''))
		for key in _good_keys:
			# branch_length = _data[key]["branch_length"]
			tax_path = _data[key]["tax_path"]
			# tax_path_sn = _data[key]["tax_path_sn"]
			descendants = _data[key]["descendants"]
			# ancestor = _data[key]["ancestor"]
			for descendant in descendants:
				if descendant in bad_keys:
					descendants.pop(descendants.index(descendant))  # get rid of the bad descendants
			bad_indicies = []  # find the indicies of the bad ones in the tax path
			for tax_id in tax_path:
				if tax_id in bad_keys:
					index = tax_path.index(tax_id)
					bad_indicies.append(index)  # store all the bad indicies
			bad_indicies.sort(reverse=True)  # in reverse order
			# branch_length += len(bad_indicies)  # increment the branch_lengths accordingly (might not work for good|bad|good|bad
			for index in bad_indicies:
				tax_path[index] = ''  # remove the bad tax_ids
			# fix the branch lengths and find the ancestors
			if len(tax_path) >= 2:
				ancestor = tax_path[-2]
				_data[key]["branch_length"] = 1
				i = -3
				while ancestor is "" or ancestor == key:  # if it's a blank or repeated, go up until finding ancestor
					if i < -len(tax_path):  # Path is all the way full with bad tax_ids, connect to root
						_data[key]["branch_length"] += 1
						ancestor = "-1"
						_data["-1"]["descendants"].append(key)  # note this is now a descendant of the root
						break
					else:
						ancestor = tax_path[i]
						_data[key]["branch_length"] += 1
						i -= 1
				_data[key]["ancestor"] = ancestor
				if ancestor in _data:
					if key not in _data[ancestor]["descendants"]:
						_data[ancestor]["descendants"].append(key)  # Note this is now a descendant of it's ancestor
		return

	def _populate_missing_dont_use(self):
		# Unfortunately, some of the profile files may be missing intermediate ranks,
		# so we should manually populate them here
		# This will only really fix one missing intermediate rank....
		# INSTEAD!!! Let's just delete the missing intermediate ranks, if the tax_id isn't a key, delete it from the tax_id_path
		# This is really only useful for MetaPhlAn profiles, which due to the infrequently updated taxonomy, contain
		# many missing internal nodes
		_data = self._data
		for key in _data.keys():
			if "abundance" not in _data[key]:  # this is a missing intermediate rank
				# all the descendants *should* be in there, so leverage this info
				if "descendants" not in _data[key]:
					print("You're screwed, malformed profile file with rank %s" % key)
					raise Exception
				else:
					descendant = _data[key]["descendants"][0]
					to_populate_key = descendant  # just need the first one since the higher up path will be the same
					to_populate = copy.deepcopy(_data[to_populate_key])
					tax_path = to_populate["tax_path"]
					tax_path_sn = to_populate["tax_path_sn"]
					descendant_pos = tax_path.index(descendant)
					for i in range(len(tax_path) - 1, descendant_pos - 1, -1):
						tax_path.pop(i)
						tax_path_sn.pop(i)
					to_populate["branch_length"] = 1
					if "rank" in to_populate:
						rank = to_populate["rank"]
						if rank == "strain":
							rank = "species"
						elif rank == "species":
							rank = "genus"
						elif rank == "genus":
							rank = "family"
						elif rank == "family":
							rank = "order"
						elif rank == "order":
							rank = "class"
						elif rank == "class":
							rank = "phylum"
						elif rank == "phylum":
							rank = "superkingdom"
						else:
							print('Invalid rank')
							raise Exception
					else:
						print('Missing rank')
						raise Exception
					to_populate["ancestor"] = tax_path[-2]
					to_populate["rank"] = rank
					# Now go through and sum up the abundance for which this guy is the ancestor
					to_populate["abundance"] = 0
					to_populate["descendants"] = []
					for temp_key in _data.keys():
						if "ancestor" in _data[temp_key]:
							if _data[temp_key]["ancestor"] == key:
								to_populate["abundance"] += _data[temp_key]["abundance"]
								to_populate["descendants"].append(temp_key)
					# Make sure this guy is listed as a descendant to his ancestor
					if key not in _data[to_populate["ancestor"]]["descendants"]:
						_data[to_populate["ancestor"]]["descendants"].append(key)
					_data[key] = to_populate
		return

	def write_file(self, out_file_name=None):
		if out_file_name is None:
			raise Exception
		_data = self._data
		keys = _data.keys()
		# This will be annoying to keep things in order...
		# Let's iterate on the length of the tax_path since we know that will be in there
		tax_path_lengths = max([len(_data[key]["tax_path"]) for key in keys])
		fid = open(out_file_name, 'w')
		# Write the header
		for head in self._header:
			fid.write("%s\n" % head)

		# Loop over length of tax_path and write data
		# always make the output tax_id, rank, tax_path, tax_path_sn, abundance in that order
		for path_length in xrange(1, tax_path_lengths + 1):
			for key in keys:
				if len(_data[key]["tax_path"]) == path_length and _data[key]["abundance"] > self._eps:
					line_data = _data[key]
					fid.write("%s\t" % key)
					if self._rank_pos is not None:
						fid.write("%s\t" % line_data["rank"])
					fid.write("%s\t" % "|".join(line_data["tax_path"]))
					if self._tax_path_sn_pos is not None:
						fid.write("%s\t" % "|".join(line_data["tax_path_sn"]))
					fid.write("%f\n" % line_data["abundance"])
		fid.close()
		return

	def threshold(self, threshold=None):
		if threshold is None:
			raise Exception
		_data = self._data
		keys = _data.keys()
		for key in keys:
			if _data[key]["abundance"] < threshold:
				_data[key]["abundance"] = 0
		return

	def _subtract_down(self):
		# helper function to push all the weights up by subtracting
		# NOTE: when subtracting, need to start at root and go down
		# NOTE: when adding, need to start at leaves and go up
		_data = self._data
		keys = _data.keys()
		# This will be annoying to keep things in order...
		# Let's iterate on the length of the tax_path since we know that will be in there
		tax_path_lengths = max([len(_data[key]["tax_path"]) for key in keys])
		for path_length in range(1, tax_path_lengths):  # eg tax_path_lengths = 5, use 1,2,3,4 since we stop at leaves
			for key in keys:
				if len(_data[key]["tax_path"]) == path_length:
					descendants = _data[key]["descendants"]  # get all descendants
					for descendant in descendants:
						_data[key]["abundance"] -= _data[descendant]["abundance"]  # subtract the descendants abundance

	def _add_up(self):
		# helper function to push all the weights up by subtracting
		# NOTE: when subtracting, need to start at root and go down
		# NOTE: when adding, need to start at leaves and go up
		_data = self._data
		keys = _data.keys()
		keys = list(set(keys) - {'-1'})
		# This will be annoying to keep things in order...
		# Let's iterate on the length of the tax_path since we know that will be in there
		tax_path_lengths = max([len(_data[key]["tax_path"]) for key in keys])
		for path_length in range(tax_path_lengths, 1, -1):  # eg tax_path_lengths = 5, use 5,4,3,2, since we stop at roots
			for key in keys:
				if len(_data[key]["tax_path"]) == path_length:
					ancestor = _data[key]["ancestor"]
					if ancestor in _data:  # don't do anything if this is a/the root node
						_data[ancestor]["abundance"] += _data[key]["abundance"]  # add the descendants abundance

	def normalize(self):
		# Need to really push it up while subtracting, then normalize, then push up wile adding
		# self._push_up(operation="subtract")
		self._subtract_down()
		_data = self._data
		keys = _data.keys()
		total_abundance = 0
		for key in keys:
			total_abundance += _data[key]["abundance"]
		# print(total_abundance)
		for key in keys:
			if total_abundance > 0:
				_data[key]["abundance"] /= total_abundance
				_data[key]["abundance"] *= 100  # make back into a percentage
		# self._push_up(operation="add")
		self._add_up()
		return

	def merge(self, other):
		# Warning: not checking for taxonomic consistency
		if not isinstance(other, Profile):
			print("Only works with other Profiles")
			raise Exception
		if self._merged_flag is False:
			self._header.insert(0, "# This is a merged file, ignore files in headers below")
			self._merged_flag = True
		_data = self._data
		_other_data = other._data
		other_keys = _other_data.keys()
		for key in other_keys:
			if key in _data:
				_data[key]["abundance"] += _other_data[key]["abundance"]  # if already in there, add abundances
			else:
				_data[key] = copy.copy(_other_data[key])  # otherwise use the whole thing

	def make_unifrac_input_and_normalize(self, other):
		if not isinstance(other, Profile):
			raise Exception
		_data = self._data
		_other_data = other._data

		_data_keys = _data.keys()
		tax_path_lengths1 = max([len(_data[key]["tax_path"]) for key in _data_keys])
		_other_data_keys = _other_data.keys()
		tax_path_lengths2 = max([len(_other_data[key]["tax_path"]) for key in _other_data_keys])
		tax_path_lengths = max(tax_path_lengths1, tax_path_lengths2)
		all_keys = set(_data_keys)
		all_keys.update(_other_data_keys)
		nodes_in_order = []
		for path_length in range(tax_path_lengths, 0, -1):
			for key in all_keys:
				if key in _data:
					if len(_data[key]["tax_path"]) == path_length:
						if key not in nodes_in_order:
							nodes_in_order.append(key)
				elif key in _other_data:
					if len(_other_data[key]["tax_path"]) == path_length:
						if key not in nodes_in_order:
							nodes_in_order.append(key)
		# Make the graph
		# Put the root at the very end
		if '-1' in nodes_in_order:
			nodes_in_order.pop(nodes_in_order.index('-1'))
			nodes_in_order.append('-1')
		else:
			nodes_in_order.append('-1')
		Tint = dict()
		lint = dict()
		for key in nodes_in_order:
			if key in _data:
				if "ancestor" in _data[key]:  # If ancestor is not in there, then it's an ancestor
					ancestor = _data[key]["ancestor"]
					Tint[key] = ancestor
					lint[key, ancestor] = _data[key]["branch_length"]
			elif key in _other_data:
				if "ancestor" in _other_data[key]:
					ancestor = _other_data[key]["ancestor"]
					Tint[key] = ancestor
					lint[key, ancestor] = _other_data[key]["branch_length"]
		nodes_to_index = dict(zip(nodes_in_order, range(len(nodes_in_order))))
		#return Tint, lint, nodes_in_order, nodes_to_index

		# Now need to change over to the integer-based indexing
		Tint2 = dict()
		lint2 = dict()
		nodes_in_order2 = []
		for key in nodes_in_order:
			if key in Tint:
				ancestor = Tint[key]
				Tint2[nodes_to_index[key]] = nodes_to_index[ancestor]
				if (key, ancestor) in lint:
					lint2[nodes_to_index[key], nodes_to_index[ancestor]] = lint[key, ancestor]
			nodes_in_order2.append(nodes_to_index[key])

		# Next make the probability distributions
		# Would be nice if I could find a non-destructive way to subtract up and normalize

		# Do it for P
		self._subtract_down()
		keys = _data.keys()
		total_abundance = 0
		for key in keys:
			total_abundance += _data[key]["abundance"]
		# print(total_abundance)
		for key in keys:
			if total_abundance > 0:
				_data[key]["abundance"] /= total_abundance  # Should be a fraction, summing to 1
		P = np.zeros(len(nodes_in_order))
		for key_ind in xrange(len(nodes_in_order)):
			key = nodes_in_order[key_ind]
			if key in _data:
				P[key_ind] = _data[key]["abundance"]

		# Make back into percentages and add the mass back up (effectively normalizing the vector)
		for key in keys:
			if total_abundance > 0:
				_data[key]["abundance"] *= 100
		self._add_up()

		# Next do for Q
		other._subtract_down()
		keys = _other_data.keys()
		total_abundance = 0
		for key in keys:
			total_abundance += _other_data[key]["abundance"]
		# print(total_abundance)
		for key in keys:
			if total_abundance > 0:
				_other_data[key]["abundance"] /= total_abundance  # should be a fraction, summing to 1
		Q = np.zeros(len(nodes_in_order))
		for key_ind in xrange(len(nodes_in_order)):
			key = nodes_in_order[key_ind]
			if key in _other_data:
				Q[key_ind] = _other_data[key]["abundance"]

		# Make back into percentages and add the mass back up (effectively normalizing the vector)
		for key in keys:
			if total_abundance > 0:
				_other_data[key]["abundance"] *= 100
		other._add_up()

		return Tint2, lint2, nodes_in_order2, nodes_to_index, P, Q


def test_normalize():
	profile = Profile('/home/dkoslicki/Dropbox/Repositories/CAMIProfilingTools/src/test1.profile')
	profile.write_file('/home/dkoslicki/Dropbox/Repositories/CAMIProfilingTools/src/test1.profile.import')
	profile.normalize()
	profile.write_file('/home/dkoslicki/Dropbox/Repositories/CAMIProfilingTools/src/test1.profile.normalize')
	return profile


def test_unifrac():
	sys.path.append('/home/dkoslicki/Dropbox/Repositories/EMDUnifrac/src')
	import EMDUnifrac as EMDU
	profile1 = Profile('/home/dkoslicki/Dropbox/Repositories/EMDUnifrac/data/test1.profile')
	profile2 = Profile('/home/dkoslicki/Dropbox/Repositories/EMDUnifrac/data/test2.profile')
	t0 = timeit.default_timer()
	(Tint, lint, nodes_in_order, nodes_to_index, P, Q) = profile1.make_unifrac_input_and_normalize(profile2)
	t1 = timeit.default_timer()
	#print(t1-t0)
	(Z, diffab) = EMDU.EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q)
	print(Z)
	(Tint, lint, nodes_in_order, nodes_to_index, P, Q) = profile1.make_unifrac_input_and_normalize(profile2)
	(Z2, diffab2) = EMDU.EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q)
	print(Z2)
	profile1.normalize()
	profile2.normalize()
	(Tint, lint, nodes_in_order, nodes_to_index, P, Q) = profile1.make_unifrac_input_and_normalize(profile2)
	(Z3, diffab3) = EMDU.EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q)
	print(Z3)

	print(diffab)
	print(diffab2)
	print(diffab3)


def test_real_data():
	profile1 = Profile('/home/dkoslicki/Dropbox/Repositories/CAMIProfilingTools/src/lane4-s041-indexN722-S502-ATGCGCAG-ATAGAGAG-41_M5-2_S41_L004_R1_001.fa.gz.metaphlan.profile')
	profile2 = Profile('/home/dkoslicki/Dropbox/Repositories/CAMIProfilingTools/src/lane8-s092-indexN729-S505-TCGACGTC-CTCCTTAC-91_Z0299_S92_L008_R1_001.fa.gz.metaphlan.profile')
	sys.path.append('/home/dkoslicki/Dropbox/Repositories/EMDUnifrac/src')
	import EMDUnifrac as EMDU
	(Tint, lint, nodes_in_order, nodes_to_index, P, Q) = profile1.make_unifrac_input_and_normalize(profile2)
	(Z, diffab) = EMDU.EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q)
	print(Z)
	(Tint, lint, nodes_in_order, nodes_to_index, P, Q) = profile1.make_unifrac_input_and_normalize(profile2)
	(Z2, diffab2) = EMDU.EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q)
	print(Z2)
	profile1.normalize()
	profile2.normalize()
	(Tint, lint, nodes_in_order, nodes_to_index, P, Q) = profile1.make_unifrac_input_and_normalize(profile2)
	(Z3, diffab3) = EMDU.EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q)
	print(Z3)

	print(diffab)
	print(diffab2)
	print(diffab3)
	#return profile1, profile2

