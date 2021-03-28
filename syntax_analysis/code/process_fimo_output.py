
import sys

def does_exist(cur_list, start_loc, end_loc):
	for s, e, _ in cur_list:
		if (s >= start_loc and s <= end_loc) or (start_loc >= s and start_loc <= e):
			return True
	return False

def read_fimo(fimo_input):
	fimo_dict = {}
	with open(fimo_input, "r") as f:
		for line in f:
			line = line.strip()
			if line == "":
				continue
			if line[0] == '#':
				continue
			vals = line.split('\t')
			#['Pitx2', '', 'chr2:119570822-119571033', '1', '8', '+', '11.5165', '1.79e-05', '0.671', 'ttaatccc']
			if vals[2] not in fimo_dict:
				fimo_dict[vals[2]] = []
			if not does_exist(fimo_dict[vals[2]], vals[3], vals[4]):
				fimo_dict[vals[2]].append((vals[3], vals[4], vals[6]))
	return fimo_dict

def write_output_per_site(fimo_dict, output_file):
	with open(output_file, "w") as f:
		for peak in fimo_dict:
			[chromosome, locs] = peak.split(':')
			[start_loc, end_loc] = locs.split('-')
			for (motif_start, motif_end, score) in fimo_dict[peak]:
				f.write("%s\t%s\t%s\t%s,%s,%s\n" % (chromosome, start_loc, end_loc, \
				motif_start, motif_end, score))

	return 0

def write_output_per_peak(fimo_dict, output_file):
	with open(output_file, "w") as f:
		for peak in fimo_dict:
			[chromosome, locs] = peak.split(':')
			[start_loc, end_loc] = locs.split('-')
			f.write("%s\t%s\t%s\t%d\n" % (chromosome, start_loc, end_loc, len(fimo_dict[peak])))

	return 0

def main(argv):
	if argv is None:
		argv = sys.argv
	else:
		argv = argv.split()
	
	fimo_input = argv[1]
	fimo_dict = read_fimo(fimo_input)
	output_file = argv[2]
	write_output_per_peak(fimo_dict, output_file)
	return 0

if __name__ == "__main__":
	sys.exit(main(None))
