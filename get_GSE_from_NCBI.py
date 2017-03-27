import argparse

parser = argparse.ArgumentParser(description = 'Description')
parser.add_argument('-i', '--input', help='Name of input file', required=True, type=argparse.FileType('r'), metavar='FILE')
parser.add_argument('-o', '--output', help='Name of output file', required=True, type=argparse.FileType('w'), metavar='FILE')
args = parser.parse_args()

#file_out = open(args.output)

#with open(args.input) as read_file:
for line in args.input:
    if line.startswith("Series"):
        st = line.split()
        gse = str(st[2]) + str("\n")
        args.output.write(gse)
