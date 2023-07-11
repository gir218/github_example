from fchk import Fchk
from sys import argv



def main(input_file,output_file):
	infile=Fchk(input_file,output_file)
	infile.to_pickle(output_file=output_file)
	

if __name__ == "__main__":
    main(argv[1],argv[2])
