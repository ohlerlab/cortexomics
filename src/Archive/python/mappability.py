import sys
import pdb
kmerdict = {}
sequence = ''
kmerlength = 4	
i=0
kmersdone=0

if __name__ == "__main__":
	for line in sys.stdin:

		# sys.stderr.write("DEBUG: got line: " + line)
		if(line[0] is '>'): 
			sys.stderr.write('skipping '+line)
			sequence = ''
			continue
		if(line[0] not in 'ACGTacgt'): 
			sys.stderr.write('skipping '+line)
			continue
		line = line.rstrip('\n')			
		
		sequence = sequence + line
		
		sys.stderr.write('processing '+sequence)
		while i <= (len(sequence) - kmerlength):
			kmerdict[(sequence[i:i+kmerlength])] = True
			kmerdict[hash(sequence[i:i+kmerlength])] = True
			kmersdone+=1
			i+=1

		print(len(sequence))

		print('seq ' + sequence)
		print('left ' + sequence[i:])
		sequence = sequence[i:]
		print(sequence)
		# print( kmerdict)


#        $i=0;while($i<=(length($F[3])-10)){print $F[0] ,"_",$F[1]+$i," ",substr $F[3],$i,10 ; $i = $i +1}'#


	print(sys.getsizeof(kmerdict))
	print(kmersdone)
		