import re
import sys
import array as arr
a = [32, 50, 64, 100, 128, 200, 256, 300, 400, 512, 1000, 1024, 2000, 2048, 3000, 4000, 4096, 5000, 10000]
i = 0
with open(sys.argv[1]) as file:
	for linha in file:
		if "L2 miss ratio" in linha:
			linha = re.findall("\d+\.\d+",linha)
			if len(linha) == 0:
				print("{} {}".format(a[i], 0.0))
			else:
				print("{} {}".format(a[i],float(linha[0])))
			i = i + 1		
		if "L3 bandwidth [MBytes/s]" in linha:
			#if arquivo in sys.argv[1]:
			linha = re.findall("\d+\.\d+",linha)
			if len(linha) == 0: 
				print("{} {}".format(a[i], 0.0))
			else:
				print("{} {}".format(a[i],float(linha[0])))
			i = i + 1
		if "DP MFLOP/s" in linha:
			if sys.argv[2] == 'DP':
				if not ("AVX" in linha):
					linha = re.findall("\d+\.\d+",linha)
					if len(linha) == 0:
						print("{} {}".format(a[i], 0.0))
					else:
						print("{} {}".format(a[i],float(linha[0])))
					i = i + 1
		if "AVX DP MFLOP/s" in linha:
			if sys.argv[2] == 'AVX':
				linha = re.findall("\d+\.\d+",linha)
				if len(linha) == 0:
					print("{} {}".format(a[i], 0.0))
				else:
					print("{} {}".format(a[i],float(linha[0])))		
				i = i + 1  
	
