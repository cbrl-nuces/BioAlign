import sys

if __name__ == '__main__':
	sp1 = sys.argv[1]
	sp2 = sys.argv[2]
	stage = sys.argv[3]
	for i in range(0,2):
		sum_one = 0.0
		avg_one = 0.0
		node_one = 0.0
		
		if i==0:
		    f = open("semantic_files/"+sp1+"-"+sp2+"-"+stage+"-MF.txt")
		if i==1:
		    f = open("semantic_files/"+sp1+"-"+sp2+"-"+stage+"-BP.txt")

		mf_lines = f.readlines()
		f.close()
		count=0
		for j in range(len(mf_lines)):	

		    l = mf_lines[j].strip("\n")
		    l = l.split("\t")
		    if l[-1] != 'NA':
		        sum_one+=float(l[-1])
		        node_one+=1
		        count+=1

		avg_one = sum_one/float(node_one)
		if i==0 and sp1 == 'yeast':
		    print ("MF Results for",sp1,"-",sp2,"pair")
		    print ("MF: ",round(avg_one,3), "Aligned Nodes: ", round(node_one*100/5036,2)) # max length in case of yeast-human
		elif i==1 and sp1 == 'yeast':
		    print ("BP Results for",sp1,"-",sp2,"pair")
		    print ("BP: ",round(avg_one,3), "Aligned Nodes: ", round(node_one*100/5036,2))

		if i==0 and sp1 == 'mouse':
		    print ("MF Results for",sp1,"-",sp2,"pair")
		    print ("MF: ",round(avg_one,3), "Aligned Nodes: ", round(node_one*100/744,2)) # max length in case of mouse cases
		elif i==1 and sp1 == 'mouse':
		    print ("BP Results for",sp1,"-",sp2,"pair")
		    print ("BP: ",round(avg_one,3), "Aligned Nodes: ", round(node_one*100/744,2))

